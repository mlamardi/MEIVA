#!/usr/bin/env bash
#
# run_meiva_cohort.sh — annotate a large MEI cohort with MEIVA, chromosome by
# chromosome, in parallel.
#
# WHY THIS EXISTS
#
# `meiva annotate` holds every parsed call in memory while it merges them into
# loci. That is fine for tens of samples and a problem for thousands: a parsed
# xTEA call costs roughly 5.4 kB, so a cohort of 2,500 genomes at ~1,500 calls
# each needs about 20 GB. Splitting the work by chromosome brings the peak down
# to under 2 GB per worker and lets several run at once.
#
# The split is EXACTLY equivalent to a single whole-genome run, not an
# approximation. MEIVA only ever clusters calls that share a contig and an
# element family, so no locus can span two chromosomes and no locus boundary can
# move. The only thing that changes is peak memory.
#
# The script is resumable: every stage skips work already done, so it is safe to
# re-run after an interruption.
#
# USAGE
#
#   ./run_meiva_cohort.sh --vcf-dir /path/to/vcfs --gencode gencode.v47.gtf.gz \
#       --outdir results
#
# Run with --help for all options.
#
set -euo pipefail

# ------------------------------------------------------------------ defaults --
VCF_DIR=""
VCF_LIST=""
GENCODE=""
OUTDIR=""
WORKERS=4
KEEP_SPLIT=0

# Autosomes by default. chrX and chrY are left out because ploidy differs by sex,
# so a single carrier frequency across a mixed-sex cohort is not interpretable
# for them; they belong in a separate sex-stratified analysis. Add them with
#   --chroms "chr1 ... chr22 chrX chrY"
# If your reference has no "chr" prefix, pass --chroms "1 2 3 ... 22".
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

F6_DEGS=""; F6_SAMPLES=""; F6_CAT=""
F5_ENHANCERS=""; F5_PEAKS=""; F5_PEAK_NAMES=""

usage() {
  cat << 'EOF'
Annotate a large MEI cohort with MEIVA, chromosome by chromosome, in parallel.

Required:
  --vcf-dir DIR            directory of per-sample caller VCFs (*.vcf)
  --vcf-list FILE          ...or a file listing VCF paths, one per line
  --gencode GTF            GENCODE annotation, .gtf or .gtf.gz
  --outdir DIR             where to write results and intermediates

Optional:
  --workers N              parallel chromosomes (default 4)
                           Budget ~2 GB per worker, plus ~0.5 GB for the gene
                           model and ~0.1 GB if FANTOM5 is enabled.
  --chroms "chr1 chr2 ..." contigs to process (default: autosomes)
  --keep-split             keep the per-chromosome VCFs (they are large)
  -h, --help               show this message

FANTOM6 lncRNA functional evidence (all three together, or none):
  --fantom6-degs FILE      DESeq2_genes_ASO_signif.tsv[.bz2]
  --fantom6-samples FILE   Published_sample_summary.tsv[.bz2]
  --fantom6-cat FILE       FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv[.gz]

FANTOM5 regulatory context (enhancers and promoters are independent):
  --fantom5-enhancers FILE F5.hg38.enhancers.bed[.gz]
  --fantom5-peaks FILE     hg38_fair+new_CAGE_peaks_phase1and2.bed[.gz]
  --fantom5-peak-names F   human_phase1and2_CAGE_Peak_name.txt[.gz]

Example:
  ./run_meiva_cohort.sh \
      --vcf-dir ~/cohort/vcfs \
      --gencode ~/refs/gencode.v47.annotation.gtf.gz \
      --outdir ~/cohort/meiva \
      --workers 4 \
      --fantom5-enhancers ~/refs/F5.hg38.enhancers.bed.gz
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --vcf-dir) VCF_DIR="$2"; shift 2 ;;
    --vcf-list) VCF_LIST="$2"; shift 2 ;;
    --gencode) GENCODE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --chroms) CHROMS="$2"; shift 2 ;;
    --keep-split) KEEP_SPLIT=1; shift ;;
    --fantom6-degs) F6_DEGS="$2"; shift 2 ;;
    --fantom6-samples) F6_SAMPLES="$2"; shift 2 ;;
    --fantom6-cat) F6_CAT="$2"; shift 2 ;;
    --fantom5-enhancers) F5_ENHANCERS="$2"; shift 2 ;;
    --fantom5-peaks) F5_PEAKS="$2"; shift 2 ;;
    --fantom5-peak-names) F5_PEAK_NAMES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; echo; usage; exit 2 ;;
  esac
done

# -------------------------------------------------------------- validation --
die() { echo "error: $*" >&2; exit 1; }

[ -n "$GENCODE" ] || { usage; die "--gencode is required"; }
[ -n "$OUTDIR" ] || die "--outdir is required"
[ -f "$GENCODE" ] || die "GENCODE file not found: $GENCODE"
if [ -z "$VCF_DIR" ] && [ -z "$VCF_LIST" ]; then die "one of --vcf-dir or --vcf-list is required"; fi
if [ -n "$VCF_DIR" ] && [ -n "$VCF_LIST" ]; then die "--vcf-dir and --vcf-list are mutually exclusive"; fi
command -v meiva > /dev/null || die "the 'meiva' command is not on PATH; pip install meiva"

for f in "$F6_DEGS" "$F6_SAMPLES" "$F6_CAT" "$F5_ENHANCERS" "$F5_PEAKS" "$F5_PEAK_NAMES"; do
  if [ -n "$f" ] && [ ! -f "$f" ]; then die "file not found: $f"; fi
done

F6_ARGS=""
n_f6=0
for f in "$F6_DEGS" "$F6_SAMPLES" "$F6_CAT"; do [ -n "$f" ] && n_f6=$((n_f6 + 1)); done
if [ "$n_f6" = "3" ]; then
  F6_ARGS="--fantom6-degs $F6_DEGS --fantom6-samples $F6_SAMPLES --fantom6-cat $F6_CAT"
elif [ "$n_f6" != "0" ]; then
  die "the three --fantom6-* options must be given together"
fi

F5_ARGS=""
if [ -n "$F5_ENHANCERS" ]; then F5_ARGS="$F5_ARGS --fantom5-enhancers $F5_ENHANCERS"; fi
if [ -n "$F5_PEAKS" ]; then F5_ARGS="$F5_ARGS --fantom5-peaks $F5_PEAKS"; fi
if [ -n "$F5_PEAK_NAMES" ]; then
  [ -n "$F5_PEAKS" ] || die "--fantom5-peak-names requires --fantom5-peaks"
  F5_ARGS="$F5_ARGS --fantom5-peak-names $F5_PEAK_NAMES"
fi

SPLIT="$OUTDIR/split"
GTFDIR="$OUTDIR/gtf"
OUT="$OUTDIR/per_chrom"
FINAL="$OUTDIR/meiva_cohort.tsv"
mkdir -p "$OUTDIR" "$SPLIT" "$GTFDIR" "$OUT"

# gzcat on macOS, zcat on Linux
if command -v gzcat > /dev/null; then ZCAT=gzcat; else ZCAT=zcat; fi
case "$GENCODE" in *.gz) GTF_READ="$ZCAT \"$GENCODE\"" ;; *) GTF_READ="cat \"$GENCODE\"" ;; esac

echo "=============================================================="
echo " MEIVA cohort run"
echo " output   : $OUTDIR"
echo " workers  : $WORKERS"
echo " FANTOM6  : $([ -n "$F6_ARGS" ] && echo enabled || echo skipped)"
echo " FANTOM5  : $([ -n "$F5_ARGS" ] && echo enabled || echo skipped)"
echo "=============================================================="

# ---------------------------------------------------------- 1. sample list --
MANIFEST="$OUTDIR/samples.txt"
if [ ! -s "$MANIFEST" ]; then
  echo
  echo "[1/5] building sample list"
  if [ -n "$VCF_DIR" ]; then
    find "$VCF_DIR" -maxdepth 1 -name '*.vcf' | sort > "$MANIFEST"
  else
    grep -v '^[[:space:]]*$' "$VCF_LIST" > "$MANIFEST"
  fi
else
  echo
  echo "[1/5] sample list exists, reusing"
fi
N=$(wc -l < "$MANIFEST" | tr -d ' ')
[ "$N" -gt 0 ] || die "no VCFs found"
echo "      $N samples"
while read -r v; do [ -f "$v" ] || die "listed VCF not found: $v"; done < "$MANIFEST"

# ------------------------------------------------- 2. split VCFs by contig --
echo
echo "[2/5] splitting VCFs by chromosome"
if [ -f "$SPLIT/.done" ]; then
  echo "      already split, skipping"
else
  for c in $CHROMS; do mkdir -p "$SPLIT/$c"; done
  i=0
  while read -r vcf; do
    i=$((i + 1))
    [ $((i % 250)) -eq 0 ] && printf "      %d / %d\n" "$i" "$N"
    b=$(basename "$vcf" .vcf)
    awk -v out="$SPLIT" -v b="$b" -v chroms="$CHROMS" -v logf="$OUTDIR/skipped_contigs.log" '
      BEGIN { n = split(chroms, a, " "); for (j = 1; j <= n; j++) keep[a[j]] = 1 }
      /^#/  { hdr = hdr $0 "\n"; next }
      { c = $1
        if (!(c in keep)) { skip++; next }
        f = out "/" c "/" b ".vcf"
        if (!(f in seen)) { printf "%s", hdr > f; seen[f] = 1 }
        print > f }
      END { if (skip) printf "%s: %d calls on contigs outside --chroms\n", b, skip >> logf }
    ' "$vcf"
  done < "$MANIFEST"
  touch "$SPLIT/.done"
  echo "      done"
fi
if [ -s "$OUTDIR/skipped_contigs.log" ]; then
  echo "      note: calls outside --chroms were skipped, see skipped_contigs.log"
fi

# --------------------------------------------------- 3. split GTF by contig --
echo
echo "[3/5] pre-splitting the GENCODE GTF"
if [ -f "$GTFDIR/.done" ]; then
  echo "      already split, skipping"
else
  for c in $CHROMS; do
    if [ ! -s "$GTFDIR/$c.gtf" ]; then
      eval "$GTF_READ" | awk -v c="$c" '$1 == c || /^#/' > "$GTFDIR/$c.gtf"
    fi
  done
  touch "$GTFDIR/.done"
  echo "      done"
fi

# ----------------------------------------------------------- 4. annotate ----
echo
echo "[4/5] annotating with $WORKERS workers"

# macOS xargs -I has a 255-byte replacement buffer, too small for an inline
# sh -c block, so the per-chromosome work lives in its own helper.
cat > "$OUTDIR/annotate_one.sh" << 'HELPER'
#!/usr/bin/env bash
set -euo pipefail
c="$1"
if [ -s "$OUT/$c.tsv" ]; then echo "      $c already done"; exit 0; fi
n=$(find "$SPLIT/$c" -maxdepth 1 -name '*.vcf' | wc -l | tr -d ' ')
if [ "$n" = "0" ]; then echo "      $c no VCFs, skipping"; exit 0; fi
echo "      $c starting ($n samples)"
# cd in so the argument list holds bare filenames; thousands of absolute paths
# would approach the shell's ARG_MAX limit.
cd "$SPLIT/$c"
# shellcheck disable=SC2086
meiva annotate --vcf *.vcf --gencode "$GTFDIR/$c.gtf" \
    $F6_ARGS $F5_ARGS -o "$OUT/$c.tsv.partial" 2> "$OUT/$c.log"
mv "$OUT/$c.tsv.partial" "$OUT/$c.tsv"
echo "      $c done: $(( $(wc -l < "$OUT/$c.tsv") - 1 )) loci"
HELPER
chmod +x "$OUTDIR/annotate_one.sh"

export SPLIT GTFDIR OUT F6_ARGS F5_ARGS
printf '%s\n' $CHROMS | xargs -P "$WORKERS" -n 1 bash "$OUTDIR/annotate_one.sh"

# -------------------------------------------------------- 5. concatenate ----
echo
echo "[5/5] concatenating"
first=""
for c in $CHROMS; do
  if [ -s "$OUT/$c.tsv" ]; then first="$c"; break; fi
done
[ -n "$first" ] || die "no per-chromosome output was produced; check $OUT/*.log"

head -1 "$OUT/$first.tsv" > "$FINAL"
for c in $CHROMS; do
  if [ -s "$OUT/$c.tsv" ]; then tail -n +2 "$OUT/$c.tsv" >> "$FINAL"; fi
done

# MEIVA derives cohort_size from the samples it actually observed on a given
# chromosome, so a sample with no calls there silently leaves the denominator.
# Restore the true cohort size and recompute the frequency against it.
awk -F'\t' -v OFS='\t' -v n="$N" 'NR==1{print;next} {$10=n; $11=sprintf("%.6f",$9/n); print}' \
    "$FINAL" > "$FINAL.tmp" && mv "$FINAL.tmp" "$FINAL"

if [ "$KEEP_SPLIT" = "0" ]; then
  echo "      removing per-chromosome VCFs (--keep-split to retain)"
  rm -rf "$SPLIT"
fi

LOCI=$(( $(wc -l < "$FINAL") - 1 ))
echo
echo "=============================================================="
printf " samples  : %s\n" "$N"
printf " loci     : %s\n" "$LOCI"
printf " size     : %s\n" "$(du -h "$FINAL" | cut -f1)"
printf " output   : %s\n" "$FINAL"
echo "=============================================================="
echo
echo "cohort_size check (should be a single value, $N):"
cut -f10 "$FINAL" | tail -n +2 | sort -u | head
