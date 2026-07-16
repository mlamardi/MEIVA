"""Allow ``python -m meiva`` to invoke the CLI."""

import sys

from meiva.cli import main

if __name__ == "__main__":
    sys.exit(main())
