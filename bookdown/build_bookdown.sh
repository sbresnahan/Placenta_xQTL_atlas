#!/usr/bin/env bash
set -euo pipefail

BOOK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${R_SCRIPT:-Rscript}"

exec "$R_SCRIPT" "$BOOK_DIR/render_book.R"
