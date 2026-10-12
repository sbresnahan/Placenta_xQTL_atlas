#!/usr/bin/env bash
set -euo pipefail

BOOK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$BOOK_DIR/.." && pwd)"
BRANCH="${1:-gh-pages}"
REMOTE="${REMOTE:-origin}"
R_SCRIPT="${R_SCRIPT:-Rscript}"
SITE_DIR="$BOOK_DIR/_book"

cd "$REPO_ROOT"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "ERROR: $REPO_ROOT is not a git working tree" >&2
  exit 1
fi

"$R_SCRIPT" "$BOOK_DIR/render_book.R"

if [[ ! -f "$SITE_DIR/index.html" ]]; then
  echo "ERROR: book build did not produce $SITE_DIR/index.html" >&2
  exit 1
fi

WORKTREE="$(mktemp -d "${TMPDIR:-/tmp}/placenta-xqtl-pages.XXXXXX")"
CNAME_SAVE=""
cleanup() {
  set +e
  cd "$REPO_ROOT"
  git worktree remove --force "$WORKTREE" >/dev/null 2>&1 || true
  rm -rf "$WORKTREE"
  if [[ -n "$CNAME_SAVE" ]]; then rm -f "$CNAME_SAVE"; fi
}
trap cleanup EXIT

if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git worktree add "$WORKTREE" "$BRANCH"
elif git ls-remote --exit-code --heads "$REMOTE" "$BRANCH" >/dev/null 2>&1; then
  git fetch "$REMOTE" "$BRANCH"
  git worktree add -b "$BRANCH" "$WORKTREE" "$REMOTE/$BRANCH"
else
  git worktree add --detach "$WORKTREE" HEAD
  (
    cd "$WORKTREE"
    git switch --orphan "$BRANCH"
    git rm -rf . >/dev/null 2>&1 || true
  )
fi

if [[ -f "$WORKTREE/CNAME" ]]; then
  CNAME_SAVE="$(mktemp "${TMPDIR:-/tmp}/placenta-xqtl-cname.XXXXXX")"
  cp "$WORKTREE/CNAME" "$CNAME_SAVE"
fi

# Replace the publication branch with the static book. The worktree's .git file
# is retained. CNAME is preserved if the branch already uses a custom domain.
find "$WORKTREE" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
cp -a "$SITE_DIR"/. "$WORKTREE"/
[[ -n "$CNAME_SAVE" ]] && cp "$CNAME_SAVE" "$WORKTREE/CNAME"
touch "$WORKTREE/.nojekyll"

(
  cd "$WORKTREE"
  git add -A
  if git diff --cached --quiet; then
    echo "GitHub Pages branch is already up to date."
    exit 0
  fi
  git commit -m "Update Placenta xQTL Atlas bookdown site"
  if [[ "${PUBLISH_NO_PUSH:-0}" == "1" ]]; then
    echo "Committed locally to $BRANCH; PUBLISH_NO_PUSH=1, so no push was performed."
  else
    git push -u "$REMOTE" "$BRANCH"
    echo "Published static site to $REMOTE/$BRANCH."
  fi
)
