#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# remove_cursor_coauthor.sh
# Removes all "Co-authored-by" lines referencing Cursor from every commit
# in your repository's history, then force-pushes the cleaned history.
#
# Usage:
#   cd /path/to/your/dt-fusion-illicit
#   bash remove_cursor_coauthor.sh
#
# WARNING: This rewrites git history. Coordinate with collaborators —
#          they will need to re-clone or reset their local copies afterward.
# =============================================================================

# Safety check: make sure we're in a git repo
if ! git rev-parse --is-inside-work-tree &>/dev/null; then
    echo "ERROR: Not inside a git repository. cd into your repo first."
    exit 1
fi

BRANCH=$(git rev-parse --abbrev-ref HEAD)
echo "Current branch: $BRANCH"
echo ""

# Show how many commits currently have Cursor co-author trailers
COUNT=$(git log --all --grep="Co-authored-by.*[Cc]ursor" --oneline | wc -l)
echo "Found $COUNT commit(s) with Cursor co-author trailers."
echo ""

if [ "$COUNT" -eq 0 ]; then
    echo "Nothing to do! No Cursor co-author trailers found."
    exit 0
fi

echo "This will rewrite history for ALL branches and tags."
read -p "Continue? (y/N) " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "Aborted."
    exit 0
fi

echo ""
echo "Rewriting commit messages..."

# Use git filter-branch to strip Co-authored-by lines mentioning Cursor
git filter-branch -f --msg-filter '
    sed -E "/^Co-authored-by:.*[Cc]ursor.*$/d"
' --tag-name-filter cat -- --all

echo ""
echo "Done rewriting history."
echo ""

# Clean up the backup refs that filter-branch creates
echo "Cleaning up backup refs..."
git for-each-ref --format="%(refname)" refs/original/ | while read ref; do
    git update-ref -d "$ref"
done
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo ""
echo "============================================"
echo " History rewritten successfully!"
echo "============================================"
echo ""
echo "Next steps:"
echo "  1. Verify the fix:  git log --all --grep='Co-authored-by.*Cursor' --oneline"
echo "     (should return nothing)"
echo ""
echo "  2. Force-push ALL branches:"
echo "     git push origin --force --all"
echo "     git push origin --force --tags"
echo ""
echo "  3. Tell collaborators to re-clone the repo (or run: git fetch --all && git reset --hard origin/$BRANCH)"
echo ""
echo "  4. If GitHub still shows Cursor as a contributor after a few hours,"
echo "     you may need to contact GitHub Support — it can be a caching issue."
echo ""
echo "  5. (Optional) Delete the commit_fix.sh and path/to/the/commit/ from your repo"
echo "     if those were from a previous attempt to fix this."