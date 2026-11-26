#!/usr/bin/env bash
# Rebuild the anonymized JBES replication archive (jbes_replication/ and jbes_replication.zip).
# Run from anywhere:  bash "2025 sdf paper/jbes/build_jbes_replication.sh"
#
# 1. Every file in jbes_replication/ that also exists under r_project/ (same relative path)
#    is refreshed from r_project/, so the archive cannot drift from the code that produced
#    the paper. Package-only files (README.md, run_synthetic_verification.R, output/) are kept.
#    To add a file to the archive, copy it into jbes_replication/ once; later builds refresh it.
# 2. The archive is scanned for identifying strings; the build stops on a hit.
# 3. jbes_replication.zip is rebuilt (no macOS metadata).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
pkg="$here/jbes_replication"
src="$(cd "$here/../../r_project" && pwd)"

refreshed=0
while IFS= read -r -d '' f; do
  rel="${f#"$pkg"/}"
  case "$rel" in README.md|run_synthetic_verification.R|output/*|.here) continue ;; esac
  if [ -f "$src/$rel" ] && ! cmp -s "$src/$rel" "$f"; then
    cp -p "$src/$rel" "$f"; echo "refreshed: $rel"; refreshed=$((refreshed + 1))
  fi
done < <(find "$pkg" -type f ! -name '.DS_Store' -print0)
echo "files refreshed from r_project/: $refreshed"

# Identifying strings: author names, affiliations, e-mail domains, local paths, prior review
pattern='tausch|bohnert|assetmetrix|quant-unit|hm\.edu|web\.de|gmail|christausch|/Users/|icloud|fundwise|jfec|jfinecon|referee'
if grep -rIniE "$pattern" "$pkg"; then
  echo "ERROR: identifying strings found in $pkg (see above)" >&2; exit 1
fi
echo "leak scan: clean"

cd "$here"
rm -f jbes_replication.zip
zip -qrX jbes_replication.zip jbes_replication -x '*.DS_Store'
unzip -l jbes_replication.zip | tail -1
