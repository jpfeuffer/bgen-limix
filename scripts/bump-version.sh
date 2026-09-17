#!/usr/bin/env bash
# Bumps bgen's version everywhere it's still tracked by hand, then commits
# and tags. bgen/version.h is *not* in this list on purpose -- it's
# generated from CMakeLists.txt's project(VERSION ...) at configure time
# (see cmake/version.h.in), so it can never drift.
#
# Usage: scripts/bump-version.sh X.Y.Z
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 X.Y.Z" >&2
  exit 1
fi

new_version="$1"
if [[ ! "$new_version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "error: '$new_version' is not a X.Y.Z version" >&2
  exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ -n "$(git status --porcelain)" ]]; then
  echo "error: working tree is not clean; commit or stash first" >&2
  exit 1
fi

sed_inplace() {
  # BSD sed (macOS) requires an explicit (possibly empty) backup suffix
  # after -i; GNU sed treats that argument as the pattern instead.
  if sed --version >/dev/null 2>&1; then
    sed -i "$@"
  else
    sed -i '' "$@"
  fi
}

sed_inplace -E "s/^version = \"[0-9]+\.[0-9]+\.[0-9]+\"/version = \"${new_version}\"/" pixi.toml
sed_inplace -E "s/^version = \"[0-9]+\.[0-9]+\.[0-9]+\"/version = \"${new_version}\"/" pyproject.toml
sed_inplace -E "s/^  VERSION [0-9]+\.[0-9]+\.[0-9]+$/  VERSION ${new_version}/" CMakeLists.txt

echo "Bumped to ${new_version} in:"
git diff --stat -- pixi.toml pyproject.toml CMakeLists.txt

git add pixi.toml pyproject.toml CMakeLists.txt
git commit -m "Bump version to ${new_version}"
git tag -a "v${new_version}" -m "v${new_version}"

echo
echo "Committed and tagged v${new_version} locally. Push with:"
echo "  git push origin main && git push origin v${new_version}"
