#!/bin/bash
# Build and test this recipe with conda-build against the current checkout (no release tag needed).
# Not part of the Bioconda submission (do not copy it to bioconda-recipes).
#
#   bash packaging/bioconda/build_local.sh [WORKDIR] [extra conda-build args...]
#
# It makes a tarball of the committed HEAD (`git archive`, so uncommitted edits are NOT built), copies the recipe to
# WORKDIR/recipe with `version` set to the tree's version and `url`/`sha256` pointing at that tarball, and runs
# `conda build -c conda-forge -c bioconda` on the copy. The committed meta.yaml is never modified.
# Needs `conda-build` (conda-build >= 24, libmamba solver) on PATH; `CONDA_BUILD_ARGS` adds arguments.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(git -C "$here" rev-parse --show-toplevel)
work=${1:-$(mktemp -d)}
shift || true
mkdir -p "$work"
work=$(cd "$work" && pwd)

version=$(sed -n 's/^__version__ = "\(.*\)"$/\1/p' "$repo/src/coal/__init__.py")
test -n "$version"
tarball="$work/COAL-$version.tar.gz"
git -C "$repo" archive --format=tar.gz --prefix="COAL-$version/" -o "$tarball" HEAD
sha=$(sha256sum "$tarball" | cut -d' ' -f1)

rm -rf "$work/recipe"
mkdir -p "$work/recipe"
cp "$here/meta.yaml" "$here/build.sh" "$here/test_dispatch.sh" "$here/test_protein.aln.fasta" "$here/test_dna.fasta" \
  "$work/recipe/"
# Only the two placeholder lines and the source url change.
sed -i \
  -e "s|^{% set version = \"[^\"]*\" %}.*|{% set version = \"$version\" %}|" \
  -e "s|^{% set sha256 = \"[^\"]*\" %}.*|{% set sha256 = \"$sha\" %}|" \
  -e "s|^  url: https://github.com/.*|  url: file://$tarball|" \
  "$work/recipe/meta.yaml"
grep -q "file://$tarball" "$work/recipe/meta.yaml"
grep -q "$sha" "$work/recipe/meta.yaml"

echo "recipe copy: $work/recipe (version $version, source $tarball, sha256 $sha)"
conda build -c conda-forge -c bioconda --override-channels --output-folder "$work/out" ${CONDA_BUILD_ARGS:-} "$@" \
  "$work/recipe"
