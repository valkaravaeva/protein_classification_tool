#!/usr/bin/env bash
#
# Switch which KEGG module set the pipeline runs against.
#
# pipeline.nf always reads data/kegg_modules_architecture.txt. This script copies the
# module set you want over that file, so you don't have to do the copy/rename by hand:
#
#   scripts/select_kegg_modules.sh manual     # data/manual_modules_architecture.txt
#                                              # -> data/kegg_modules_architecture.txt
#                                              # (adds the manually curated modules
#                                              # described in the publication - oxygen/
#                                              # sulfur/nitrogen reductases, archaeal
#                                              # riboflavin/heme/cobalamin variants, ...)
#
#   scripts/select_kegg_modules.sh original   # data/kegg_modules_architecture_original.txt
#                                              # -> data/kegg_modules_architecture.txt
#                                              # (unmodified KEGG modules only)
#
# Run this before `nextflow run pipeline.nf`.

set -euo pipefail

usage() {
    echo "Usage: $0 {manual|original}" >&2
    echo >&2
    echo "  manual    use the manually curated modules (as described in the publication)" >&2
    echo "  original  use the unmodified KEGG modules" >&2
    exit 1
}

if [[ $# -ne 1 ]]; then
    usage
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$script_dir/../data"
target="$data_dir/kegg_modules_architecture.txt"

case "$1" in
    manual)
        source_file="$data_dir/manual_modules_architecture.txt"
        label="manually curated modules (manual_modules_architecture.txt)"
        ;;
    original)
        source_file="$data_dir/kegg_modules_architecture_original.txt"
        label="original, unmodified KEGG modules (kegg_modules_architecture_original.txt)"
        ;;
    *)
        usage
        ;;
esac

if [[ ! -f "$source_file" ]]; then
    echo "Error: $source_file not found" >&2
    exit 1
fi

cp "$source_file" "$target"
echo "kegg_modules_architecture.txt now uses the $label"
