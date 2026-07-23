#!/usr/bin/env bash

set -eu

usage() {
    echo "Usage: $0 REFERENCE_STDOUT CANDIDATE_STDOUT [DIFF_OUTPUT]" >&2
}

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    usage
    exit 2
fi

reference=$1
candidate=$2
diff_output=${3:-}

if cmp -s "$reference" "$candidate"; then
    echo "exact"
    exit 0
fi

normalized=$(mktemp "${TMPDIR:-/tmp}/telseq-output.XXXXXX")
trap 'rm -f "$normalized"' EXIT HUP INT TERM

# A current telseq-parallel result is compatible with a stock/reference result
# when removing exactly one final K column makes the complete stdout identical.
# Validate the header (when present) and every result-row value while doing
# that normalization. Headerless output produced with -H is also supported.
if awk '
    BEGIN {
        found_header = 0
        found_row = 0
    }

    /^ReadGroup\tLibrary\tSample\tTotal\tMapped\tDuplicates\tLENGTH_ESTIMATE\t/ {
        if ($0 !~ /\tK\t$/) {
            exit 10
        }
        sub(/\tK\t$/, "\t")
        found_header = 1
        print
        next
    }

    /\t$/ {
        if ($0 !~ /\t[0-9]+\t$/) {
            exit 11
        }
        sub(/\t[0-9]+\t$/, "\t")
        found_row = 1
        print
        next
    }

    {
        print
    }

    END {
        if (!found_header && !found_row) {
            exit 12
        }
    }
' "$candidate" > "$normalized" &&
   cmp -s "$reference" "$normalized"; then
    echo "appended-k"
    exit 0
fi

if [ -n "$diff_output" ]; then
    if [ -s "$normalized" ]; then
        diff -u "$reference" "$normalized" > "$diff_output" || true
    else
        diff -u "$reference" "$candidate" > "$diff_output" || true
    fi
fi

exit 1
