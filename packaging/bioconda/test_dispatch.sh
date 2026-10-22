#!/bin/bash
# Recipe test: `coal` finds the installed coal-cpp without any environment variable and hands large inputs to it.
#   1. COAL_ENGINE=cpp with no COAL_CPP set: coal must find coal-cpp next to its own interpreter or on PATH
#      (dispatch.find_cpp), otherwise it exits 1 with "coal-cpp was not found".
#   2. Automatic mode on a >= 1 MB input (dispatch.CPP_MIN_BYTES) must run coal-cpp. Which engine ran is invisible in
#      the output (that is the point), so COAL_CPP points at a one-line wrapper that leaves a marker and execs the
#      real, installed coal-cpp; the output must equal the Python engine's.
set -euo pipefail
unset COAL_CPP COAL_ENGINE
here=$(pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

real=$(command -v coal-cpp)
test -n "$real"

# 1. Discovery without COAL_CPP.
COAL_ENGINE=cpp coal test_protein.aln.fasta test_dna.fasta > "$work/forced.fasta"

# 2. A >= 1 MB input: the 8-gene fixture repeated with unique IDs.
n=240
for i in $(seq 1 $n); do
  sed "s/^>\(.*\)/>\1_$i/" test_protein.aln.fasta
done > "$work/big_pep.fasta"
for i in $(seq 1 $n); do
  sed "s/^>\(.*\)/>\1_$i/" test_dna.fasta
done > "$work/big_dna.fasta"
size=$(( $(wc -c < "$work/big_pep.fasta") + $(wc -c < "$work/big_dna.fasta") ))
test "$size" -ge 1000000

cat > "$work/coal-cpp" <<EOF
#!/bin/sh
case "\$1" in --version) ;; *) touch "$work/cpp_ran" ;; esac
exec "$real" "\$@"
EOF
chmod +x "$work/coal-cpp"

COAL_CPP="$work/coal-cpp" coal "$work/big_pep.fasta" "$work/big_dna.fasta" > "$work/auto.fasta"
test -e "$work/cpp_ran" || { echo "coal did not hand a $size-byte input to coal-cpp" >&2; exit 1; }
rm "$work/cpp_ran"
COAL_ENGINE=python coal "$work/big_pep.fasta" "$work/big_dna.fasta" > "$work/python.fasta"
test ! -e "$work/cpp_ran"
cmp "$work/auto.fasta" "$work/python.fasta"
test "$(grep -c '^>' "$work/auto.fasta")" -eq $((8 * n))
echo "dispatch OK: $size-byte input went to $real; output identical to the Python engine"
