# Helper for the get_*.sh scripts of the COAL tutorials: sourced, not run.
#
# fetch_cds ACCESSION ASSEMBLY_NAME LOCAL_NAME SHA256
#   downloads <ACCESSION>_<ASSEMBLY_NAME>_cds_from_genomic.fna.gz from the NCBI genomes FTP site (over HTTPS)
#   into $DEST/LOCAL_NAME and checks its SHA-256. A file that is already there with the right checksum is kept.
#   RefSeq assembly versions (GCF_....N) do not change once released, so the checksums are fixed.
set -euo pipefail

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

download() {  # URL FILE
    if command -v curl >/dev/null 2>&1; then curl -fsSL --retry 3 -o "$2" "$1"
    else wget -q -O "$2" "$1"; fi
}

fetch_cds() {
    local acc=$1 asm=$2 name=$3 want=$4
    local n=${acc#GCF_}; n=${n%%.*}
    local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${n:0:3}/${n:3:3}/${n:6:3}/${acc}_${asm}/${acc}_${asm}_cds_from_genomic.fna.gz"
    local out="$DEST/$name"
    if [ -f "$out" ] && [ "$(sha256_of "$out")" = "$want" ]; then
        echo "ok (already here)  $out"; return 0
    fi
    download "$url" "$out.part"
    local got; got=$(sha256_of "$out.part")
    if [ "$got" != "$want" ]; then
        rm -f "$out.part"
        echo "error: checksum mismatch for $url" >&2
        echo "       expected $want" >&2
        echo "       got      $got" >&2
        return 1
    fi
    mv "$out.part" "$out"
    echo "ok                 $out"
}
