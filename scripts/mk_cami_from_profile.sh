#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -m mapping.tsv[.gz] -i profile.tsv -o out.profile

Inputs:
  -m  Mapping: GTDB species -> NCBI species taxid  (2 columns)
      e.g. "0-14-0-20-34-12_sp038046645<TAB>1906665"
      (gzipped or plain)

  -i  Abundance profile: sample  internal+gtdbSpecies  abundance_fraction
      e.g. "SRR...<TAB>5138_Petralouisia_sp009774455<TAB>0.03488"

  -o  Output CAMI profile (all samples concatenated)

Notes:
  * Abundances are normalized PER SAMPLE and converted to PERCENT for TaxonKit (-p).
  * Requires: taxonkit (in PATH), awk, zcat (or gzip -dc).
EOF
  exit 1
}

MAP=""
PROF=""
OUT=""

while getopts ":m:i:o:" opt; do
  case "$opt" in
    m) MAP="$OPTARG" ;;
    i) PROF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$MAP" || -z "$PROF" || -z "$OUT" ]] && usage
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not found in PATH" >&2; exit 2; }

# Clean output
: > "$OUT"

# Pick reader for mapping (gz or not)
if [[ "$MAP" =~ \.gz$ ]]; then
  MAPREAD="zcat --quiet --force"
else
  MAPREAD="cat"
fi

# Loop all samples in the profile file
cut -f1 "$PROF" | LC_ALL=C sort -u | while read -r S; do
  # Build TaxId\tPercent table on the fly for this sample and feed to taxonkit
  awk -F'\t' -v S="$S" '
    BEGIN{ OFS="\t" }
    # ---------- load mapping (first input stream via process substitution) ----------
    NR==FNR {
      core=$1; sub(/^[^_]+_/, "", core);         # drop leading internal token if any
      taxid=$2
      map_core[core]=taxid
      if (match(core, /_sp[0-9]+$/)) {
        suf=substr(core, RSTART+1)               # e.g. sp038046645
        map_sp[suf]=taxid
      }
      next
    }
    # ---------- process profile rows for this sample ----------
    $1==S {
      sp=$2; sub(/^[^_]+_/, "", sp)              # strip internal token from species
      suf=""
      if (match(sp, /_sp[0-9]+$/)) suf=substr(sp, RSTART+1)

      taxid = (sp in map_core) ? map_core[sp] : ((suf in map_sp) ? map_sp[suf] : "")
      if (taxid=="") next                        # skip unmapped

      abd = $3 + 0                               # fraction (0..1)
      sum += abd
      agg[taxid] += abd
    }
    END {
      for (t in agg) {
        # Normalize to percent for taxonkit -p
        printf "%s\t%.6f\n", t, 100.0*agg[t]/sum
      }
    }
  ' <($MAPREAD "$MAP") "$PROF" \
  | taxonkit profile2cami -p -s "$S" >> "$OUT"
done

echo "Done: wrote $(wc -l < "$OUT") lines to $OUT"

