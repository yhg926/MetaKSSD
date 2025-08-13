#!/usr/bin/env bash
set -euo pipefail

usage(){ echo "Usage: $0 -m mapping.tsv[.gz] -i profile.tsv -o out.profile [-k]"; exit 1; }

KEEP_UNMAPPED=0
MAP=""; PROF=""; OUT=""
while getopts ":m:i:o:k" opt; do
  case "$opt" in
    m) MAP="$OPTARG" ;;
    i) PROF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    k) KEEP_UNMAPPED=1 ;;
    *) usage ;;
  esac
done
[[ -z "$MAP" || -z "$PROF" || -z "$OUT" ]] && usage
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not found" >&2; exit 2; }

: > "$OUT"

# reader for mapping
if [[ "$MAP" =~ \.gz$ ]]; then MAPREAD="zcat --quiet --force"; else MAPREAD="cat"; fi

# iterate samples
cut -f1 "$PROF" | LC_ALL=C sort -u | while read -r S; do
  tmp="$(mktemp)"
  awk -F'\t' -v S="$S" -v KEEP="$KEEP_UNMAPPED" '
    BEGIN{ OFS="\t" }
    # Load mapping: exact key is the species name (no internal id)
    NR==FNR {
      map[$1]=$2; next;
    }
    # Process rows for this sample
    $1==S {
      sp=$2; sub(/^[^_]+_/,"",sp);   # keep everything after first "_"
      taxid = (sp in map) ? map[sp] : ""
      abd = $3 + 0
      sum_all += abd
      if (taxid=="") {
        unm[sp] += abd
        if (KEEP) { taxid="0"; kept[taxid]+=abd; sum_kept+=abd }
      } else {
        agg[taxid] += abd
        sum_mapped += abd
      }
    }
    END{
      # nothing to output?
      norm = KEEP ? (sum_mapped+sum_kept) : sum_mapped
      if (norm <= 0) exit
      # mapped
      for (t in agg) printf "%s\t%.6f\n", t, 100.0*agg[t]/norm
      # optional kept (TaxId=0)
      if (KEEP) for (t in kept) printf "%s\t%.6f\n", t, 100.0*kept[t]/norm

      # brief unmapped log (top few)
      n=0; for (u in unm){ arr[++n]=u }
      if (n>0){
        lim = (n<15?n:15)
        print "[INFO] Sample " S ": unmapped (top " lim " by abundance):" > "/dev/stderr"
        # simple top-k: partial selection
        for(i=1;i<=lim;i++){
          imax=i; for(j=i+1;j<=n;j++) if (unm[arr[j]]>unm[arr[imax]]) imax=j
          tmp=arr[i]; arr[i]=arr[imax]; arr[imax]=tmp
          printf "  %s\t%.6f%%\n", arr[i], (unm[arr[i]]/(sum_all>0?sum_all:1))*100.0 > "/dev/stderr"
        }
      }
    }
  ' <($MAPREAD "$MAP") "$PROF" > "$tmp"

  if [[ -s "$tmp" ]]; then
    taxonkit profile2cami -p -s "$S" "$tmp" >> "$OUT"
  else
    echo "WARN: no mapped taxa for sample $S; skipping" >&2
  fi
  rm -f "$tmp"
done

echo "Done: wrote $(wc -l < "$OUT") lines to $OUT"

