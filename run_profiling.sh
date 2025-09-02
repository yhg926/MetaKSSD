#!/bin/bash

citations=$(cat <<'EOF'
In keeping with academic convention, please cite the works on which your article is based.

If you use MetaKSSD, please cite:

Yi, H., Lu, X., & Chang, Q. (2025). MetaKSSD: boosting the scalability of the reference
taxonomic marker database and the performance of metagenomic profiling using sketch operations.
Nature Computational Science. https://doi.org/10.1038/s43588-025-00855-0

If you also refer to the K-mer Substring Space Sampling and Decomposition (KSSD) method, please
additionally cite:

Yi, H., Lin, Y., Lin, C., et al. (2021). Kssd: sequence dimensionality reduction by k-mer substring
space sampling enables real-time large-scale datasets analysis. Genome Biology, 22, 84.
https://doi.org/10.1186/s13059-021-02303-4
EOF
)

echo "$citations"

[ -d "$1" ] && [ -f "$2" ] || { echo "USAGE: $0 <MarkerDB> <sample1.fq> ..." ; exit 1; }
PRO_D=$METAKSSD_PATH
KSSD="$PRO_D/bin/metakssd"
SCRIPT_D="$PRO_D/scripts"
SCRIPT_NAME="possion.kssd2out.pl"
SHUF_F="$PRO_D/shuf_files/L3K11.shuf"

MARKERDB=`dirname $1"/any"` 

[ -f $KSSD ] || { echo "$KSSD does not exists" ; exit 1 ;}
[ -f $SCRIPT_D/$SCRIPT_NAME ] || { echo "$SCRIPT_D/$SCRIPT_NAME does not exists" ; exit 1; }
[ -f $SHUF_F ] || { echo "$SHUF_F does not exists" ; exit 1 ;}
[ -d $MARKERDB ] || { echo "$MARKERDB does not exists" ; exit 1; }

TMP_D="./tmp_dir"
mkdir -p $TMP_D
OUT_SKTCH="$TMP_D/Samples_L3K11_sketch"
RAW_RSLT_F="$TMP_D/species_coverage.tsv"
OUT_PROFILE_F=$2"_metakssd.profile.tsv"
echo ">> MetaKSSD sektching  ..."
start_time=$(date +%s)
shift 
$KSSD dist -L $SHUF_F -A -o $OUT_SKTCH $@
[ -d $OUT_SKTCH ] || { echo "aborted!"; exit 1 ;}
end_time=$(date +%s)
echo "... Samples sketch created. Elapsed time: $((end_time - start_time)) seconds. <<"
echo ">> MetaKSSD profiling ..."
$KSSD composite -r $MARKERDB -q $OUT_SKTCH > $RAW_RSLT_F
perl $SCRIPT_D/$SCRIPT_NAME $RAW_RSLT_F 18  > $OUT_PROFILE_F
end_time=$(date +%s)
echo "... MetaKSSD profile created!: $OUT_PROFILE_F. Elapsed time: $((end_time - start_time)) seconds."
#rm -r $TMP_D



