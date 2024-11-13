#!/bin/bash
[ -d "$1" ] && [ -f "$2" ] || { echo "USAGE: $0 <MarkerDB> <sample1.fq> ..." ; exit 1; }
PRO_D="."
KSSD="$PRO_D/bin/metakssd"
SCRIPT_D="$PRO_D/scripts"
SCRIPT_NAME="possion.kssd2out.pl"
SHUF_F="$PRO_D/shuf_files/L3K11.shuf"

MARKERDB=`dirname $1"/any"` 
shift

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
$KSSD dist -L $SHUF_F -A -o $OUT_SKTCH $@
[ -d $OUT_SKTCH ] || { echo "aborted!"; exit 1 ;}
end_time=$(date +%s)
echo "... Samples sketch created. Elapsed time: $((end_time - start_time)) seconds. <<"
echo ">> MetaKSSD profiling ..."
$KSSD composite -r $MARKERDB -q $OUT_SKTCH > $RAW_RSLT_F
perl $SCRIPT_D/$SCRIPT_NAME $RAW_RSLT_F 18  > $OUT_PROFILE_F
end_time=$(date +%s)
echo "... MetaKSSD profile created!: $OUT_PROFILE_F. Elapsed time: $((end_time - start_time)) seconds."
rm -r $TMP_D



