#!/bin/bash

 [ -d "$1" ] || { echo "USAGE: $0 <all_gtdb_genomes_dir>" ; exit 1; }

PRO_D="."
KSSD="$PRO_D/bin/metakssd"
SCRIPT_D="$PRO_D/scripts"
SCRIPT_NAME="genome_species_labeling.pl"
SHUF_F="$PRO_D/shuf_files/L3K11.shuf"
GENOMES_D=`dirname $1"/any"`

[ -f $KSSD ] || { echo "$KSSD does not exists" ; exit 1 ;}
[ -f $SCRIPT_D/$SCRIPT_NAME ] || { echo "$SCRIPT_D/$SCRIPT_NAME does not exists" ; exit 1; }
[ -f $SHUF_F ] || { echo "$SHUF_F does not exists" ; exit 1 ;}
[ -d $GENOMES_D ] || { echo "$GENOMES_D does not exists" ; exit 1; }

OUT_D=$GENOMES_D"_L3K11_sketch"
PAN_D=$OUT_D"_pan"
SPSP_D=$PAN_D"_union_sp"
MARKERDB=$OUT_D"_markerdb"

TMP_D="./tmp_dir"
mkdir -p $TMP_D

GLIST_F="./$TMP_D/g.list"

G2TAXONOMY_GZ="$PRO_D/data/*_taxonomy_r226.tsv.gz"
G2TAXONOMY_F="./$TMP_D/genome2taxonomy.tsv"
GROUPING_F="./$TMP_D/tax_group.tsv"

echo ">> Running Step 1. KSSD sektching $GENOMES_D ..."
start_time=$(date +%s)
$KSSD dist -L $SHUF_F -o $OUT_D $GENOMES_D
[ -d $OUT_D ] || { echo "Step 1 aborted"; exit 1 ;}

end_time=$(date +%s)
echo "... All genomes sketch created. Elapsed time: $((end_time - start_time)) seconds. <<"

echo ">> Running Step 2. Labeling GTDB taxonomy ..."
$KSSD set -P $OUT_D > $GLIST_F
gunzip -c $G2TAXONOMY_GZ > $G2TAXONOMY_F || { echo "$G2TAXONOMY_GZ does not exists!" ; exit 1; }
perl $SCRIPT_D/$SCRIPT_NAME $GLIST_F $G2TAXONOMY_F >$GROUPING_F

[ -f $GROUPING_F ] || { echo "Step 2 aborted!" ; exit 1; }

end_time=$(date +%s)
echo ">> ... GTDB taxonomy labeled. Elapsed time: $((end_time - start_time)) seconds. << "

echo ">> Running Step 3. Generating pangenome sketch ..."
$KSSD set -g $GROUPING_F -o $PAN_D $OUT_D

[ -d $PAN_D ] || { echo "Step 3 aborted" ;exit 1; }

end_time=$(date +%s) 
echo ">> ... All pangenome sketches obtained. Elapsed time: $((end_time - start_time)) seconds. << "

echo ">> Running Step 4. Generating species-specific sketch ..."
end_time=$(date +%s)
$KSSD set -q -o $SPSP_D $PAN_D 
$KSSD set -i $SPSP_D -o $MARKERDB $PAN_D

[ -d $MARKERDB ] || { echo "Step 4 aborted"; exit 1; }


echo ">> ... MarkerDB $MARKERDB constructed. Elapsed time: $((end_time - start_time)) seconds. << "

rm -r $TMP_D
