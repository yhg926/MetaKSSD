KSSD="./metakssd"
SCRIPT_D="./src"
SCRIPT_NAME="genome_species_labeling.pl"

SHUF_F="./shuf_files/L3K11.shuf"
OUT_D=$1"_L3K11_sketch"
PAN_D=$OUT_D"_pan"
SPSP_D=$PAN_D"_union_sp"
MARKERDB=$OUT_D"_markerdb"

GLIST_F="./g.list"
GENOMES_D=$1
G2TAXONOMY_F=$2
GROUPING_F="./tax_group.tsv"

$KSSD dist -L $SHUF_F -o $OUT_D $GENOMES_D
$KSSD set -P $OUT_D > $GLIST_F
perl $SCRIPT_D/$SCRIPT_NAME $GLIST_F $G2TAXONOMY_F >	$GROUPING_F

$KSSD set -g $GROUPING_F -o $PAN_D $OUT_D
$KSSD set -q -o $SPSP_D $PAN_D 
$KSSD set -i $SPSP_D -o $MARKERDB $PAN_D




