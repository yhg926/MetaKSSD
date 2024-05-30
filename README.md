#  Instantaneous Metagenome Taxonomic Profiling with MetaKSSD

MetaKSSD is version 2 of KSSD (K-mer substring space sampling/shuffling Decomposition).

K-mer substring space decomposition (KSSD) facilitates highly efficient genome sketching and enables lossless sketch operations including union, intersection and subtraction [doi.org/10.1186/s13059-021-02303-4]. Building upon the KSSD framework, MetaKSSD introduce a novel feature that tracks k-mer counts within the sketch. Leveraging these foundational functionalities, MetaKSSD further innovates in the methods of taxonomic marker database (MarkerDB) construction, metagenome taxonomic profiling and profile searching. 

# 1. Installation 
```
git clone https://github.com/yhg926/MetaKSSD.git &&
cd MetaKSSD &&
make
```
# 2. Metagenome profiling
make sure you already have an <markerdb>. If not, skip to [build custom MarkerDB](#5-build-custom-MarkerDB) to create one.
```
#sketching with k-mer counts tracking
metakssd dist -L <*.shuf> -A -o <sample1_sketch> <sample1.fastq>
#profiling
kssd composite -r <markerdb> -q <sample1_sketch> > <species_coverage.tsv>
# abundance normalization
perl possion.kssd2out.pl < species_coverage.tsv > <minimum overlapped k-mer > > <species relative abundance profile>
```
# 3. Abundance Vector Searching 
Make sure you already have an indexed abundance vector database. If not, skip to [Index abundance vector database](#4-Index-abundance-vector-database) to create one
```
#generate abundance vector in a given path 
metakssd composite -r <markerdb> -q <metagenome sketch> -b -o <path>
```
To retrieve abundance vectors similar to an abundance vector "input.abv" from the database "markerdb", we utilized the following command:

```
metakssd composite -r <markerdb> -s<0 or 1> <path/input.abv>
```
Here, the options -s0 and -s1 enable searching based on L1 norm and cosine similarity, respectively.

# 4. Index abundance vector database 
```
#make folder named abundance_Vec under your markerdb path
mkdir -p <markerdb path>/abundance_Vec
#collect all *.abv to the folder
cp *.abv <markerdb path>/abundance_Vec
#index 
metakssd composite -r <markerdb path> -i
```

# 5. build custom MarkerDB
```
# sketching reference genomes
metakssd dist -L <*.shuf> -o <L3K11_sketch> <all genomes Dir>
# print genome name
metakssd set -P <L3K11_sketch> > <genome_name.txt>

```
The genome names within the sketch were print line by line and direct the output to a file named "genome_name.txt" .
Then a grouping file named "group_name.txt" were prepared as follow: The "group_name.txt" file should contain the same lines as "genome_name.txt", where each line specifies the species name of the genome corresponding to the line in "genome_name.txt". Each line follows the format: "ID<TAB>species_name\n", where the ID is a non-negative integer uniquely labeling the species_name, and the species name of the genome can be found in the metadata files ".*_metadata_r214.tar.gz". For example, "1 Escherichia_coli\n" represents the species name "Escherichia coli" labeled with ID 1. The ID 0 is reserved for excluding genomes from the resulting sketch. Once "group_name.txt" is prepared, genomes can be grouped by their originating species using the following command:
```
metakssd set -g <group_name.txt> -o <L3K11_pan-sketch> <L3K11_sketch>
```
Here, "L3K11_pan-sketch" represents the consolidated ‘pangenome’ sketches for all species.

Subsequently, the union sketch of all species-specific marker, named "L3K11_union_sp-sketch", was obtained using this command:

```
metakssd set -q -o <L3K11_union_sp-sketch> <L3K11_pan-sketch>
```
Finally, the MarkerDB named “markerdb_L3K11” was generated by overlapping "L3K11_union_sp-sketch" with "L3K11_pan-sketch" using this command:

```
metakssd set -i <L3K11_union_sp-sketch> -o <markerdb_L3K11> <L3K11_pan-sketch>
```









