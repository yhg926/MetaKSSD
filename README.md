#  Instantaneous Metagenome Taxonomic Profiling with MetaKSSD

MetaKSSD is the second version of KSSD (K-mer Substring Space Sampling/Shuffling Decomposition), designed for instantaneous metagenome taxonomic profiling using WGS fastq data.

K-mer Substring Space Decomposition (KSSD) facilitates highly efficient genome sketching and enables lossless sketch operations, including union, intersection, and subtraction [doi.org/10.1186/s13059-021-02303-4]. Building upon the KSSD framework, MetaKSSD introduces a novel feature that tracks k-mer counts within the sketch. Leveraging these foundational functionalities, MetaKSSD further innovates methods for constructing a taxonomic marker database (MarkerDB), metagenome taxonomic profiling, and profile searching.

# 1. Installation 
```
git clone https://github.com/yhg926/MetaKSSD.git &&
cd MetaKSSD &&
make
```
## 1.1 (Optional) Get pre-built MarkerDB (L3K11)

```
wget http://www.genomesketchub.com/download/markerdb.L3K11_gtdb_r214.tar.gz
tar xf markerdb.L3K11_gtdb_r214.tar.gz
```
This step is only needed when you do not have a MarkerDB. 
You can also prepare your own MarkerDB, 
see [build custom MarkerDB](#5-build-custom-MarkerDB).

## 1.2 (Optional) Prepare gtdbr214 to ncbi taxonomy convertion tables 
```
gunzip -d data/*.gz;
```
These files are only needed when you have to convert gtdb to ncbi taxonomy.

## 1.3 (Optional) Get pre-build Abundance Vector Database (L3K11)
```
wget http://www.genomesketchub.com/download/markerdb.abvdb231227.L3K11_gtdb_r214.tar.gz
tar xf markerdb.abvdb231227.L3K11_gtdb_r214.tar.gz
```
The Abundance Vector Database is only needed when you perform abundance vector searching.
You can also prepare your own Abundance Vector Database, 
see [Index abundance vector database](#4-Index-abundance-vector-database).

# 2. Metagenome profiling

```
#sketching with k-mer counts tracking
metakssd dist -L shuf_files/L3K11.shuf -A -o <sample1_sketch> <sample1.fastq>
#generate raw profile
metakssd composite -r <markerdb> -q <sample1_sketch> > <species_coverage.tsv>
#abundance normalization
perl src/possion.kssd2out.pl <species_coverage.tsv> <minimum overlapped k-mer S (default:18)> > <species relative abundance profile>
```
If need to covert species abundaces to CAMI format profile, using 
```
#format convertion 
perl src/possion.kssdcomposite2taxonomy_profilefmt.pl <species_coverage.tsv> data/best.gtdbr214_psid2ncbi_specid.tsv data/scienficaname.ncbitaxid_rank_parentnode_name.gtdbr214_pseudoidrelated.tsv 18 > sample1.profile
```

# 3. Abundance Vector Searching 
## 3.1 Generate your abundance vector in a given path

```
metakssd composite -r <markerdb> -q <metagenome sketch> -b -o <path>
```
## 3.2 Abundance Vector Searching 
To retrieve abundance vectors similar to an abundance vector "input.abv" from the database "markerdb":

```
metakssd composite -r <markerdb> -s<0 or 1> <path/input.abv>
```
Here, the options -s0 and -s1 enable searching based on L1 norm and cosine similarity, respectively.

# 4. Index abundance vector database 
Suppose you have many  *.abv generated following [section 3.1](#31-Generate-your-abundance-vector-in-a-given-path).
Then the abundance vector database could be indexed as follow:
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

# 6. MetaKSSD benchmarking results

OPAL benchmarking results on five datasets are available:
1.	https://yhg926.github.io/KSSD2/OPAL/mouse_gut/
2.	https://yhg926.github.io/KSSD2/OPAL/marine/
3.	https://yhg926.github.io/KSSD2/OPAL/strain_madness/
4.	https://yhg926.github.io/KSSD2/OPAL/rhizosphere/
5.	https://yhg926.github.io/KSSD2/OPAL/new_released/







