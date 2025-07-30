#  Instantaneous Metagenomic Taxonomic Profiling with MetaKSSD

MetaKSSD is the second version of KSSD (K-mer Substring Space Sampling/Shuffling Decomposition), designed for instantaneous metagenome taxonomic profiling using WGS fastq data.

K-mer Substring Space Decomposition (KSSD) facilitates highly efficient genome sketching and enables lossless sketch operations, including union, intersection, and subtraction [doi.org/10.1186/s13059-021-02303-4]. Building upon the KSSD framework, MetaKSSD introduces a novel feature that tracks k-mer counts within the sketch. Leveraging these foundational functionalities, MetaKSSD further innovates methods for constructing a taxonomic marker database (MarkerDB), metagenome taxonomic profiling, and profile searching.

For users not familiar with linux command-line, we provide one-click metagenomic analysis app:
[MacOS MetaKSSD Clients](http://www.genomesketchub.com/download/MetaKSSD_Mac.zip ), see [tutorial video](https://youtu.be/-JyctdOiMO4)
[Windows OS MetaKSSD Clients] (https://youtu.be/ck5af1ewX4w), see [tutorial video](http://www.genomesketchub.com/download/MetaKSSD_Windows.exe)



The philosophy behind MetaKSSD: https://github.com/yhg926/MetaKSSD/wiki/Philosophy-behind-MetaKSSD 

# 1. Installation 
```
git clone https://github.com/yhg926/MetaKSSD.git &&
cd MetaKSSD &&
make
export PATH=$PATH:$(pwd)/bin
```
## 1.1 (Optional) Get pre-built MarkerDB (L3K11)

```
#latest version
wget https://zenodo.org/records/16317275/files/GTDBr226_genomes_L3K11_sketch_markerdb.tar.gz
# older version
wget https://zenodo.org/records/14353354/files/markerdb.gtdb_r214v2_L3K11_sketch.tar.gz
tar xf *markerdb*.tar.gz
```
This step is only needed when you do not have a MarkerDB. 
You can also prepare your own MarkerDB, 
see [build custom MarkerDB](#5-build-custom-MarkerDB).

## 1.2 (Optional) Prepare gtdb to ncbi taxonomy convertion tables 
```
gunzip -d data/best.gtdbr[214|226]_psid2ncbi_specid.tsv.gz;
gunzip -d data/scienficaname.ncbitaxid_rank_parentnode_name.gtdbr[214|226]_pseudoidrelated.tsv.gz
```
These files are only needed when you have to convert gtdb to ncbi taxonomy.

## 1.3 (Optional) Get pre-build Abundance Vector Database (L3K11)
```
#only r214 are avialble currently 
wget https://zenodo.org/records/11437234/files/markerdb.abvdb231227.L3K11_gtdb_r214.tar.gz
tar xf markerdb.abvdb231227.L3K11_gtdb_r214.tar.gz
```
The Abundance Vector Database is only needed when you perform abundance vector searching.
You can also prepare your own Abundance Vector Database, 
see [Index abundance vector database](#4-Index-abundance-vector-database).

# 2. Metagenome profiling

One stand profiling (GTDB taxonomy only)
```
./run_profiling.sh <MarkerDB> <sample1.fq> ...
```
Profiling breakdown (for user customization)
```
#sketching with k-mer counts tracking
metakssd dist -L shuf_files/L3K11.shuf -A -o <sample1_sketch> <sample1.fastq>
#generate raw profile
metakssd composite -r <markerdb> -q <sample1_sketch> > <species_coverage.tsv>
#abundance normalization
perl src/possion.kssd2out.pl <species_coverage.tsv> <minimum overlapped k-mer S (default:18)> > <species_relative_abundance_profile>
```
If need to covert species abundaces to full gtdb taxonomy profile, using 
```
perl scripts/kssd2out2gtdb_taxonomy_profile.pl <species_relative_abundance_profile> data/gtdbr[214|226]_psid2krona_taxonomy.tsv
```

If need to covert species_coverage.tsv to CAMI format profile with NCBI taxonomy, using 
```
#format convertion 
perl scripts/possion.kssdcomposite2taxonomy_profilefmt.pl <species_coverage.tsv> data/best.gtdbr[214|226]_psid2ncbi_specid.tsv data/scienficaname.ncbitaxid_rank_parentnode_name.gtdbr[214|226]_pseudoidrelated.tsv 18 > sample1.profile
```
If need to covert species_coverage.tsv to Krona format profile, using
```
perl scripts/kssdcomposite2gtdb_tax_kronafmt.pl <species_coverage.tsv> data/gtdbr[214|]_psid2krona_taxonomy.tsv <outdir>
```

# 3. Abundance Vector Searching 
## 3.1 Generate your abundance vector in a given path

```
metakssd composite -r <markerdb> -q <metagenome sketch> -b -o <path>
```
## 3.2 Abundance Vector Searching 
To retrieve abundance vectors similar to an abundance vector "input.abv" from a markerdb with indexed abv:

```
metakssd composite -r <markerdb.abvdb> -s<0 or 1> <path/input.abv>
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
Alternatively, we have provided an all-in-one workflow for MarKerDB construction based on GTDB r226/r214 genomes.
```
# build gtdbr226
sh build_MarkerDB_gtdbr226.sh <all_gtdbr226_genomes_dir>
# gtdbr214
sh build_MarkerDB.sh <all_gtdbr214_genomes_dir>
```

# 6. MetaKSSD benchmarking results

OPAL benchmarking results on five datasets are available (based on r214 markerdb):
1.	https://yhg926.github.io/KSSD2/OPAL/mouse_gut/
2.	https://yhg926.github.io/KSSD2/OPAL/marine/
3.	https://yhg926.github.io/KSSD2/OPAL/strain_madness/
4.	https://yhg926.github.io/KSSD2/OPAL/rhizosphere/
5.	https://yhg926.github.io/KSSD2/OPAL/new_released/


# 7. Related papers
Yi, H. MetaKSSD: Boosting the Scalability of Reference Taxonomic Marker Database and the Performance of Metagenomic Profiling Using Sketch Operations. bioRxiv 2024.06.21.600011 (2024) doi:10.1101/2024.06.21.600011.








