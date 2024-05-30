#  Instantaneous Metagenome Taxonomic Profiling with MetaKSSD

MetaKSSD is version 2 of KSSD (K-mer substring space sampling/shuffling Decomposition).
K-mer substring space decomposition (KSSD) facilitates highly efficient genome sketching and enables lossless sketch operations including union, intersection and subtraction [doi.org/10.1186/s13059-021-02303-4]. Building upon the KSSD framework, MetaKSSD introduce a novel feature that tracks k-mer counts within the sketch. Leveraging these foundational functionalities, MetaKSSD further innovates in the methods of taxonomic marker database (MarkerDB) construction, metagenome taxonomic profiling and profile searching. 

# 1 Installation 
```
git clone https://github.com/yhg926/MetaKSSD.git &&
cd MetaKSSD &&
make
```
# 2 Metagenome profiling
```
#sketching 
metakssd dist -L <*.shuf> -A -o <sample1_sketch> <sample1.fastq>
#profiling
kssd composite -r <markerdb_L3K11> -q <sample1_sketch> > <species_coverage.tsv>
# abundance normalization
perl possion.kssd2out.pl < species_coverage.tsv > <minimum overlapped k-mer > > <species relative abundance profile>
```

