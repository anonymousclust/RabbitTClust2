# `RabbitTClust v.2.3.0`

## Installation
`RabbitTClust v.2.3.0` can only support 64-bit Linux Systems.




### Install from source code
#### Dependancy
* cmake v.3.0 or later
* c++14
* [zlib](https://zlib.net/)
* icc

#### Compile and install
```bash
git clone https://github.com/anonymousclust/RabbitTClust2.git -b fastmpi
cd RabbitTClust2
./install.sh
```
## Usage
```bash
The MPI version only support clust-mst --fast -l -i  and clust-mst --presketched for large-scale genome datasets clustering
# clust-mst, minimum-spanning-tree-based module for RabbitTClust
Usage: ./clust-mst [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -c,--containment INT        use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress  -k,--kmer-size INT          set the kmer size
  -s,--sketch-size INT        set the sketch size for Jaccard Index and Mash distance, default 1000
  -l,--list                   input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches or MST
  -d,--threshold FLOAT        set the distance threshold for clustering
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT Excludes: --append
                              set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
  --fast                      use the kssd algorithm for sketching and distance computing for clust-mst

```

## Example:
```bash

mpirun -np 2 ./clust-mst --fast -l -i 10000.list -o 32.clust -t 32


mpirun -np 4 ./clust-mst --fast --presketched 2025_08_01_05-00-04/ -o 32.clust -t 8


```
## Output
The output file is in a CD-HIT output format and is slightly different when running with or without `-l` input option.  
When using the `-l` option, the input is expected to be a FASTA file list, with each file representing a genome. Without the `-l` option, the input should be a single FASTA file, with each sequence representing a genome.

#### Output format for a FASTA file list input
With `-l*` option, the tab-delimited values in the lines beginning with tab delimiters are:
* local index in a cluster
* global index of the genome
* genome length
* genome file name (including genome assembly accession number)
* sequence name (first sequence in the genome file)
* sequence comment (remaining part of the line)

**Example:**
```txt
the cluster 0 is:
    0   0   14782125nt  bacteria/GCF_000418325.1_ASM41832v1_genomic.fna     NC_021658.1     Sorangium cellulosum So0157-2, complete sequence
    1   1   14598830nt  bacteria/GCF_004135755.1_ASM413575v1_genomic.fna    NZ_CP012672.1   Sorangium cellulosum strain So ce836 chromosome, complete genome

the cluster 1 is:
    0   2   14557589nt  bacteria/GCF_002950945.1_ASM295094v1_genomic.fna    NZ_CP012673.1   Sorangium cellulosum strain So ce26 chromosome, complete genome

the cluster 2 is:
    0   3   13673866nt  bacteria/GCF_019396345.1_ASM1939634v1_genomic.fna   NZ_JAHKRM010000001.1    Nonomuraea guangzhouensis strain CGMCC 4.7101 NODE_1, whole genome shotgun sequence

......
```

#### Output format for a single FASTA file input
Without `-l` option, the tab-delimited values in the lines beginning with tab delimiters are:
* local index in a cluster
* global index of the genome
* genome length
* sequence name 
* sequence comment (remaining part of this line)

**Example:**
```txt
the cluster 0 is:
    0   0   11030030nt  NZ_GG657755.1   Streptomyces  himastatinicus ATCC 53653 supercont1.2, whole genome shotgun sequence
    1   1   11008137nt  NZ_RIBZ01000339.1   Streptomyces  sp. NEAU-LD23 C2041, whole genome shotgun sequence

the cluster 1 is:
    0   2   11006208nt  NZ_KL647031.1   Nonomuraea  candida strain NRRL B-24552 Doro1_scaffold1, whole genome shotgun sequence
    
the cluster 2 is:
    0   3   10940472nt  NZ_VTHK01000001.1   Amycolatopsis anabasis strain EGI 650086 RDPYD18112716_A.Scaf1, whole genome shotgun sequence

......
```

#### Output the newick tree format (v.2.2.1 or latter)
When the `--newick-tree` option is used, an additional output file will be generated in the Newick tree format with a suffix name of ".newick.tree".


# Bug Report
We highly appreciate all bug reports, comments, and suggestions from our users.  
Please feel free to raise any concerns or feedback with us without hesitation by `issue`. 

