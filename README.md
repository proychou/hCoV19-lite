# Pipeline for assembly of SARS-CoV-2 sequences
*Pavitra Roychoudhury*

This "lite" version is purely reference mapping-based. Compared to the [full version](https://github.com/proychou/hCoV19), no de novo assembly is performed, no annotation, no read filtering. SARS-CoV-2 sequences submitted to Genbank are automatically annotated with NCBI's pipeline.

Originally written to run on Fred Hutch HPC infrastructure with dependencies loaded via module load, but can be adapted to other platforms with the relevant dependencies installed (all open source). 

Dependencies: 
- Bowtie2/2.4.1
- FastQC/0.11.9
- R/4.0.0
- SAMtools/1.10
- BBMap/38.44

One-time steps: 
```bash
bowtie2-build './refs/NC_045512.2.fasta' ./refs/NC_045512.2
```

Example usage: 
For paired-end library
```bash 
covid_wgs_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz \
							-s samplename -aq -p ./refs/swift_primers.fasta
```
For single-end library
```bash	
covid_wgs_pipeline.sh -u yourreads.fastq.gz -s samplename -aq -p ./refs/swift_primers.fasta
```

Options:

* `-1` and `-2` : to specify R1 and R2 files for paired-end fastqs. Doesn't currently support interleaved fastqs. 
* `-u` : single-end fastq 
* `-s` : sample name
* `-a` : trim all Illumina adapters using list provided in bbduk
* `-q` : quality trimming 
* `-p primerlist.fasta` : remove primers from ends of reads with BBDuk using the list of primer sequences provided. Reads shorter than 75bp post-trim are discarded. Swift Biosciences primer list is provided in the references folder. 