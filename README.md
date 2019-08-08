# SoggyCas
CRISPR/CAS gRNA design tool for Walnut

SoggyCas is a bioinformatics tool to aid selection of guide RNAs targeting specific genomic sites and to identify potential off-target sites in Walnut. SoggyCas is developed using [Perl](https://www.perl.org/).

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

## Prerequisites

Perl Modules

[BioPerl](https://metacpan.org/release/BioPerl)

Bio::SeqFeature::Tools::Unflattener

Bio::SeqIO

Bio::Restriction::EnzymeCollection


## Download gene model & fasta file

Organism : [Juglans regia](https://www.ncbi.nlm.nih.gov/genome/?term=juglans+regia)

1) GBFF : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/411/555/GCF_001411555.1_wgs.5d/GCF_001411555.1_wgs.5d_genomic.gbff.gz


2) Fasta : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/411/555/GCF_001411555.1_wgs.5d/GCF_001411555.1_wgs.5d_genomic.fna.gz


3) gene2refseq : ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz


## Usage

### Step1 : Extract gene features

Parse GBFF and gene2refseq files to extract gene features

```bash
perl STEP1_GBFFParser.pl -f GBFFfile -r gene2refseqfile
```

### Step2 : Generate tables for annotation

Extract exon,intron,UTR information

```bash
perl SETP2_GenerateTables.pl
```

### Step3 : Extract gRNAs for a given gene

Parse GBFF and gene2refseq files to extract gene features

```bash
perl STEP3_FindgRNA.pl -i geneid
```



## SoggyCas Workflow

![Bioinformatics Workflow](https://github.com/dandekarlab/SoggyCas/blob/master/Workflow.png)


