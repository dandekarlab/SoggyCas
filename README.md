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


4) GFF File : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/411/555/GCF_001411555.1_wgs.5d/GCF_001411555.1_wgs.5d_genomic.gff.gz


5) Protein Table : https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=GetFeatures4Grid&amp;download=on&amp;type=Proteins&amp;genome_id=17683&genome_assembly_id=255039&mode=2&is_locus=1&is_locus_tag=0&optcols=1,0,,0,0

Extract these files in the working directory.


## Usage

### Step1 : Extract gene features

Parse GBFF and gene2refseq files to extract gene features

```bash
perl STEP1_GBFFParser.pl -f GBFFfile -r gene2refseqfile
```

The result file can also download from the link : https://drive.google.com/open?id=1YKD37g7Sv7JKYSeCVtwBg0wA2whIx-y1

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


