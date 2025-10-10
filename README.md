# Micrite

> [!WARNING]  
> This repo is in early development and not yet ready for use

**micrite** is a collection of tools and workflows detecting and characterising microbes from cancer sequencing data.

micrite is a highly opinionated framework providing a stable API over several existing tools that do the real heavylifting (e.g. kraken2, blast, bwa, etc).

Configuration has been optimised for detecting microbial DNA/RNA from cancer biopsy samples characterised by NGS. In these samples, most will be human, with only small amounts of microbial DNA present. The goal is typically to detect strong signals of infective viruses/bacteria.

## Why another pipeline

mitrite provides a single interface for DNA and total-RNA based analyses, and provides a variety of approaches depending what strategy your upstream sequencing / alignment has employed and what your goals are (quick detection, slow detection).

## Installation

The easiest way to get started with micrite pipelines is just run the [nextflow pipeline](https://github.com/selkamand/nf-micrite).

To install locally (not using nextflow) just download the latest release for your OS and add to your path.

To install the dev version 

```
cargo install --git https://github.com/selkamand/micrite
```

## Quick Start

We break operations into two groups. **Screen** for microbial detection, **Sleuth** for validation and **Subtype** for subtyping.

```
# display help
micrite --help

# Screen a bam for microbes in the supplied database. If any hits are found, extract reads 
micrite screen --extract-reads 10376 --db <kraken_database> <bam>


# Search bam for microbes using kraken database
micrite screen full kraken --db "plusPF" --taxids 10376,10407 <bam>

# Search bam for microbes using metagenome realignment (where metagenome is created by downloading taxid genomes from NCBI)
micrite screen full t2tplusbugs --taxids 10376,10407 <bam>

# Search bam for microbes using metagenome realignment (where metagenome is supplied by user
micrite screen full t2tplusbugs --metagenome metagenome.fa <bam>

# Ensure microbe is really the genome we expect
micrite sleuth --taxid 10376 <putative_ebv_reads.fa>


# Using the results
micrite subtype --microbe EBV <putative_ebv_reads.fa>
```

## Detailed Description

### Screen

Rapidly scan an alignment file for microbes

**Types of Screens**

There are 3 different screens depending on your situation (and desired speed)

#### 1. Microbial reference genomes included during alignment

In this situation you can use the _superfast_ OR _quick_ subcommands

**superfast**: an idxstats based screen

Requirements:

- The reference genome includes your microbe of interest

**quick**: Based on filtering reads alinged to mappable parts of microbial reference genomes

Requirements:

- The reference genome includes your microbe of interest.
- Micrite includes a bedfile of hard-to-map regions to exclude

#### microbial reference genomes not included during alignment

**full**:

Collect reads from bam which are

1. unmapped (or partly unmapped?)
2. aligned to decoy contigs
3. mapped to regions of genome that show high homology to common viruses / bacteria (e.g. from common EBV decoy contigs)

Once we have these unmapped reads do 1 of two processes.

1. (kraken) Realign to chm13 (telomere to telomere reference genome) pull the unmapped reads again, and run through a kraken2.
2. (t2tplusbugs) Realign to a metagenome comprised of chm13 plus a small collection of microbes we're interested in charactising.

### Sleuth:

**Sleuth** carefully characterises suspected 'microbial' to confirm the microbe is really present. This is a slow process and should only be done once a **screen** has thrown up a putative hit in a microbe we care about.

### Subtype

Some microbes, e.g. EBV, can be subtyped if coverage is sufficient.

### Sift

Pull reads belonging to a specific microbe. 

### Testing Data

micrite comes packaged with a mini test kraken database ([krakendb](testfiles/krakendb)) comprised of human chromosome 2 from the T2T assembly + 3 viral genomes (Human gammaherpesvirus 8, Human gammaherpesvirus 4 (EBV), and Human papillomavirus).

We also include a mini test [deacon host-depletion database](testfiles/deacon/chm13v2.chr2.deacon.k57w131.idx) comprised of human chromosome 2 form the T2T assembly. 


## Other Tools

**[virus-interpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter)**

The best option for detecting viruses from tumour bams particularly if you're running the hmftools pipeline

**[viralrecon](https://nf-co.re/viralrecon/2.6.0)**

Great if you have just 1 virus of interest and want to know if its in your fastqs. Nextflow pipeline.

[**Hecatomb**](https://hecatomb.readthedocs.io/en/latest/)

Great for detecting many viruses. Very modular, and less opinionated than micrite. Snakemate pipeline.

### For Maintainers

A collection of quick commands to test micrite functions are working appropriately

**1. Run kraken screen**

```{r}
cargo run -- --outdir outdir screen --db-host testfiles/deacon/chm13v2.chr2.deac
on.k57w131.idx --db-kraken testfiles/database/krakendb/ --report-zero-counts --bam t
estfiles/humanGRCh38_9000_ebv_1000_hpv16_1000_hpylori_1000.grch38_noalt.bam
```

**1. Extract reads from a particular species**

Must have run kraken screen first

```{r}
cargo run -- sift \
--taxid 10376 \
--prefix ebv_reads \
--kout outdir/humanGRCh38_9000_ebv_1000_hpv16_1000_hpylori_1000.grch38_noalt.unmapped.hostdepleted.kout.tsv \
--fasta outdir/humanGRCh38_9000_ebv_1000_hpv16_1000_hpylori_1000.grch38_noalt.unmapped.hostdepleted.fasta \
--kreport outdir/humanGRCh38_9000_ebv_1000_hpv16_1000_hpylori_1000.grch38_noalt.unmapped.hostdepleted.kreport \
--outdir outdir/
```