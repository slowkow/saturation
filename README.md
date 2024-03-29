# saturation.R

<a href="https://zenodo.org/badge/latestdoi/611932244"><img src="https://zenodo.org/badge/611932244.svg" alt="DOI"></a>

Table of contents:
* [Introduction](#introduction)
* [Getting started](#getting-started)
* [Usage examples](#usage-examples)
  * [GEX](#gex)
  * [VDJ](#vdj)
  * [ADT](#adt)
* [Running parallel jobs with rush](#running-parallel-jobs-with-rush)

## Introduction

Here is an R script [`saturation.R`](saturation.R) for estimating sequencing saturation from a
GEX, VDJ, or ADT dataset from the 10x Genomics platform.

The script uses the binomial distribution to downsample the reads and estimate
a saturation curve. This can be helpful to determine if a sequencing experiment
has enough reads.

[10xgenomics.com][1] gives us this formula for sequencing saturation:

```
Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
```

[1]: https://kb.10xgenomics.com/hc/en-us/articles/115003646912

Here is my illustration of the relationship between sequencing saturation and reads per [unique molecular identifier](https://dnatech.genomecenter.ucdavis.edu/faqs/what-are-umis-and-why-are-they-used-in-high-throughput-sequencing/) (UMI):

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/231836618-a6cb746b-be7f-4605-8c84-c419cf48aab6.png">
</p>

We can compute reads per UMI from the saturation, and vice versa:

```r
d <- data.frame(sat = seq(0, 1, length.out = 1001))
d$rpu <- 1 / (1 - d$sat)
```

Learn more about sequencing saturation from the 10xgenomics.com documentation:

- [What is sequencing saturation?](https://kb.10xgenomics.com/hc/en-us/articles/115005062366)
- [How is sequencing saturation calculated?][1]
- [How much sequencing saturation should I aim for?](https://kb.10xgenomics.com/hc/en-us/articles/115002474263)

File formats:

- [What is the molecule_info.h5 file?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) (GEX)
- [What is the all_contig_annotations.csv file?](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/overview) (VDJ)
- [What is the sample.stat.csv.gz file?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/antibody) (ADT)

## Getting started

Install the dependencies:

```r
install.packages(
  c("data.table", "ggplot2", "ggtext", "glue", "optparse", "pbapply", "scales", "stringr", "BiocManager")
)
BiocManager::install("rhdf5")
```

See [output](output) for an example of the output files.


## Usage examples

### GEX

```bash
Rscript saturation.R --out output --file molecule_info.h5
```
```
Reading molecule_info.h5
Estimating GEX saturation for 519882 barcodes
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
Writing output/saturation-gex.tsv
Writing output/total_reads-vs-saturation-gex.pdf
```

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224188198-19f808f7-cbe9-4c21-a88a-8bbf9eb4cec7.png">
</p>


### VDJ

```bash
Rscript saturation.R --out output/tcr --file all_contig_annotations.csv
```
```
Reading all_contig_annotations.csv
INFO: Removing '-1' from the end of each barcode
Estimating VDJ saturation for 20233 barcodes
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
Writing output/tcr/saturation-vdj.tsv
Writing output/tcr/total_reads-vs-saturation-vdj.pdf
```

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/225680974-4607e8df-9925-486d-ba00-4aebe51da22f.png">
</p>

### ADT

```bash
Rscript saturation.R --out output --file Batch_1A_ADT.stat.csv.gz
```
```
Reading Batch_1A_ADT.stat.csv.gz
Estimating ADT saturation curve for 729372 barcodes
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
Writing output/saturation-adt.tsv
Writing output/total_reads-vs-saturation-adt.pdf
Writing output/saturation-adt-feature.tsv
Writing output/histogram-saturation-adt-feature.pdf
```

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224188298-b136b303-fcf1-4c73-a764-6e64ee178d0c.png">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224188339-428709d8-c885-4d28-85f1-16159b58821f.png">
</p>

## Running parallel jobs with rush

Install [rush](https://github.com/shenwei356/rush) by Wei Shen:

```bash
go install github.com/shenwei356/rush
```

Then make a list of input files and pass it to rush:

```bash
ls /project/cellranger_output/*/{molecule_info.h5,all_contig_annotations.csv,*.stat.csv.gz} > files.txt

# Run 16 jobs in parallel, capture outputs from each job in one file
rush -i files.txt -o rush-saturation.txt -j16 'Rscript saturation.R --out out/{/%} --file {}'
```

