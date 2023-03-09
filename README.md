# saturation.R

Here is an R script `saturation.R` for estimating sequencing saturation from a
GEX, VDJ, or ADT dataset from the 10xgenomics platform.

The script uses the binomial distribution to downsample the reads and estimate
a saturation curve. This can be helpful to determine if a sequencing experiment
has enough reads.

[10xgenomics.com][1] gives us the formula for sequencing saturation:

```
Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
```

[1]: https://kb.10xgenomics.com/hc/en-us/articles/115003646912

Learn more about sequencing saturation from the 10xgenomics.com documentation:

- [What is sequencing saturation?](https://kb.10xgenomics.com/hc/en-us/articles/115005062366)
- [How is sequencing saturation calculated?][1]
- [How much sequencing saturation should I aim for?](https://kb.10xgenomics.com/hc/en-us/articles/115002474263)

File formats:

- [What is the molecule_info.h5 file?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) (GEX)
- [What is the all_contig_annotations.csv file?](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/overview) (VDJ)
- [What is the sample.stat.csv.gz file?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/antibody) (ADT)

## Getting started

See [output](output) for an example of the output files.

Install the dependencies:

```r
install.packages(
  c("data.table", "ggplot2", "glue", "optparse", "pbapply", "scales", "stringr", "BiocManager")
)
BiocManager::install("rhdf5")
```

## Usage examples

### GEX

```bash
Rscript saturation.R --out output --file molecule_info.h5
```
```
Reading molecule_info.h5
Estimating GEX saturation for 519882 barcodes
Writing output/saturation-gex.tsv
Writing output/total_reads-vs-saturation-gex.pdf
```

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224153589-ef7b1580-c29e-43e6-938a-d4d9e66af541.png">
</p>


### VDJ

```bash
Rscript saturation.R --out output/tcr --file all_contig_annotations.csv
```
```
Reading all_contig_annotations.csv
Estimating VDJ saturation for 20233 barcodes
Writing output/tcr/saturation-vdj.tsv
Writing output/tcr/total_reads-vs-saturation-vdj.pdf
```

<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224153652-933585f5-3e48-4a35-8514-cc2f3b7d339e.png">
</p>

### ADT

```bash
Rscript saturation.R --out output --file sample.stat.csv.gz
```
```
Reading Batch_1A_ADT.stat.csv.gz
Estimating ADT saturation
Writing output/saturation-adt.tsv
Writing output/total_reads-vs-saturation-adt.pdf
Writing output/saturation-adt-feature.tsv
Writing output/histogram-saturation-adt-feature.pdf
```
<p align="center">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224153734-6fe76c17-0aef-487c-8e3f-de70a1540c86.png">
<img width="50%" src="https://user-images.githubusercontent.com/209714/224153888-63613f19-0840-4d79-8a2f-5ebca0ffafb3.png">
</p>
