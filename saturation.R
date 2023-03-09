#!/usr/bin/env Rscript
#
# saturation.R
#
# MIT License
# 
# Copyright (c) 2023 Kamil Slowikowski
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Usage:
#   Rscript saturation.R --out outdir --file molecule_info.h5
#   Rscript saturation.R --out outdir --file all_contig_annotations.csv
#   Rscript saturation.R --out outdir --file all_contig_annotations.csv --barcodes barcodes.txt 
#   Rscript saturation.R --out outdir --file sample.stat.csv.gz
#
# The input file must be one of:
#
#   molecule_info.h5            for GEX data
#   all_contig_annotations.csv  for VDJ data
#   *.stat.csv.gz               for ADT data
#
# The --barcodes file is optional. It does not make a big difference, so don't
# worry if you don't have one.
#
# The following output files will be written:
#
#   GEX
#     saturation-gex.tsv
#     total_reads-vs-saturation-gex.pdf
#   
#   VDJ
#     saturation-vdj.tsv
#     total_reads-vs-saturation-vdj.pdf
#
#   VDJ (with --barcodes)
#     saturation-vdj.tsv
#     total_reads-vs-saturation-vdj.pdf
#     saturation-pct_cdr3.tsv
#     total_reads-vs-pct_cdr3.pdf
#   
#   ADT
#     saturation-adt.tsv
#     saturation-adt-feature.tsv
#     total_reads-vs-saturation-adt.pdf
#     histogram-saturation-adt-feature.pdf
#
# Install the dependencies:
#
#   install.packages(c("data.table", "ggplot2", "glue", "optparse", "pbapply", "scales", "stringr", "BiocManager"))
#   BiocManager::install("rhdf5")

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
  library(optparse)
  library(pbapply)
  library(rhdf5)
  library(scales)
  library(stringr)
})
status <- function(x) {
  message(Sys.time(), '\t', x)
}
theme_kamil <- theme_classic(
  base_size = 16,
  base_line_size = 0.3,
  base_rect_size = 0.3
) +
theme(
  panel.spacing    = unit(2, "lines"),
  panel.border     = element_rect(linewidth = 0.5, fill = NA),
  axis.ticks       = element_line(linewidth = 0.4),
  axis.line        = element_blank(),
  strip.background = element_blank(),
  plot.title       = element_text(size = 16),
  plot.title.position = "plot",
  plot.subtitle    = element_text(size = 16),
  plot.caption     = element_text(size = 16),
  strip.text       = element_text(size = 16),
  legend.text      = element_text(size = 16),
  legend.title     = element_text(size = 16),
  axis.text        = element_text(size = 16),
  axis.title       = element_text(size = 16)
)
theme_set(theme_kamil)

# Show a progress bar
pboptions(type = "txt")

option_list <- list(
  make_option(
    c("-b", "--barcodes"), type = "character", default = NULL,
    help = "csv or tsv file with cell barcodes in the first column",
    metavar = "character"
  ),
  make_option(
    c("-f", "--file"), type = "character", default = NULL,
    help = "file name (molecule_info.h5,all_contig_annotations.csv,*.stat.csv.gz",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"), type = "character", default = NULL,
    help = "output folder, same as input file by default",
    metavar = "character"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

# opt$out <- "output"
# opt$barcodes <- "out2/barcodes.txt"
# opt$file <- "/projects/irae_blood/cellranger_output/C1_CD45A_gex/molecule_info.h5"
# opt$file <- "/projects/irae_blood/cellranger_output/C1_CD45A_tcr/all_contig_annotations.csv"
# opt$file <- "/projects/irae_blood/cellranger_output/irHepatitis_Batch_1A_ADT/irHepatitis_Batch_1A_ADT.stat.csv.gz"

if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
}
if (!file.exists(opt$file)) {
  status(glue("ERROR: File does not exist: '{opt$file}'"))
  quit(status = 1, save = "no")
}
file_base <- basename(opt$file)
file_dir <- basename(dirname(opt$file))
file_dir <- str_remove(file_dir, "^irHepatitis_")
accepted_files <- c(
  "GEX" = "^molecule_info.h5$",
  "VDJ" = "^all_contig_annotations.csv$",
  "ADT" = "^.+\\.stat\\.csv\\.gz$"
)
type <- names(accepted_files)[str_detect(file_base, accepted_files)]
if (length(type) == 0) {
  status(glue("ERROR: File is not valid: '{file_base}'"))
  status(glue("Accepted file formats:\n{paste(accepted_files, collapse = '\n')}"))
  quit(status = 1, save = "no")
}

barcodes <- NULL
if (!is.null(opt$barcodes)) {
  if (!file.exists(opt$barcodes)) {
    status(glue("ERROR: File does not exist: '{opt$barcodes}'"))
    quit(status = 1, save = "no")
  }
  barcodes <- unique(fread(opt$barcodes)[[1]])
  status(glue("INFO: Removing '-1' from the end of each barcode"))
  barcodes <- str_remove(barcodes, "-1$")
}

saturation <- function(reads, probs = seq(0.1, 1, by = 0.1), keep = NULL) {
  rbindlist(pblapply(probs, function(prob) {
    my_reads <- rbinom(n = length(reads), size = reads, prob = prob)
    if (!is.null(keep)) {
      my_reads <- my_reads[keep]
    }
    my_deduped_reads <- sum(my_reads > 0)
    my_total_reads <- sum(my_reads)
    list(
      prob = prob,
      sat = 1 - my_deduped_reads / my_total_reads,
      deduped_reads = my_deduped_reads,
      total_reads = my_total_reads
    )
  }))
}

if (type == "GEX") {

  status(glue("Reading {opt$file}"))
  gex_reads <- h5read(opt$file, "count")
  gex_barcodes <- h5read(opt$file, "barcode_idx")
  gex_barcodes <- h5read(opt$file, "barcodes")[gex_barcodes + 1]
  if (str_detect(gex_barcodes[1], "-1$")) {
    status(glue("INFO: Removing '-1' from the end of each barcode"))
    gex_barcodes <- str_remove(gex_barcodes, "-1$")
  }
  if (!is.null(barcodes)) {
    keep <- which(gex_barcodes %in% barcodes)
  } else {
    keep <- seq_len(length(gex_reads))
  }
  status(glue("Estimating GEX saturation for {length(unique(gex_barcodes[keep]))} barcodes"))
  sat <- saturation(gex_reads, keep = keep)
  total_sat <- 1 - sum(gex_reads[keep] > 0) / sum(gex_reads[keep])
  sat_file <- file.path(opt$out, "saturation-gex.tsv")
  fwrite(sat, sat_file, sep = "\t")

  my_caption <- glue("{signif(100 * total_sat, 2)}% saturation ({signif(mean(gex_reads[keep]), 2)} reads per UMI)")
  p <- ggplot(sat) +
    aes(x = total_reads, y = sat) +
    geom_hline(yintercept = total_sat, linewidth = 0.3, color = "red") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.3) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format()
    ) +
    scale_x_continuous(
      limits = c(0, max(sat$total_reads)),
      labels = label_number(scale_cut = cut_short_scale()),
      breaks = pretty_breaks(5)
    ) +
    labs(
      x = "Total reads",
      y = "Saturation",
      title = file_dir,
      caption = my_caption
    )
  p_file <- glue("{opt$out}/total_reads-vs-saturation-gex.pdf")
  status(glue("Writing {p_file}"))
  ggsave(p_file, p, width = 4, height = 3)

} else if (type == "VDJ") {

  status(glue("Reading {opt$file}"))
  d_vdj <- fread(opt$file)
  if (str_detect(d_vdj$barcode[1], "-1$")) {
    status(glue("INFO: Removing '-1' from the end of each barcode"))
    d_vdj$barcode <- str_remove(d_vdj$barcode, "-1$")
  }

  my_probs <- c(10 ^ seq(-6, -1, by = 0.5), seq(0.2, 1, length.out = 10))

  if (!is.null(barcodes)) {
    keep <- which(d_vdj$barcode %in% barcodes)
  } else {
    keep <- seq_len(nrow(d_vdj))
  }
  status(glue("Estimating VDJ saturation for {length(unique(d_vdj$barcode[keep]))} barcodes"))
  sat <- saturation(d_vdj$reads, my_probs, keep = keep)
  total_sat <- 1 - sum(d_vdj$reads[keep] > 0) / sum(d_vdj$reads[keep])
  sat_file <- file.path(opt$out, "saturation-vdj.tsv")
  fwrite(sat, sat_file, sep = "\t")

  my_caption <- glue("{signif(100 * total_sat, 2)}% saturation ({signif(mean(d_vdj$reads[keep]), 2)} reads per UMI)")
  p <- ggplot(sat) +
    aes(x = total_reads, y = sat) +
    geom_hline(yintercept = total_sat, linewidth = 0.3, color = "red") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.3) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format()
    ) +
    annotation_logticks(sides = "b", size = 0.3) +
    scale_x_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      trans = "log10"
    ) +
    labs(
      x = "Total reads",
      y = "Saturation",
      title = file_dir,
      caption = my_caption
    )
  p_file <- glue("{opt$out}/total_reads-vs-saturation-vdj.pdf")
  status(glue("Writing {p_file}"))
  ggsave(p_file, p, width = 4, height = 3)

  if (!is.null(barcodes)) {
    status("Estimating VDJ pct_cdr3 saturation")

    sat_cdr3 <- rbindlist(pblapply(my_probs, function(prob) {
      my_reads <- rbinom(n = nrow(d_vdj), size = d_vdj$reads, prob = prob)
      my_ix <- my_reads > 0 & with(d_vdj, productive & nchar(cdr3) > 0)
      n_cdr3 <- sum(barcodes %in% d_vdj$barcode[my_ix])
      list(
        prob = prob,
        pct_cdr3 = n_cdr3 / length(barcodes),
        total_reads = sum(d_vdj$reads) * prob
      )
    }))
    sat_cdr3_file <- file.path(opt$out, "saturation-pct_cdr3.tsv")
    fwrite(sat_cdr3, sat_cdr3_file, sep = "\t")

    p <- ggplot(sat_cdr3) +
      aes(x = total_reads, y = pct_cdr3) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.3) +
      scale_y_continuous(
        limits = c(0, 1),
        labels = percent_format()
      ) +
      annotation_logticks(sides = "b", size = 0.3) +
      scale_x_continuous(
        labels = label_number(scale_cut = cut_short_scale()),
        trans = "log10"
      ) +
      labs(
        x = "Total reads",
        y = "Percent",
        title = "Percent of cells with CDR3",
        subtitle = file_dir
      )
    p_file <- glue("{opt$out}/total_reads-vs-pct_cdr3.pdf")
    status(glue("Writing {p_file}"))
    ggsave(p_file, p, width = 4, height = 3)
  }

} else if (type == "ADT") {

  status(glue("Reading {opt$file}"))
  d_adt <- fread(opt$file)
  if (str_detect(d_adt$Barcode[1], "-1$")) {
    status(glue("INFO: Removing '-1' from the end of each barcode"))
    d_adt$Barcode <- str_remove(d_adt$Barcode, "-1$")
  }

  if (!is.null(barcodes)) {
    keep <- which(d_adt$Barcode %in% barcodes)
    status(glue("{length(unique(d_adt$Barcode[keep]))} ADT barcodes found in --barcodes"))
  } else {
    keep <- seq_len(nrow(d_adt))
  }

  if (length(keep) <= 1) {
    status("ERROR: ADT barcodes are not found in --barcodes")
    quit(status = 1, save = "no")
  }

  total_sat <- 1 - length(keep) / sum(d_adt$Count[keep])

  status(glue("Estimating ADT saturation"))
  sat <- saturation(d_adt$Count, keep = keep)
  sat_file <- file.path(opt$out, "saturation-adt.tsv")
  fwrite(sat, sat_file, sep = "\t")

  my_caption <- glue("{signif(100 * total_sat, 2)}% saturation ({signif(mean(d_adt$Count[keep]), 2)} reads per UMI)")
  p <- ggplot(sat) +
    aes(x = total_reads, y = sat) +
    geom_hline(yintercept = total_sat, linewidth = 0.3, color = "red") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.3) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format()
    ) +
    scale_x_continuous(
      limits = c(0, max(sat$total_reads)),
      labels = label_number(scale_cut = cut_short_scale()),
      breaks = pretty_breaks(5)
    ) +
    labs(
      x = "Total reads",
      y = "Saturation",
      title = file_dir,
      caption = my_caption
    )
  p_file <- glue("{opt$out}/total_reads-vs-saturation-adt.pdf")
  status(glue("Writing {p_file}"))
  ggsave(p_file, p, width = 4, height = 3)

  feature_sat <- d_adt[keep, .(saturation = 1 - .N / sum(Count)), by = c("Feature")]
  feature_sat_file <- file.path(opt$out, "saturation-adt-feature.tsv")
  fwrite(feature_sat, feature_sat_file, sep = "\t")
  p <- ggplot(feature_sat) +
    aes(x = saturation) +
    geom_histogram(bins = 71) +
    geom_vline(xintercept = total_sat, color = "red", linewidth = 0.3) +
    scale_x_continuous(
      labels = percent_format(),
      limits = c(-0.01, 1),
      expand = expansion(c(0.05, 0.08))
    ) + 
    annotate(
      geom = "text",
      x = total_sat, hjust = ifelse(total_sat < 0.5, -0.05, 1.05),
      y = Inf, vjust = 2,
      label = sprintf("%.1f%% (total)", 100 * total_sat)
    ) +
    labs(
      x = "PCR duplication rate (saturation)",
      y = "Count",
      title = file_dir,
      subtitle = "Histogram of saturation for each ADT feature"
    )
  p_file <- glue("{opt$out}/histogram-saturation-adt-feature.pdf")
  status(glue("Writing {p_file}"))
  ggsave(p_file, p, width = 5, height = 3.5)

}

