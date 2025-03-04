# adapted from https://github.com/umccr/sigrap/blob/main/R/mutationalpatterns.R
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(assertthat)

# functions
sig_contribution <- function(mut_mat, signatures) {
  # Fit mutation matrix to cancer signatures
  fit_res <-
    MutationalPatterns::fit_to_signatures(mut_mat, signatures)$contribution |>
    tibble::as_tibble(rownames = "sig") |>
    dplyr::rename(contr = 2) |>
    dplyr::filter(.data$contr > 0)
  
  if (nrow(fit_res) == 0) {
    fit_res <- tibble::tribble(
      ~sig, ~contr,
      "No Signatures found!", 0
    )
  }
  
  fit_res_contr <- fit_res |>
    dplyr::mutate(
      contr = round(.data$contr, 0),
      RelFreq = round(.data$contr / sum(.data$contr), 2),
      Rank = as.integer(base::rank(-.data$contr))
    ) |>
    dplyr::select("Rank",
                  Signature = "sig",
                  Contribution = "contr", "RelFreq"
    ) |>
    dplyr::arrange(.data$Rank)
  
  fit_res_contr
}

sig_contribution_table <- function(contr, type, outdir = NULL) {
  available_types <- c(
    "Sig" = "v2_2015/Sig",
    "SBS" = "v3.2_2021-march/SBS",
    "DBS" = "v3.2_2021-march/DBS",
    "ID" = "v3.2_2021-march/ID"
  )
  assertthat::assert_that(length(type) == 1, type %in% names(available_types))
  assertthat::assert_that(all(colnames(contr) == c("Rank", "Signature", "Contribution", "RelFreq")))
  
  #sig_dir <- system.file(file.path("extdata/sigs", available_types[type]), package = "sigrap")
  sig_dir <- '~/Downloads/'
  outdir <- '~/Desktop/signatures'
  # if an outdir is specified, copy out the COSMIC plots to be linked to
  if (!is.null(outdir)) {
    img_cp_dir <- file.path(outdir, "sig_plots", type)
    fs::dir_create(img_cp_dir)
    sig_table <-
      readr::read_tsv(file = file.path(sig_dir, "description.tsv"), col_types = "cc") |>
      dplyr::mutate(
        Plot_original = file.path(sig_dir, paste0(.data$signature, ".png")),
        Plot_copy = file.path(img_cp_dir, paste0(.data$signature, ".png")),
        signature = paste0(type, .data$signature)
      ) |>
      dplyr::rename(Signature = .data$signature)
    
    d <- contr |>
      dplyr::left_join(sig_table, by = "Signature") |>
      dplyr::rowwise() |>
      dplyr::mutate(cp = file.copy(from = .data$Plot_original, to = .data$Plot_copy)) |>
      dplyr::mutate(plot_in_md = paste0("![](", .data$Plot_copy, ")")) |>
      dplyr::select(
        "Rank", "Signature", "Contribution", "RelFreq",
        Description = "description", Plot = "plot_in_md"
      )
  } else { # don't copy out any COSMIC plots
    sig_table <-
      readr::read_tsv(file = file.path(sig_dir, "description.tsv"), col_types = "cc") |>
      dplyr::mutate(
        signature = paste0(type, .data$signature)
      ) |>
      dplyr::rename(Signature = .data$signature)
    d <- contr |>
      dplyr::left_join(sig_table, by = "Signature") |>
      dplyr::rename(Description = .data$description)
  }
  return(d)
}

sig_count_snv <- function(vcf_gr, ref_genome) {
  gr_snv <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "snv")
  snv_counts <- MutationalPatterns::mut_matrix(vcf_list = gr_snv, ref_genome = ref_genome)
  list(
    snv_counts = snv_counts,
    gr_snv = gr_snv
  )
}

sig_plot_snv <- function(gr_snv, snv_counts, ref_genome) {
  mut_to <- MutationalPatterns::mut_type_occurrences(
    vcf_list = gr_snv, ref_genome = ref_genome
  )
  
  mut_mat_ext_context <- MutationalPatterns::mut_matrix(
    vcf_list = gr_snv, ref_genome = ref_genome, extension = 2
  )
  
  p_spectrum <- MutationalPatterns::plot_spectrum(
    type_occurrences = mut_to, CT = TRUE,
    condensed = TRUE, error_bars = "none"
  ) +
    ggplot2::theme(legend.position = "top")
  p_96_profile <- MutationalPatterns::plot_96_profile(mut_matrix = snv_counts, condensed = TRUE)
  p_heatmap <- MutationalPatterns::plot_profile_heatmap(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")
  p_river <- MutationalPatterns::plot_river(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")
  
  list(
    p_heatmap = p_heatmap,
    p_river = p_river,
    p_96_profile = p_96_profile,
    p_spectrum = p_spectrum
  )
}

sig_count_indel <- function(vcf_gr, ref_genome) {
  gr_indel <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "indel")
  gr_indel <- MutationalPatterns::get_indel_context(vcf_list = gr_indel, ref_genome = ref_genome)
  indel_counts <- MutationalPatterns::count_indel_contexts(vcf_list = gr_indel)
  indel_counts
}

sig_plot_indel <- function(indel_counts) {
  p_indel_main <- MutationalPatterns::plot_main_indel_contexts(counts = indel_counts) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p_indel_cont <- MutationalPatterns::plot_indel_contexts(counts = indel_counts, condensed = TRUE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1),
      legend.position = "top"
    )
  
  list(
    p_indel_main = p_indel_main,
    p_indel_cont = p_indel_cont
  )
}

sig_count_dbs <- function(vcf_gr) {
  # dbs need additional editing due to https://github.com/UMCUGenetics/MutationalPatterns/issues/56
  gr = vcf_gr[[1]]
  dbs_i = elementNROWS(gr$REF) == 2 & sapply(gr$ALT, function(x) nchar(x) == 2)
  dbs_predetermined_gr = gr[dbs_i]
  other_gr = gr[!dbs_i]
  gr_dbs = MutationalPatterns::get_mut_type(other_gr, type = "dbs")
  #dbs_gr = c(dbs_merged_gr, dbs_predetermined_gr)
  #gr_dbs <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "dbs")
  gr_dbs <- MutationalPatterns::get_dbs_context(vcf_list = gr_dbs)
  dbs_counts <- MutationalPatterns::count_dbs_contexts(vcf_list = gr_dbs)
  dbs_counts
}

sig_plot_dbs <- function(dbs_counts) {
  p_dbs_main <- MutationalPatterns::plot_main_dbs_contexts(counts = dbs_counts)
  p_dbs_cont <- MutationalPatterns::plot_dbs_contexts(counts = dbs_counts, condensed = TRUE)
  
  list(
    p_dbs_main = p_dbs_main,
  )
}

save_plot_list <- function(pl, outdir) {
  assertthat::assert_that(is.list(pl), length(pl) > 0)
  assertthat::assert_that(all(sapply(pl, inherits, "ggplot")))
  fs::dir_create(outdir)
  for (i in seq_len(length(pl))) {
    nm <- names(pl)[i]
    fn <- file.path(outdir, paste0(nm, ".png"))
    plot_obj <- pl[[i]]
    ggplot2::ggsave(filename = fn, plot = plot_obj)
  }
}

write_jsongz <- function(x, path, ...) {
  assertthat::assert_that(endsWith(path, ".gz"))
  fs::dir_create(dirname(path))
  gz <- gzfile(path, open = "w")
  jsonlite::write_json(x = x, path = gz, ...)
  close(gz)
}


vcf_file <- "/Users/kanwals/UMCCR/data/projects/duplex-seq/01SkNC1D9/01SkNC1D9.1.consensus.variant-calls.vcf.gz"
sample_name <- "01SkNC1D9.1"

# read vcf
grl <- MutationalPatterns::read_vcfs_as_granges(vcf_file, 
                                                sample_name, 
                                                ref_genome, 
                                                group = "auto+sex",
                                                type = "all")

#---- SBS ----#
# plots
snv_counts <- sig_count_snv(vcf_gr = grl, ref_genome = ref_genome)
p_snv <- sig_plot_snv(
  gr_snv = snv_counts$gr_snv, snv_counts = snv_counts$snv_counts,
  ref_genome = ref_genome
)

# signature contributions (2020)
sigs_snv_2020 <-
  MutationalPatterns::get_known_signatures(
    muttype = "snv",
    incl_poss_artifacts = TRUE
  ) |>
  {
    \(sigs) sig_contribution(mut_mat = snv_counts$snv_counts, signatures = sigs)
  }() |>
  sig_contribution_table(type = "SBS")

#---- DBS ----#
# plots
dbs_counts <- sig_count_dbs(vcf_gr = grl)
p_dbs <- sig_plot_dbs(dbs_counts = dbs_counts)

# signature contributions
sigs_dbs <-
  MutationalPatterns::get_known_signatures(muttype = "dbs") |>
  {
    \(sigs) sig_contribution(mut_mat = dbs_counts, signatures = sigs)
  }() |>
  sig_contribution_table(type = "DBS")

#---- Indels ----#
# plots
indel_counts <- sig_count_indel(vcf_gr = gr, ref_genome = ref_genome)
p_indel <- sig_plot_indel(indel_counts = indel_counts)

# signature contributions
sigs_indel <-
  MutationalPatterns::get_known_signatures(muttype = "indel") |>
  {
    \(sigs) sig_contribution(mut_mat = indel_counts, signatures = sigs)
  }() |>
  sig_contribution_table(type = "ID")

# save results
save_plot_list(p_snv, file.path(outdir, "plot/snv"))
save_plot_list(p_dbs, file.path(outdir, "plot/dbs"))
save_plot_list(p_indel, file.path(outdir, "plot/indel"))
write_jsongz(x = sigs_snv_2020, path = file.path(outdir, "sigs/snv2020.json.gz"))
write_jsongz(x = sigs_dbs, path = file.path(outdir, "sigs/dbs.json.gz"))
write_jsongz(x = sigs_indel, path = file.path(outdir, "sigs/indel.json.gz"))
