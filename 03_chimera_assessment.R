# Prepare and aggregate chimera assessments


library(data.table)
library(Biostrings)

datasets <- c("ultralong", "full_its")


###############################################################################
# Prepare abundance-annotated OTU representatives for per-sample UCHIME runs
###############################################################################

prepare_denovo_input <- function(dataset) {
  step2_dir <- file.path("results", "nextits", paste0("step2_", dataset),
                         "04.PooledResults")
  otu_table <- readRDS(file.path(step2_dir, "OTU_table_long.RData"))
  setDT(otu_table)

  sequences <- readDNAStringSet(file.path(step2_dir, "OTUs.fa.gz"))
  sequence_table <- data.table(
    OTU = sub(";.*$", "", names(sequences)),
    Sequence = as.character(sequences))

  output_dir <- file.path("results", "chimeras", "input", dataset)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_tables <- split(otu_table, by = "SampleID", keep.by = TRUE)
  for (sample_id in names(sample_tables)) {
    sample_table <- merge(
      sample_tables[[sample_id]], sequence_table,
      by = "OTU", all.x = TRUE, sort = FALSE)
    stopifnot(!anyNA(sample_table$Sequence))
    setorder(sample_table, -Abundance, OTU)

    sample_sequences <- DNAStringSet(sample_table$Sequence)
    names(sample_sequences) <- paste0(
      sample_table$OTU, ";size=", sample_table$Abundance
    )

    writeXStringSet(
      sample_sequences,
      filepath = file.path(output_dir, paste0(sample_id, ".fasta.gz")),
      compress = TRUE,
      format = "fasta",
      width = 19999)
  }

  invisible(output_dir)
}

invisible(lapply(datasets, prepare_denovo_input))


###############################################################################
# Aggregate VSEARCH de novo and reference-based chimera results
# Run after the UCHIME sections in 02_region_and_quality_processing.sh
###############################################################################

read_denovo_results <- function(dataset) {
  input_dir <- file.path("results", "chimeras", "denovo", dataset)
  files <- list.files(input_dir, pattern = "[.]tsv$", full.names = TRUE)
  files <- files[file.info(files)$size > 0]

  if (!length(files)) {
    return(list(
      per_sample = data.table(
        OTU = character(),
        DeNovoChimera_Score = numeric(),
        SampleID = character()
      ),
      summary = data.table(
        OTU = character(),
        DeNovoChimera_N_samples = integer(),
        DeNovoChimera_ScoreMax = numeric(),
        DeNovoChimera_ScoreMean = numeric()
      )
    ))
  }

  per_sample <- rbindlist(lapply(files, function(file) {
    fread(
      file,
      sep = "\t",
      header = FALSE,
      col.names = c("OTU", "DeNovoChimera_Score", "SampleID")
    )
  }), use.names = TRUE)

  summary <- per_sample[, .(
    DeNovoChimera_N_samples = uniqueN(SampleID),
    DeNovoChimera_ScoreMax = max(DeNovoChimera_Score, na.rm = TRUE),
    DeNovoChimera_ScoreMean = mean(DeNovoChimera_Score, na.rm = TRUE)
  ), by = OTU]

  list(per_sample = per_sample, summary = summary)
}

read_reference_fasta <- function(file, status) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    return(data.table(
      OTU = character(),
      RefChimera_Status = character(),
      RefChimera_Score = numeric()
    ))
  }

  sequences <- readDNAStringSet(file)
  headers <- names(sequences)
  scores <- sub("^.*uchime_ref=([^;]+).*$", "\\1", headers)

  data.table(
    OTU = sub(";.*$", "", headers),
    RefChimera_Status = status,
    RefChimera_Score = as.numeric(scores))
}

read_reference_results <- function(dataset) {
  input_dir <- file.path("results", "chimeras", "reference", dataset)
  rbindlist(list(
    read_reference_fasta(
      file.path(input_dir, "chimeras.fasta"), "chimera"
    ),
    read_reference_fasta(
      file.path(input_dir, "borderline.fasta"), "borderline"
    )
  ), use.names = TRUE, fill = TRUE)
}

chimera_results <- setNames(lapply(datasets, function(dataset) {
  denovo <- read_denovo_results(dataset)
  reference <- read_reference_results(dataset)

  combined <- merge(
    denovo$summary,
    reference,
    by = "OTU",
    all = TRUE,
    sort = FALSE)

  list(
    de_novo_per_sample = denovo$per_sample,
    summary = combined)
}), datasets)

dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
saveRDS(chimera_results, file.path("results", "tables", "chimera_results.rds"))
