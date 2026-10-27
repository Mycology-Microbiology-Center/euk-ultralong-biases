# Combine NextITS results, sequence-quality metrics, chimera scores, and BLAST tables

library(data.table)
library(arrow)
library(dplyr)
library(Biostrings)
library(openxlsx)

datasets <- c("ultralong", "full_its")


###############################################################################
# Treatment metadata encoded in the demultiplexed-read filenames
###############################################################################

add_treatment_metadata <- function(table) {
  table <- copy(table)
  table[, RunID := tstrsplit(SampleID, "__", fixed = TRUE, keep = 1L)]

  table[, PrimerPair := fcase(
    grepl("^3NDL_",    SampleID), "3NDF-21R",
    grepl("^EUK002F_", SampleID), "Euk002F-Euk5493R",
    grepl("^EUK575F_", SampleID), "Euk575F-21Rngs",
    grepl("^Parabas_", SampleID), "Euk575F-21Rparabas",
    grepl("^ITS9mun_", SampleID), "ITS9mun-ITS4ngsUni",
    default = NA_character_
  )]

  table[, Polymerase := fcase(
    grepl("_[AFGJKN]$", SampleID), "Solis",
    grepl("_[BH]$", SampleID),     "TaKaRa",
    grepl("_[CI]$", SampleID),     "Platinum",
    grepl("_[DL]$", SampleID),     "Fusion",
    grepl("_[EM]$", SampleID),     "Phusion",
    default = "Solis"
  )]

  table
}


###############################################################################
# Input helpers
###############################################################################

read_quality_data <- function(dataset, otu_ids) {
  step1_dir <- file.path("results", "nextits", paste0("step1_", dataset))
  quality_files <- list.files(
    step1_dir,
    pattern = "SeqQualities[.]parquet$",
    recursive = TRUE,
    full.names = TRUE)

  if (!length(quality_files)) {
    stop("No SeqQualities.parquet files found under ", step1_dir)
  }

  open_dataset(quality_files, format = "parquet") |>
    filter(Hash %in% otu_ids) |>
    select(Hash, AvgPhredScore, MaxEE, MEEP) |>
    collect() |>
    as.data.table()
}

read_blast_tables <- function(dataset) {
  input_dir <- file.path("data", "blast", dataset)
  filenames <- c(
    ITS = "BLAST_ITS.RData",
    SSU = "BLAST_SSU.RData",
    LSU = "BLAST_LSU.RData",
    SSU_first_500 = "BLAST_SSU_f500.RData",
    LSU_last_500 = "BLAST_LSU_l500.RData")

  paths <- setNames(file.path(input_dir, filenames), names(filenames))
  if (any(!file.exists(paths))) {
    stop(
      "Missing precomputed BLAST tables: ",
      paste(paths[!file.exists(paths)], collapse = ", ")
    )
  }

  lapply(paths, function(path) {
    table <- readRDS(path)
    setDT(table)
    table
  })
}

align_blast_table <- function(table, otu_ids) {
  if (anyDuplicated(table$QueryName)) {
    stop("BLAST table contains duplicated QueryName values")
  }

  index <- match(otu_ids, table$QueryName)
  aligned <- table[index]
  missing <- is.na(index)
  aligned[missing, QueryName := otu_ids[missing]]
  aligned
}


###############################################################################
# Compile one dataset
###############################################################################

compile_dataset <- function(dataset, chimera_results) {
  step2_dir <- file.path(
    "results", "nextits", paste0("step2_", dataset), "04.PooledResults"
  )

  otu_table <- readRDS(file.path(step2_dir, "OTU_table_long.RData"))
  setDT(otu_table)
  otu_table <- add_treatment_metadata(otu_table)

  otu_summary <- otu_table[, .(
    Occurrence = uniqueN(SampleID),
    TotalAbundance = sum(Abundance, na.rm = TRUE)
  ), by = OTU]
  setorder(otu_summary, -TotalAbundance, -Occurrence, OTU)

  representative_sequences <- readDNAStringSet(file.path(step2_dir, "OTUs.fa.gz"))
  representative_table <- data.table(
    OTU = sub(";.*$", "", names(representative_sequences)),
    RepresentativeSequence = as.character(representative_sequences)
  )
  representative_table[, SequenceLength := nchar(RepresentativeSequence)]

  otu_summary <- merge(
    otu_summary,
    representative_table,
    by = "OTU",
    all.x = TRUE,
    sort = FALSE)

  # Quality of the representative sequence:
  # the highest-quality read among the identical sequences of the OTU representative
  quality_data <- read_quality_data(dataset, otu_summary$OTU)
  setorder(quality_data, Hash, -AvgPhredScore)
  quality_data <- unique(quality_data, by = "Hash")
  quality_summary <- quality_data[, .(
    OTU = Hash,
    AvgPhredScore,
    MaxEE,
    MEEP
  )]

  otu_summary <- merge(
    otu_summary,
    quality_summary,
    by = "OTU",
    all.x = TRUE,
    sort = FALSE)

  dataset_chimeras <- chimera_results[[dataset]]$summary
  otu_summary <- merge(
    otu_summary,
    dataset_chimeras,
    by = "OTU",
    all.x = TRUE,
    sort = FALSE)

  otu_summary[, Incomplete_lt2600 := SequenceLength < 2600]
  otu_summary[, LowQualityCandidate :=
    TotalAbundance < 10 & AvgPhredScore < 25]
  otu_summary[, SingletonOrDoubleton := TotalAbundance < 3]
  otu_summary[, DeNovoScore_gt1 := DeNovoChimera_ScoreMax > 1]
  otu_summary[, DeNovoScore_gt2 := DeNovoChimera_ScoreMax > 2]
  otu_summary[, DeNovoScore_gt10 := DeNovoChimera_ScoreMax > 10]

  otu_summary[, `:=`(
    ManualTaxonomy = NA_character_,
    ManualQualityClass = NA_character_,
    ManualNotes = NA_character_
  )]
  setorder(otu_summary, -TotalAbundance, -Occurrence, OTU)

  otu_wide <- dcast(
    otu_table,
    OTU ~ SampleID,
    value.var = "Abundance",
    fun.aggregate = sum,
    fill = 0)
  otu_wide <- otu_wide[match(otu_summary$OTU, OTU)]

  blast_tables <- read_blast_tables(dataset)
  blast_tables <- lapply(
    blast_tables,
    align_blast_table,
    otu_ids = otu_summary$OTU)

  c(
    list(
      OTU_summary = otu_summary,
      OTU_table_long = otu_table,
      OTU_table_wide = otu_wide
    ),
    blast_tables
  )
}


###############################################################################
# Compile and export the two analysis datasets
###############################################################################

chimera_results <- readRDS(
  file.path("results", "tables", "chimera_results.rds")
)

dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

compiled <- setNames(lapply(datasets, function(dataset) {
  result <- compile_dataset(dataset, chimera_results)

  saveRDS(
    result,
    file.path("results", "tables", paste0(dataset, "_analysis_tables.rds")))
  
  write.xlsx(
    result,
    file.path("results", "tables", paste0(dataset, "_analysis_tables.xlsx")),
    colNames = TRUE,
    overwrite = TRUE)

  result
}), datasets)
