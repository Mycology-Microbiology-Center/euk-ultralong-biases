# Prepare mock-community SSU reads for error assessment and estimate index-switching

library(data.table)
library(arrow)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(openxlsx)

selected_otus <- c(
  "545e673b15cfd513a0d0f087939015dcaf3f0ebc",
  "8d18beb8591353f51999af1b651861e12dd9f1ac",
  "48cc2e65961f0fe9ded3b15f07943a918c51e6b3",
  "4e2dd4dde821b7b4c67b945cc3e3e01fc8524ccc",
  "f5b035a3f2eb5f74a83abda01474894ff8fd8a3f"
)

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
# Extract all mock-community reads assigned to the five selected OTUs
###############################################################################

step1_dir <- file.path("results", "nextits", "step1_ultralong")
step2_dir <- file.path(
  "results", "nextits", "step2_ultralong", "04.PooledResults"
)

uc_table <- read_parquet(file.path(step2_dir, "UC_Pooled.parquet")) |>
  as.data.table()
uc_table <- uc_table[OTU %in% selected_otus]

quality_files <- list.files(
  step1_dir,
  pattern = "SeqQualities[.]parquet$",
  recursive  = TRUE,
  full.names = TRUE)

if (!length(quality_files)) {
  stop("No SeqQualities.parquet files found under ", step1_dir)
}

mock_reads <- open_dataset(quality_files, format = "parquet") |>
  filter(Hash %in% uc_table$SeqID) |>
  collect() |>
  as.data.table()

mock_reads <- mock_reads[grepl("mock", SampleID, ignore.case = TRUE)]
mock_reads <- merge(
  mock_reads,
  uc_table[, .(SeqID, OTU)],
  by.x  = "Hash",
  by.y  = "SeqID",
  all.x = TRUE,
  sort  = FALSE)

mock_reads <- add_treatment_metadata(mock_reads)

dir.create(file.path("results", "mock"), recursive = TRUE, showWarnings = FALSE)
saveRDS(mock_reads, file.path("results", "mock", "mock_sequences.rds"))

unique_reads <- unique(mock_reads, by = "Sequence")
full_sequences <- DNAStringSet(unique_reads$Sequence)
names(full_sequences) <- unique_reads$Hash

writeXStringSet(
  full_sequences,
  filepath = file.path("results", "mock", "mock_sequences_full.fasta.gz"),
  compress = TRUE,
  format   = "fasta",
  width    = 9999)


###############################################################################
# Parse Infernal output and extract the Euk575F-to-SSU-terminus segment
# Run after the cmsearch section in 02_region_and_quality_processing.sh
###############################################################################

parse_cmsearch_tblout <- function(file) {
  column_names <- c(
    "target_name", "target_accession", "query_name", "query_accession",
    "model", "model_from", "model_to", "seq_from", "seq_to", "strand",
    "trunc", "pass", "gc", "bias", "score", "E_value", "inc", "description")

  lines <- readLines(file)
  lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]

  # The last column (description) may contain spaces
  rows <- lapply(strsplit(trimws(lines), "[[:space:]]+"), function(fields) {
    c(fields[1:17], paste(fields[-(1:17)], collapse = " "))
  })

  result <- as.data.table(do.call(rbind, rows))
  setnames(result, column_names)

  integer_columns <- c("model_from", "model_to", "seq_from", "seq_to", "pass")
  numeric_columns <- c("gc", "bias", "score", "E_value")
  result[, (integer_columns) := lapply(.SD, as.integer), .SDcols = integer_columns]
  result[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
  result
}

infernal <- parse_cmsearch_tblout( file.path("results", "mock", "Infernal_SSU.tblout") )
infernal <- infernal[inc == "!" & strand == "+"]
setorder(infernal, target_name, -score)
infernal <- unique(infernal, by = "target_name")

mock_reads <- readRDS(file.path("results", "mock", "mock_sequences.rds"))
setDT(mock_reads)
mock_reads <- merge(
  mock_reads,
  infernal[, .(target_name, seq_to)],
  by.x = "Hash",
  by.y = "target_name",
  all  = FALSE,
  sort = FALSE)

mock_reads[, SSU_sequence := substr(Sequence, 1L, seq_to)]
mock_reads[, SSU_quality  := substr(Quality, 1L, seq_to) ]

ssu_fastq <- ShortReadQ(
  sread = DNAStringSet(mock_reads$SSU_sequence),
  quality = FastqQuality(mock_reads$SSU_quality),
  id = BStringSet(mock_reads$PacBioID))

writeFastq(
  ssu_fastq,
  file = file.path("results", "mock", "mock_sequences_SSU.fastq.gz"),
  mode = "w",
  full = FALSE,
  compress = TRUE)

saveRDS(mock_reads, file.path("results", "mock", "mock_sequences_SSU.rds"))


###############################################################################
# Add SSU Phred scores and retain up to 100 maximum-quality reads per treatment
# Run after the phredsort section in 02_region_and_quality_processing.sh
###############################################################################

mock_reads <- readRDS(
  file.path("results", "mock", "mock_sequences_SSU.rds")
)
setDT(mock_reads)

ssu_quality <- fread(
  file.path("results", "mock", "mock_sequences_SSU_quality.tsv"),
  sep = "\t",
  header = FALSE,
  col.names = c("PacBioID", "SSU_AvgPhredScore"))

mock_reads <- merge(
  mock_reads,
  ssu_quality,
  by = "PacBioID",
  all.x = TRUE,
  sort = FALSE)

stopifnot(!anyNA(mock_reads$SSU_AvgPhredScore))

ssu_summary <- mock_reads[, .(
  NumReads = .N,
  NumUniqueFullSequences = uniqueN(Sequence),
  FullSequencePhred_Mean = mean(AvgPhredScore),
  FullSequencePhred_Min  = min(AvgPhredScore),
  FullSequencePhred_Max  = max(AvgPhredScore),
  SSUPhred_Mean          = mean(SSU_AvgPhredScore),
  SSUPhred_Min           = min(SSU_AvgPhredScore),
  SSUPhred_Max           = max(SSU_AvgPhredScore),
  SSU_Length             = nchar(SSU_sequence[1L])
), by = .(OTU, PrimerPair, Polymerase, SampleID, SSU_sequence)]

set.seed(1L)
alignment_reads <- mock_reads[SSU_AvgPhredScore == 40]
alignment_reads <- alignment_reads[, {
  retained <- sample.int(.N, min(100L, .N))
  .SD[retained]
}, by = .(OTU, PrimerPair, Polymerase)]

alignment_dir <- file.path("results", "mock", "alignment_input")
dir.create(alignment_dir, recursive = TRUE, showWarnings = FALSE)

# One alignment per species (OTU); treatment is encoded in the read names
alignment_groups <- split(alignment_reads, by = "OTU", keep.by = TRUE)
for (group in alignment_groups) {
  sequences <- DNAStringSet(group$SSU_sequence)
  names(sequences) <- paste(
    group$PacBioID, group$PrimerPair, group$Polymerase,
    sep = "|")
  writeXStringSet(
    sequences,
    filepath = file.path(alignment_dir, paste0(group$OTU[1L], ".fasta")),
    format = "fasta",
    width = 9999)
}

write.xlsx(
  list(SSU_summary = ssu_summary),
  file.path("results", "mock", "mock_SSU_sequence_summary.xlsx"),
  colNames  = TRUE,
  overwrite = TRUE)


###############################################################################
# UNCROSS-style index-switch estimate recorded in the working analysis
###############################################################################

uncross_score <- function(abundance, total_abundance,
  number_of_samples,
  expected_rate = 0.01, minimum_score = 0.1,
  exponent = 1) {

  expected_abundance <- expected_rate * total_abundance / number_of_samples
  score <- 2 / (1 + exp(abundance / expected_abundance)^exponent)

  data.table(
    IndexSwitchScore = score,
    IndexSwitch = score >= minimum_score
  )
}

otu_table <- readRDS(
  file.path(
    "results", "nextits", "step2_ultralong", "04.PooledResults",
    "OTU_table_long.RData"))

setDT(otu_table)
otu_table <- add_treatment_metadata(otu_table)
otu_table[, TotalAbundance := sum(Abundance, na.rm = TRUE), by = OTU]

index_scores <- uncross_score(
  abundance = otu_table$Abundance,
  total_abundance = otu_table$TotalAbundance,
  number_of_samples = uniqueN(otu_table$SampleID),
  expected_rate = 0.01,
  minimum_score = 0.1,
  exponent = 1
)
otu_table <- cbind(otu_table, index_scores)

# This treatment contained only one sample and was excluded from summaries.
otu_table <- otu_table[RunID != "Parabas_LONG2"]

index_summary <- otu_table[, .(
  NumOTUs                   = uniqueN(OTU),
  TotalReads_Unfiltered     = sum(Abundance, na.rm = TRUE),
  TotalReads_Filtered       = sum(Abundance[!IndexSwitch], na.rm = TRUE),
  NumIndexSwitchOccurrences = sum(IndexSwitch, na.rm = TRUE),
  NumIndexSwitchReads       = sum(Abundance[IndexSwitch], na.rm = TRUE)
), by = .(PrimerPair, Polymerase, SampleID)]

dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(per_occurrence = otu_table, per_sample = index_summary),
  file.path("results", "tables", "index_switching.rds"))

write.xlsx(
  list(Index_switch_summary = index_summary),
  file.path("results", "tables", "index_switching.xlsx"),
  colNames = TRUE,
  overwrite = TRUE)
