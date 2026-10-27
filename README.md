<a href="https://mmc.ut.ee"><img src="assets/MycologyMicrobiologyCenter_logo.png" alt="Mycology and Microbiology Center" width="120" align="right"></a>

# Polymerase and primer biases in eukaryote ultra-long-read metabarcoding

Code accompanying Tedersoo L, Flores-Almaraz V, Puusepp R, Mikryukov V. "Polymerase and primer biases in eukaryote ultra-long-read metabarcoding" (Molecular Ecology Resources, 2026).

Raw reads are available from the European Nucleotide Archive (ENA) under project [PRJEB123023](https://www.ebi.ac.uk/ena/browser/view/PRJEB123023). FASTQ files, reference databases, BLAST results, or generated results are not included in this repository.  

## Input data

Demultiplexed reads can be downloaded from ENA with their submitted filenames using the first section of `01_nextits_processing.sh`. The ultra-long input contains all files except `ITS9mun_*`; the full-ITS input contains all files, because ultra-long reads were also trimmed to the ITS9mun-ITS4ngsUni region.

In the original analysis, NextITS was run on a larger set of demultiplexed files from the same sequencing runs. These included samples from an unrelated study and two negative controls, which together contained ~33% of ultra-long reads assigned to OTUs. The ENA submission contains only the samples analysed in the manuscript. Re-running the workflow from ENA data might therefore produce similar but not identical OTU sets, abundances, and index-switching scores. The OTU hashes of the five mock-community species in `05_mock_error_and_index_switching.R` may also need updating.

## Expected layout

```
data/
  reads/
    all/                       FASTQ files downloaded from ENA
    ultralong/                 symlinks to ultra-long FASTQ files
    full_its/                  symlinks to all FASTQ files
  blast/
    ultralong/                 tables with BLAST results for ITS, SSU, and LSU regions
    full_its/                  tables with BLAST results for ITS, SSU, and LSU regions
database/
  EUKARYOME_2.0_ITS.fasta
  EUKARYOME_2.0_ITS.udb
  RF01960.cm
results/                       generated locally; ignored by Git
```

BLAST tables are expected as serialized R objects named `BLAST_ITS.RData`, `BLAST_SSU.RData`, `BLAST_LSU.RData`, `BLAST_SSU_f500.RData`, and `BLAST_LSU_l500.RData`. Each table must contain `QueryName` and the ranked hit columns produced for the study. BLAST searches and manual taxonomic curation are outside the scope of this repository.

Input filenames must retain the primer-pair prefix and polymerase suffix used in the study, for example `EUK575F_<run>__<sample>_A.fq.gz`. Primer prefixes are `3NDL`, `EUK002F`, `EUK575F`, `Parabas`, and `ITS9mun`. Polymerase codes are A/F/G/J/K/N = Solis, B/H = TaKaRa, C/I = Platinum, D/L = Fusion, E/M = Phusion, and an absent suffix denotes Solis.


## Reference databases

EUKARYOME 2.0 reference databases are available from the [official archive](https://eukaryome.org/archive/). Download and extract the required references.

To download the covariance model from Rfam, you may use:
```bash
curl -fL https://rfam.org/family/RF01960/cm -o database/RF01960.cm
```

## Workflow

Scripts in this repository:

1. `01_nextits_processing.sh` - NextITS-based processing and OTU clustering  
2. `02_region_and_quality_processing.sh` - ITSx, UCHIME, Infernal, and sequence quality estimation commands  
3. `03_chimera_assessment.R` - per-sample UCHIME inputs  
4. `04_prepare_analysis_tables.R` - aggregate quality metrics and prepare analysis tables  
5. `05_mock_error_and_index_switching.R` - mock error and index switching assessment  

The workflow used Nextflow v25.10.4, NextITS v1.2.0, VSEARCH v2.31.0, SeqKit v2.13.0, ITSx v1.1.3, Infernal v1.1.5, phredsort v1.3.0 (https://github.com/vmikk/phredsort), MAFFT v7, and R v4.6.0. Required R packages are listed at the start of each R file.

Chimera and taxonomic assignments that required manual inspection remain described in the manuscript (the scripts prepare the relevant sequences and tables but do not replace those decisions).
