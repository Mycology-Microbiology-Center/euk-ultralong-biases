# euk-ultralong-biases

Code accompanying Tedersoo L, Flores-Almaraz V, Puusepp R, Mikryukov V. "Polymerase and primer biases in eukaryote ultra-long-read metabarcoding" (Molecular Ecology Resources, 2026).

Raw reads are available from the European Nucleotide Archive (ENA) under project [PRJEB123023](https://www.ebi.ac.uk/ena/browser/view/PRJEB123023). FASTQ files, reference databases, BLAST results, or generated results are not included in this repository.  


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
