# euk-ultralong-biases

Code accompanying Tedersoo L, Flores-Almaraz V, Puusepp R, Mikryukov V. "Polymerase and primer biases in eukaryote ultra-long-read metabarcoding" (Molecular Ecology Resources, 2026).

Raw reads are available from the European Nucleotide Archive (ENA) under project [PRJEB123023](https://www.ebi.ac.uk/ena/browser/view/PRJEB123023). FASTQ files, reference databases, BLAST results, or generated results are not included in this repository.  


## Reference databases

EUKARYOME 2.0 reference databases are available from the [official archive](https://eukaryome.org/archive/). Download and extract the required references.

To download the covariance model from Rfam, you may use:
```bash
curl -fL https://rfam.org/family/RF01960/cm -o database/RF01960.cm
```

