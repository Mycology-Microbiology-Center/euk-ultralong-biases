# Region extraction and quality assessment
# (precomputed BLAST tables are expected under data/blast/)

threads=8

mkdir -p results/regions/ultralong results/regions/full_its
mkdir -p results/chimeras results/mock


###############################################################################
# Extract ITS, SSU, and LSU regions from 99% OTU representatives with ITSx 1.1.3
###############################################################################

seqkit seq -w 0 \
  results/nextits/step2_ultralong/04.PooledResults/OTUs.fa.gz \
  > results/regions/ultralong/OTUs.fasta

ITSx \
  -i results/regions/ultralong/OTUs.fasta \
  -o results/regions/ultralong/ultralong \
  --complement T \
  --save_regions all \
  --graphical F \
  --detailed_results T \
  --positions T \
  --not_found T \
  -E 1e-1 \
  -t all \
  --partial 0 \
  --cpu "$threads" \
  --preserve T

seqkit seq -w 0 \
  results/nextits/step2_full_its/04.PooledResults/OTUs.fa.gz \
  > results/regions/full_its/OTUs.fasta

ITSx \
  -i results/regions/full_its/OTUs.fasta \
  -o results/regions/full_its/full_its \
  --complement T \
  --save_regions all \
  --graphical F \
  --detailed_results T \
  --positions T \
  --not_found T \
  -E 1e-1 \
  -t all \
  --partial 0 \
  --cpu "$threads" \
  --preserve T


###############################################################################
# Prepare terminal SSU and LSU fragments used for taxonomic/chimera review
###############################################################################

seqkit subseq -r 1:500 -w 0 \
  results/regions/ultralong/ultralong.SSU.fasta \
  > results/regions/ultralong/SSU_first_500.fasta

seqkit subseq -r -500:-1 -w 0 \
  results/regions/ultralong/ultralong.LSU.fasta \
  > results/regions/ultralong/LSU_last_500.fasta

seqkit subseq -r 1:500 -w 0 \
  results/regions/full_its/full_its.SSU.fasta \
  > results/regions/full_its/SSU_first_500.fasta

seqkit subseq -r -500:-1 -w 0 \
  results/regions/full_its/full_its.LSU.fasta \
  > results/regions/full_its/LSU_last_500.fasta


###############################################################################
# De novo chimera detection with VSEARCH 2.31.0
# First run the input-preparation section in 03_chimera_assessment.R
###############################################################################

detect_denovo_chimeras() {
  input_file="$1"
  output_file="$2"
  sample_id=$(basename "$input_file" .fasta.gz)

  vsearch \
    --uchime_denovo "$input_file" \
    --abskew 2.0 \
    --dn 1.4 \
    --mindiffs 3 \
    --mindiv 0.8 \
    --minh 0.28 \
    --xn 8.0 \
    --threads 1 \
    --qmask dust \
    --sizein \
    --xsize \
    --fasta_width 0 \
    --fasta_score \
    --chimeras - \
  | seqkit seq --name \
  | awk -v sample="$sample_id" 'BEGIN {OFS="\t"} {gsub(/;\+/, ";"); score=$0; sub(/^.*uchime_denovo=/, "", score); sub(/;.*/, "", score); otu=$0; sub(/;.*/, "", otu); print otu, score, sample}' \
  > "$output_file"
}

for dataset in ultralong full_its; do
  mkdir -p "results/chimeras/denovo/${dataset}"
  while IFS= read -r input_file; do
    sample_id=$(basename "$input_file" .fasta.gz)
    detect_denovo_chimeras \
      "$input_file" \
      "results/chimeras/denovo/${dataset}/${sample_id}.tsv"
  done < <(find "results/chimeras/input/${dataset}" -name '*.fasta.gz' -type f | sort)
done


###############################################################################
# Reference-based chimera detection against EUKARYOME 2.0 ITS
###############################################################################

vsearch \
  --makeudb_usearch database/EUKARYOME_2.0_ITS.fasta \
  --output database/EUKARYOME_2.0_ITS.udb

for dataset in ultralong full_its; do
  mkdir -p "results/chimeras/reference/${dataset}"

  vsearch \
    --uchime_ref "results/nextits/step2_${dataset}/04.PooledResults/OTUs.fa.gz" \
    --db database/EUKARYOME_2.0_ITS.udb \
    --selfid \
    --fasta_width 0 \
    --threads "$threads" \
    --sizein \
    --sizeout \
    --fasta_score \
    --chimeras "results/chimeras/reference/${dataset}/chimeras.fasta" \
    --borderline "results/chimeras/reference/${dataset}/borderline.fasta"
done


###############################################################################
# Delimit SSU in selected mock-community reads with Infernal 1.1.5
# First run the extraction section in 05_mock_error_and_index_switching.R
###############################################################################

seqkit seq -w 0 \
  results/mock/mock_sequences_full.fasta.gz \
  > results/mock/mock_sequences_full.fasta

cmsearch -g \
  --tblout results/mock/Infernal_SSU.tblout \
  --noali \
  --cpu "$threads" \
  database/RF01960.cm \
  results/mock/mock_sequences_full.fasta


###############################################################################
# Calculate SSU mean Phred scores
# First run the SSU-extraction section in 05_mock_error_and_index_switching.R
###############################################################################

phredsort \
  --metric avgphred \
  --header avgphred \
  --in results/mock/mock_sequences_SSU.fastq.gz \
  --out - \
| seqkit seq --name \
| sed 's/ avgphred=/\t/' \
> results/mock/mock_sequences_SSU_quality.tsv


###############################################################################
# Align the retained reads separately for each species (OTU)
# Read names are formatted as <PacBioID>|<PrimerPair>|<Polymerase>
# Inspect each alignment in AliView, remove chimeras, split it by primer and polymerase,
# and count substitutions and indels relative to the consensus as described in the manuscript
# This manual step is not automated here
###############################################################################

mkdir -p results/mock/alignments

for input_file in results/mock/alignment_input/*.fasta; do
  output_file="results/mock/alignments/$(basename "$input_file" .fasta).aligned.fasta"
  mafft --auto "$input_file" > "$output_file"
done

