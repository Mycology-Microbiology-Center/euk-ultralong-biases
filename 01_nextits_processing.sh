# NextITS processing used for the common-region comparisons

# Expected inputs:
#   data/reads/ultralong/*.fq.gz
#   data/reads/full_its/*.fq.gz

threads=8   # the number of CPUs to use. NB! adjust to your system

mkdir -p results/nextits results/work


###############################################################################
# Download demultiplexed reads from ENA (PRJEB123023)
###############################################################################

mkdir -p data/reads/all data/reads/ultralong data/reads/full_its

curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB123023&result=read_run&fields=submitted_ftp,submitted_md5&format=tsv" \
  | tail -n +2 \
  > data/reads/ena_files.tsv

while IFS=$'\t' read -r _ ftp md5; do
  curl -fsS -o "data/reads/all/$(basename "$ftp")" "https://${ftp}"
done < data/reads/ena_files.tsv

awk -F'\t' '{n = split($2, p, "/"); print $3 "  data/reads/all/" p[n]}' data/reads/ena_files.tsv \
  | md5sum -c --quiet

## Ultra-long reads (all primer pairs except ITS9mun) are trimmed to Euk575F-21Rngs
## all reads, including ultra-long ones, are trimmed to ITS9mun-ITS4ngsUni
for f in data/reads/all/*.fq.gz; do
  ln -sf "../all/$(basename "$f")" data/reads/full_its/
  case "$(basename "$f")" in
    ITS9mun_*) ;;
    *) ln -sf "../all/$(basename "$f")" data/reads/ultralong/ ;;
  esac
done


###############################################################################
# Step 1: ultra-long reads trimmed to Euk575F-21Rngs
###############################################################################

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  -qs "$threads" \
  --step "Step1" \
  --input "data/reads/ultralong" \
  --demultiplexed true \
  --primer_forward "TASCYGYGGTAAYWCCAGC" \
  --primer_reverse "GAGACGAGGCATTTGGCTAC" \
  --its_region "none" \
  --chimera_methods "none" \
  --hp false \
  --tj false \
  --outdir  "results/nextits/step1_ultralong" \
  -work-dir "results/work/step1_ultralong"


###############################################################################
# Step 1: full-ITS reads trimmed to ITS9mun-ITS4ngsUni
###############################################################################

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  -qs "$threads" \
  --step "Step1" \
  --input "data/reads/full_its" \
  --demultiplexed true \
  --primer_forward "GTACACACCGCCCGTCG" \
  --primer_reverse "CCTSCSCTTANTDATATGC" \
  --its_region "none" \
  --chimera_methods "none" \
  --hp false \
  --tj false \
  --outdir  "results/nextits/step1_full_its" \
  -work-dir "results/work/step1_full_its"


###############################################################################
# Step 2: ultra-long reads, VSEARCH clustering at 99% similarity
###############################################################################

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  -qs "$threads" \
  --step "Step2" \
  --chunking_n 20 \
  --chunking_id 0.7 \
  --clustering "vsearch" \
  --otu_id 0.99 \
  --merge_replicates false \
  --max_MEEP 999 \
  --max_ChimeraScore 999 \
  --lulu false \
  --data_path "results/nextits/step1_ultralong" \
  --outdir    "results/nextits/step2_ultralong" \
  -work-dir   "results/work/step2_ultralong"


###############################################################################
# Step 2: full-ITS reads, VSEARCH clustering at 99% similarity
# (chunking was not used for this smaller dataset)
###############################################################################

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  -qs "$threads" \
  --step "Step2" \
  --clustering "vsearch" \
  --otu_id 0.99 \
  --merge_replicates false \
  --max_MEEP 999 \
  --max_ChimeraScore 999 \
  --lulu false \
  --data_path "results/nextits/step1_full_its" \
  --outdir    "results/nextits/step2_full_its" \
  -work-dir   "results/work/step2_full_its"

