#!/bin/bash

# Set directories
SOFTWARE_DIR="path_to_software"
REFERENCE_DIR="path_to_reference"
FASTQ_DIR="path_to_fastq"
GENOME_NAME="G.gynandra"

# Generate reference
$SOFTWARE_DIR/cellranger mkref --genome=$GENOME_NAME \
  --fasta=Gynandropsis_gynandra.fa \
  --genes=Gynandropsis_gynandra.gtf

# Define samples
declare -A SAMPLES
SAMPLES=(
  [Gy_nuc_sample1_force]="Gyn_nuc_1-SCI7T075-SCI5T075_HNHGLDSX5,Gyn_nuc_1-SCI7T075-SCI5T075_HTY3LDSX5"
  [Gy_nuc_sample2_force]="Gyn_nuc_2-SCI7T087-SCI5T087_HNHGLDSX5,Gyn_nuc_2-SCI7T087-SCI5T087_HNHMFDSX5,Gyn_nuc_2-SCI7T087-SCI5T087_HTY3LDSX5"
  [Gy_proto_sample1]="Gyn_proto_1-SCI7T004-SCI5T004_HNHGLDSX5,Gyn_proto_1-SCI7T004-SCI5T004_HNJ5MDSX5,Gyn_proto_1-SCI7T004-SCI5T004_HTWT7DSX5"
  [Gy_proto_sample2]="Gyn_proto_2-SCI7T016-SCI5T016_HNHGLDSX5,Gyn_proto_2-SCI7T016-SCI5T016_HNJ5MDSX5,Gyn_proto_2-SCI7T016-SCI5T016_HTWT7DSX5"
)

# Run Cell Ranger count for each sample
for ID in "${!SAMPLES[@]}"; do
  if [[ "$ID" == "Gy_nuc_sample1_force" || "$ID" == "Gy_nuc_sample2_force" ]]; then
    $SOFTWARE_DIR/cellranger count \
      --id=$ID \
      --force-cells 10000 \
      --transcriptome=$REFERENCE_DIR/$GENOME_NAME \
      --fastqs=$FASTQ_DIR \
      --sample="${SAMPLES[$ID]}"
  else
    $SOFTWARE_DIR/cellranger count \
      --id=$ID \
      --transcriptome=$REFERENCE_DIR/$GENOME_NAME \
      --fastqs=$FASTQ_DIR \
      --sample="${SAMPLES[$ID]}"
  fi

  echo "Completed processing: $ID"
done