#!/bin/bash

# Set directories
SOFTWARE_DIR="path_to_software"
REFERENCE_DIR="path_to_reference"
FASTQ_DIR="path_to_fastq"
GENOME_NAME="F.bidentis"

# Generate reference
$SOFTWARE_DIR/cellranger mkref --genome=$GENOME_NAME \
  --fasta=Flaveria_bidentis.fa \
  --genes=Flaveria_bidentis.gtf

# Define samples
declare -A SAMPLES
SAMPLES=(
  [nuclei_apr_fb_v7]="Fb_nucleiunspun_1-SCI7T031-SCI5T031_HL3G7DSX3,Fb_nucleiunspun_2-SCI7T056-SCI5T056_HL3G7DSX3"
  [nuclei_july_fb_v7]="S1_Fb_unspun_nuc-SCI7T001-SCI5T001_HYVCWDSX3,S4_Fb_unspun_nuc-SCI7T002-SCI5T002_HYVCWDSX3"
  [Fb_proto_Jul_sample1]="S2_Fb_proto_1-SCI7T013-SCI5T013_H2CNKDSX5,S5_Fb_proto_1-SCI7T014-SCI5T014_H2CNFDSX5,S5_Fb_proto_1-SCI7T014-SCI5T014_H2CNKDSX5"
  [Fb_proto_Jul_sample2]="S3_Fb_proto_2-SCI7T025-SCI5T025_H2CNKDSX5,S6_Fb_proto_2-SCI7T026-SCI5T026_H2CNFDSX5,S6_Fb_proto_2-SCI7T026-SCI5T026_H2CNKDSX5"
  [Fb_proto_Jul_sample3]="Fb_protomanni_1-SCI7T044-SCI5T044_HL3G7DSX3,Fb_protomanni_2-SCI7T092-SCI5T092_HL3G7DSX3"
)

# Run Cell Ranger count for each sample
for ID in "${!SAMPLES[@]}"; do
  $SOFTWARE_DIR/cellranger count \
    --id=$ID \
    --transcriptome=$REFERENCE_DIR/$GENOME_NAME \
    --fastqs=$FASTQ_DIR \
    --sample="${SAMPLES[$ID]}"
  fi
  echo "Completed processing: $ID"
done