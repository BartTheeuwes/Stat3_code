#!/bin/bash
# -------------------------------------------------------
#title           :run_cellranger_count.sh
#author          :hpb29
#date            :20180404
#version         :0.5    
#description     :Controls the submition of 10x libraries
#                 via SLURM to 'cellranger count'
#usage           :Edit variables below accordingly and run
# --------------------------------------------------------




# ================== PROBABLY THERE WILL BE NO NEED TO CHANGE THESE ==================

declare -A library
refdir=/rds/project/rds-SDzz0CATGms/references/10x/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=mm10_tomato


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/SLX-21143_fastq/
# libraries and respective expected sizes
# ---------------------------------------
#library[SLX-20795_SITTG10_HKTG2DRXY]=4465
library[SLX-21143_SITTA2_HTJH3DSX2]=13784
library[SLX-21143_SITTA4_HTJH3DSX2]=946
library[SLX-21143_SITTB2_HTJH3DSX2]=8573
library[SLX-21143_SITTB4_HTJH3DSX2]=4561
library[SLX-21143_SITTC4_HTJH3DSX2]=4745
library[SLX-21143_SITTD4_HTJH3DSX2]=1513
library[SLX-21143_SITTD5_HTJH3DSX2]=3463
library[SLX-21143_SITTE3_HTJH3DSX2]=3463
library[SLX-21143_SITTE5_HTJH3DSX2]=3463
library[SLX-21143_SITTF3_HTJH3DSX2]=3463
library[SLX-21143_SITTF5_HTJH3DSX2]=3463
library[SLX-21143_SITTG1_HTJH3DSX2]=3463
library[SLX-21143_SITTG3_HTJH3DSX2]=3463
library[SLX-21143_SITTG5_HTJH3DSX2]=3463
library[SLX-21143_SITTH1_HTJH3DSX2]=3463
library[SLX-21143_SITTH3_HTJH3DSX2]=3463

# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/SLX-21143_fastq/cellranger_count.sh
# ========================================================================




echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
