#!/bin/bash
module load singularity
mkdir pump/$1
cd pump/$1
/bin/bash /home/mcgaugheyd/git/monorail-external/singularity/run_recount_pump.sh \
	/data/mcgaugheyd/projects/nei/bharti/metaRPE/recount-rs5_1.0.6.sif \
	$1 \
	local \
	hg38 \
	6 \
	/data/mcgaugheyd/projects/nei/bharti/metaRPE/references \
	/data/mcgaugheyd/projects/nei/OGVFB_rna_seq/human/$1_R1_001.fastq.gz \
	/data/mcgaugheyd/projects/nei/OGVFB_rna_seq/human/$1_R2_001.fastq.gz \
	$2
