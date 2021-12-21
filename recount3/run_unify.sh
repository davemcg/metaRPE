#!/bin/bash
module load singularity
/bin/bash /home/mcgaugheyd/git/monorail-external/singularity/run_recount_unify.sh \
		recount-unify_1.0.9.sif \
		hg38 \
		/data/mcgaugheyd/projects/nei/bharti/metaRPE/references \
		/data/mcgaugheyd/projects/nei/bharti/metaRPE/unify_output \
		/data/mcgaugheyd/projects/nei/bharti/metaRPE/pump_output \
		/data/mcgaugheyd/projects/nei/bharti/metaRPE/unify_output/recount_sample_metadata.tsv \
		6 \
		metaRPE:110
