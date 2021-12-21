# NIH biowulf2

working dir: /data/mcgaugheyd/projects/nei/bharti/metaRPE

1. cd /data/mcgaugheyd/projects/nei/bharti/metaRPE
2. cp /home/mcgaugheyd/metaRPE/recount3/get_image.sh . ; bash get_image.sh # get both pump and unify singularity images
3. bash ~/git/monorail-external/get_unify_refs.sh hg38 # unify refs
4. bash ~/git/monorail-external/get_human_ref_indexes.sh
5. mkdir references; mv hg38 references/; mv hg38_unify references/ # mv pump and unify ref folders to subfolder
6. cp /home/mcgaugheyd/metaRPE/recount3/run_pump.sh . ;cp /home/mcgaugheyd/metaRPE/recount3/pump_commands.sh . 
7. bash pump_commands.sh # runs all pump jobs (invoking the run_pump.sh script)
  - outputs in pump/[lane_file] so you can simultaneously run the pump commands and avoid snakemake rage
  - wait until the job finishes (the slower ones take 3-4 hours)
8. mkdir pump_output; rsync --progress -rav pump/*/output/ pump_output # consolidate the pump outputs to one directory
9. mkdir unify_output; cd unify_output; cp /home/mcgaugheyd/git/metaRPE/data/recount_sample_metadata.tsv . ; cp /home/mcgaugheyd/metaRPE/recount3/run_unify.sh .
10. sbatch --cpus-per-task 6 --mem=32G run_unify.sh
