# WTF are you doing?

This is a RNA-seq "meta" analysis where I'm taking many internal RNA-seq datasets and glomming
together so we can ask some new questions. As we likely wish to
ask some questions that involve bringing in datasets from external sources (e.g. some of [GTEx](http://gtexportal.org)) and
likely other bits and pieces across the SRA universe it seems logical to use the [recount3](http://rna.recount.bio) which has,
briefly, reprocessed all of SRA (human/mouse as of early 2021??). The key bit is that their processing pipeline
([monorail](https://github.com/langmead-lab/monorail-external)) has been munged into a container so you
can, in relatively few steps, run their exact quant pipeline (which is [STAR](https://github.com/alexdobin/STAR) /
[salmon](https://combine-lab.github.io/salmon/) based). Thus you can compare your custom RNA-seq data against
public RNA-seq and at least the compute part is identical.

# How to run monorail pump and unify
## Done on NIH biowulf2

## Observations, some unfounded and poorly understand
1. Do not have fastq files (local) that have a non-character / non-number within two of the end of the "core name"
        - e.g. A1_1_R1.fastq.gz and A1_2_R2.fastq.gz seem to cause issues with `unify` as `pump` uses the last two of the name (in this case _1) as a subfolder. Which *seems* to cause issues with some concatenation steps later
2. Just make the project names one case, maybe lower? I had some weird issues when I have project names like "Bob" which were perhaps resolved when I changed them to "bob"?
3. **VERY IMPORTANT:** `unify` will use *ALL THE FOLDERS/SAMPLES* in the `pump` output as input for `unify`. So you cannot "mix and match" custom unify output by just messing with the sample metadata file. 

## Workflow
1. `cd /data/mcgaugheyd/projects/nei/bharti/metaRPE`
2. `cp /home/mcgaugheyd/metaRPE/recount3/get_image.sh . ; bash get_image.sh # get both pump and unify singularity images`
3. `bash ~/git/monorail-external/get_unify_refs.sh hg38 # unify refs`
4. `bash ~/git/monorail-external/get_human_ref_indexes.sh`
5. `mkdir references; mv hg38 references/; mv hg38_unify references/ # mv pump and unify ref folders to subfolder`
6. `cp /home/mcgaugheyd/git/metaRPE/recount3/run_pump.sh . ;cp /home/mcgaugheyd/git/metaRPE/recount3/pump_commands.sh . `
7. `bash pump_commands.sh # runs all pump jobs (invoking the run_pump.sh script)`
    - outputs in pump/[lane_file] so you can simultaneously run the pump commands and avoid overlapping snakemake job rage
    - wait until the job finishes (the slower ones take 2-3 hours)
8. `mkdir pump_output; rsync --progress -rav pump/*/output/ pump_output # consolidate the pump outputs to one directory`
9. `mkdir unify_output; mv recount-unify_1.0.9.sif unify_output; cd unify_output; cp /home/mcgaugheyd/git/metaRPE/data/recount_sample_metadata.tsv . ; cp /home/mcgaugheyd/metaRPE/recount3/run_unify.sh . `
10. `sbatch --cpus-per-task 6 --mem=32G run_unify.sh`

# How to munge monorail-unify output into recount3 format

Why bother? Well, the unify output is **not** straight counts, but rather the sum of the base pair level coverage. [recount3](http://rna.recount.bio)
has tooling to do the conversion. I briefly looked into whether I could "directly" do the transform but after a
few minutes of traversing through their code, it seemed easier to just make the *RangedSummerizedExperiment* (RSE)
data structure by moving the unify outputs around. Plus the RSE is
a more portable and useful format should someone else want to use the data.

OK, so *super* briefly you have to move the unify outputs around a bit (to my local computer in my example) so recount3 can *import* the data and make the RSE
and then you can run (in R)`recount3::transform_counts` to get the straight counts for downstream use.

Their example directory is here: http://snaptron.cs.jhu.edu/data/temp/recount3test/

Fairly complete instructions: https://github.com/langmead-lab/monorail-external#loading-custom-unifier-runs-into-recount3

*My* directory structure (gaps are where I've deleted lines to make a touch shorter to view):
```
└── human
    ├── annotations
    │   ├── exon_sums
    │   │   ├── human.exon_sums.G026.gtf.gz
    │   │   └── human.exon_sums.G029.gtf.gz
    │   └── gene_sums
    │       ├── human.gene_sums.G026.gtf.gz
    │       └── human.gene_sums.G029.gtf.gz
    ├── data_sources
    │   └── metaRPE
    │       ├── base_sums
    │       ├── exon_sums
    │       │   ├── de
    │       │   │   └── davide
    │       │   │       ├── metaRPE.exon_sums.davide.ERCC.gz
    │       │   │       ├── metaRPE.exon_sums.davide.F006.gz
    │       │   │       ├── metaRPE.exon_sums.davide.G026.gz
    │       │   │       ├── metaRPE.exon_sums.davide.G029.gz
    │       │   │       ├── metaRPE.exon_sums.davide.R109.gz
    │       │   │       └── metaRPE.exon_sums.davide.SIRV.gz
    │       │   ├── hi
    │       │   │   └── ruchi
    │       │   │       ├── metaRPE.exon_sums.ruchi.ERCC.gz
    │       │   │       ├── metaRPE.exon_sums.ruchi.F006.gz


    │       ├── gene_sums
    │       │   ├── de
    │       │   │   └── davide
    │       │   │       ├── metaRPE.gene_sums.davide.ERCC.gz
    │       │   │       ├── metaRPE.gene_sums.davide.F006.gz
    │       │   │       ├── metaRPE.gene_sums.davide.G026.gz
    │       │   │       ├── metaRPE.gene_sums.davide.G029.gz


    │       ├── junctions
    │       │   ├── de
    │       │   │   └── davide
    │       │   │       ├── metaRPE.junctions.davide.ALL.ID.gz
    │       │   │       ├── metaRPE.junctions.davide.ALL.MM.gz
    │       │   │       ├── metaRPE.junctions.davide.ALL.RR.gz

    │       └── metadata
    │           ├── de
    │           │   └── davide
    │           │       ├── metaRPE.metaRPE.davide.MD.gz
    │           │       ├── metaRPE.recount_project.davide.MD.gz
    │           │       └── metaRPE.recount_qc.davide.MD.gz

    │           ├── in
    │           │   └── qin
    │           │       ├── metaRPE.metaRPE.qin.MD.gz
    │           │       ├── metaRPE.recount_project.qin.MD.gz
    │           │       └── metaRPE.recount_qc.qin.MD.gz
    │           ├── la
    │           │   └── karla
    │           │       ├── metaRPE.metaRPE.karla.MD.gz
    │           │       ├── metaRPE.recount_project.karla.MD.gz
    │           │       └── metaRPE.recount_qc.karla.MD.gz
    │           └── metaRPE.recount_project.MD.gz
    └── homes_index
```

## Notes:
1. Top level folder is organism (human or mouse)
2. Fill `annotation` folder by `wget`ing the files from Langmead and co
    - Example for G026 gene sums: `wget http://duffel.rail.bio/recount3/human/new_annotations/gene_sums/human.gene_sums.G026.gtf.gz`
3. `metaRPE` is my project name. monorail/recount tend to use `sra` as theirs.
4. The `data_sources` folder is filled by `rsync`ing various unify output folders
    - in my case it is `/data/mcgaugheyd/projects/nei/bharti/metaRPE/unify_output`
    - `base_sums` is currently empty as recount3 doesn't need the bigWig files (which are in the **pump** outputs)
    - `exon_sums_per_study` renamed to `exon_sums`
    - `gene_sums_per_study` renamed to `gene_sums`
    - `junction_counts_per_study` renamed to `junctions`
    - `metadata` doesn't need to be renamed
        - **HOWEVER** you need to "hand make" `metaRPE.recount_project.MD.gz`
          - My janky ass code that I ran in the unify output folder (see path above):
          - `zcat metadata/*/*/*recount_project.* | head -n 1 | gzip > metadata/metaRPE.recount_project.MD.gz # this grabs just the header`
          - `zcat metadata/*/*/*recount_project.* | grep -v rail_id  | gzip >> metadata/metaRPE.recount_project.MD.gz # copy the rest of the meta sans the headers`
5. The `homes_index` file is just a text file with `data_sources/metaRPE` in it
    - replace `metaRPE` with whatever your "project" name is (again, `sra` is commonly used by the monorail/recount team)


# tldr add new files
1. rsync the fastq to `/data/mcgaugheyd/projects/nei/OGVFB_rna_seq/[organism]`

        - `rsync -rav  --progress /data/mcgaugheyd/projects/nei/unicorns/fastq/* /data/mcgaugheyd/projects/nei/OGVFB_rna_seq/human/`
2. `cd /data/mcgaugheyd/projects/nei/bharti/metaRPE` (step 1 of workflow)
3. Run `pump` a la: (step 7 of workflow)

        - `sbatch --mem=35G --cpus-per-task 6 --time=8:00:00 run_pump.sh fastq_name_except_ending project_name`
        - see `pump_commands.sh` 
        - or investigate my prototype snakefile at ~/git/monorail-external/Snakerail
4. rsync --progress -rav pump/*/output/ pump_output # consolidate the pump outputs to one directory`
6. `mv unify_output unify_output_OLD; mkdir unify_output; mv unify_output_OLD/recount-unify_1.0.9.sif unify_output; cd unify_output; cp /home/mcgaugheyd/git/metaRPE/data/recount_sample_metadata.tsv . ; cp /home/mcgaugheyd/git/metaRPE/recount3/run_unify.sh . `
7. `sbatch --cpus-per-task 6 --mem=32G run_unify.sh`
8. Set up RSE for recount3:

        - `cp ~/git/metaRPE/src/make_rse.sh .`
        - `bash make_rse.sh bharti`
9. Add metadata a la in `analysis/scavenging.Rmd`
10. Tweak and run `src/monorail_to_rse.R` to add the metadata to the RSE
