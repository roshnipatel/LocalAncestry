# LocalAncestry
Roshni Patel (roshnipatel@berkeley.edu)

## Description
Calling local ancestry tracts with RFMix v1.5.4 and (optionally) validating global ancestry fractions with ADMIXTURE.

### NCBI build
This pipeline was designed to call local ancestry tracts on hg38 VCFs. In practice, it can be used to call local ancestry tracts on hg37 VCFs, but you will need to obtain hg37 versions of the files in `maps/` (e.g. [here](https://github.com/armartin/ancestry_pipeline)).

### Cluster
This pipeline is designed to be run on a cluster that uses the Slurm workload manager. If this is not the case, `sm_script.sh` and `sm_slurm_config.json` will need to be reworked.

## Required files
See `scripts/snakemake_variables.py` for more details on naming conventions and directory structure.
### Programs
* [RFMix v1.5.4](https://sites.google.com/site/rfmixlocalancestryinference/)
* [ADMIXTURE](http://software.genetics.ucla.edu/admixture/) (optional)
### Data
* Phased VCF containing admixed data
* Phased VCF containing reference data
* Metadata file for admixed individuals (mapping ancestry to VCF sample ID)
* Metadata file for reference population (mapping ancestry to VCF sample ID)

## Workflow
1. Obtain required files
2. Determine whether your admixed and reference VCFs have the same chromosome naming system; otherwise you might want to uncomment rule `annotate_admix` in the Snakefile and modify accordingly to suit your purposes.
3. Edit variables in `scripts/snakemake_variables.py` as needed. (In all likelihood, only the variables under Data and Programs will need to be edited.)
4. Edit `sm_slurm_config.json` with the partitions you'll be running jobs on.
5. Edit the top of `sm_script.sh` and the top of the Snakefile to reflect the location/type of shell you're using.
6. Edit rule `all` in Snakefile depending on whether you want to run both RFMix and ADMIXTURE or only the former.
7. Make the conda environments in `envs/`
8. Unzip `maps/plink.GRCh38.genetic_map.zip`
9. Run the pipeline with `./sm_script.sh`
