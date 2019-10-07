"""
Variables and filepaths for Snakefile
"""

### Directories ###

# Directory in which to find admixed VCFs and store output files
DATA_DIR = "data/"

# Directory in which to find reference VCFs
REF_DIR = "data/1KGenomes_hg38/"

# Directory in which to find genetic maps
MAP_DIR = "data/maps/"



### Data ###

# Basename of phased VCF for admixed and reference data
ADMIX_DATA = "mesa.1331samples.genotypes.pass.phased.maf01.vcf.gz"
REF_DATA = "GRCh38.filt.1kg.phase3.v5a.biSNPs.vcf.gz"

# Files containing metadata for admixed and reference data. Must be .csv, and
# first column should be sample ID in VCF.
ADMIX_METADATA = "MESA_sample_info.csv"
REF_METADATA = "20130606_sample_info.csv"

# Column specifying race/ancestry in metadata file
REF_POP_COL = "Population"
ADMIX_POP_COL = "race1c"

# Race/ancestry values in metadata file to filter for
REF_POP_VAL = "CEU.YRI"
ADMIX_POP_VAL = "1"


### Scripts ###

# Script for filtering sample IDs based on a specific ancestry
FILTER_SCRIPT = "scripts/ancestry_filter.py"

# Script for creating RFMix input files (classes, alleles, snp_locations)
RFMIX_INPUT_SCRIPT = "scripts/make_rfmix_input.py"

# Script for collapsing RFMix output into local ancestry tracts; adapted from
# armartin on GitHub
RFMIX_OUTPUT_SCRIPT = "scripts/collapse_ancestry.py"

# Script for plotting karyograms from local ancestry tracts; adapted from Alicia
# Martin's GitHub
KARYOGRAM_SCRIPT = "scripts/plot_karyogram.py"

# Script for inferring global ancestry proportion from local ancestry tracts;
# adapted from Alicia Martin's GitHub
GLOBAL_INF_SCRIPT = "scripts/lai_global.py"

# Script for combining ancestry calls between RFMix and ADMIXTURE for each
# individual
COMB_SCRIPT = "scripts/combine_ancestry_calls.py"

# Script for combining global ancestry estimates across chromosomes
COMB_CHR_SCRIPT = "scripts/combine_chrs.py"



### Programs ###

# Path to RFMix script
RFMIX = "bin/rfmix/RFMix_v1.5.4/RunRFMix.py"

# Path to ADMIXTURE binary
ADMIXTURE = "bin/admixture"



### Variables ###

CHROMS = [str(i) for i in range(1, 23)]

# Partition chromosomes into acrocentric (no centromere) and non-acrocentric.
# For the latter, we infer local ancestry separately on each arm of the
# chromosome to decrease computational time - this is possible because there is
# no recombination at centromeres.
ACROCENTRIC = ['13', '14', '15', '21', '22']
NON_ACRO = [c for c in CHROMS if c not in ACROCENTRIC]

# Number of EM iterations to run with RFMix; can optionally be a list (will take
# much longer, though!)
EM_ITER = 0

# r2 value used for pruning ADMIXTURE input
ADMIX_R2 = "0.1"
RFMIX_R2 = "0.5"

# Number of reference populations to use
NPOP = len(REF_POP_VAL.split('.'))

# Arms of chromosome
ARMS = ['p', 'q']
