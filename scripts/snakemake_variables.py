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
REF_POP_VAL = "CEU,YRI"
ADMIX_POP_VAL = "3"


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

# Number of reference populations to use
NPOP = len(REF_POP_VAL.split(','))

# Arms of chromosome
ARMS = ['p', 'q']

# Sets of haplotype
HAPL = ['A', 'B']

# Sample IDs of all individuals for which local ancestry analysis should be
# performed (generally all admixed individuals)
INDIV = ['NWD163203', 'NWD838454', 'NWD447217', 'NWD506902', 'NWD516805', \
        'NWD598769', 'NWD276731', 'NWD857305', 'NWD573145', 'NWD378807', \
        'NWD773408', 'NWD862131', 'NWD696182', 'NWD601219', 'NWD445484', \
        'NWD245122', 'NWD574025', 'NWD962792', 'NWD140824', 'NWD347140', \
        'NWD388542', 'NWD409874', 'NWD641240', 'NWD620567', 'NWD754070', \
        'NWD338378', 'NWD660071', 'NWD452926', 'NWD233967', 'NWD878543', \
        'NWD484275', 'NWD913991', 'NWD643775', 'NWD453328', 'NWD865962', \
        'NWD524336', 'NWD728777', 'NWD652755', 'NWD346975', 'NWD795057', \
        'NWD630602', 'NWD467083', 'NWD610801', 'NWD158237', 'NWD600992', \
        'NWD381518', 'NWD710489', 'NWD220560', 'NWD263184', 'NWD343349', \
        'NWD998816', 'NWD897019', 'NWD612175', 'NWD606337', 'NWD404517', \
        'NWD486871', 'NWD170091', 'NWD902298', 'NWD402569', 'NWD251680', \
        'NWD256042', 'NWD242199', 'NWD141713', 'NWD816583', 'NWD487726', \
        'NWD462712', 'NWD916599', 'NWD856383', 'NWD576327', 'NWD815786', \
        'NWD571531', 'NWD890143', 'NWD901710', 'NWD171365', 'NWD782813', \
        'NWD993431', 'NWD847600', 'NWD716639', 'NWD293024', 'NWD574916', \
        'NWD862869', 'NWD557722', 'NWD299581', 'NWD135441', 'NWD345273', \
        'NWD546564', 'NWD774561', 'NWD969537', 'NWD449948', 'NWD716470', \
        'NWD945886', 'NWD417668', 'NWD362144', 'NWD671182', 'NWD281909', \
        'NWD210347', 'NWD998609', 'NWD353250', 'NWD781750', 'NWD426008', \
        'NWD587844', 'NWD133198', 'NWD149502', 'NWD975921', 'NWD181781', \
        'NWD661828', 'NWD712674', 'NWD478709', 'NWD337903', 'NWD288321', \
        'NWD158184', 'NWD559886', 'NWD737180', 'NWD419692', 'NWD479130', \
        'NWD797979', 'NWD352494', 'NWD752087', 'NWD390054', 'NWD479314', \
        'NWD626329', 'NWD521616', 'NWD878668', 'NWD268401', 'NWD138795', \
        'NWD159679', 'NWD609750', 'NWD392992', 'NWD730299', 'NWD993659', \
        'NWD108116', 'NWD787458', 'NWD657624', 'NWD105648', 'NWD311932', \
        'NWD640412', 'NWD866219', 'NWD687607', 'NWD788130', 'NWD519039', \
        'NWD678520', 'NWD212434', 'NWD687323', 'NWD523161', 'NWD730873', \
        'NWD828733', 'NWD168541', 'NWD164446', 'NWD930303', 'NWD457909', \
        'NWD363273', 'NWD547021', 'NWD314678', 'NWD950794', 'NWD677387', \
        'NWD542871', 'NWD879683', 'NWD795363', 'NWD346624', 'NWD262638', \
        'NWD419335', 'NWD177635', 'NWD462567', 'NWD553021', 'NWD563349', \
        'NWD822723', 'NWD944827', 'NWD959607', 'NWD681170', 'NWD490768', \
        'NWD685703', 'NWD907433', 'NWD712032', 'NWD677776', 'NWD599306', \
        'NWD684484', 'NWD826684', 'NWD428381', 'NWD174966', 'NWD938526', \
        'NWD328083', 'NWD318497', 'NWD731252', 'NWD412116', 'NWD640560', \
        'NWD639242', 'NWD599371', 'NWD720391', 'NWD431959', 'NWD584394', \
        'NWD344310', 'NWD517673', 'NWD527323', 'NWD418773', 'NWD857221', \
        'NWD836765', 'NWD396724', 'NWD125656', 'NWD912070', 'NWD597815', 'NWD512818']
