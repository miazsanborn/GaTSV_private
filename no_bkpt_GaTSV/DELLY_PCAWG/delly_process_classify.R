## stop test from exploding...

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)


.libPaths("/data/beroukhim1/oumayma/Rlibs")

require(BiocGenerics)
require(caTools)
require(data.table)
require(e1071)
require(GenomeInfoDb)
require(GenomicRanges)
require(gUtils)
require(IRanges)
require(parallel)
require(rlang)
require(ROCR)
library(rstudioapi)
require(S4Vectors)
require(stats4)
require(stringr)

print("Setting Working Directory")
setwd("/data/beroukhim1/oumayma/GaTSV")

source('/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/annotate_pcawg.R')

#ANNOTATE 
cat('Loading reference files...\n')
gnomad_hg38 = readRDS('data/gnomAD.v4.hg38.rds')
gnomad_hg19 = readRDS('data/gnomAD.v4.hg19.liftover.rds')
LINE_dt_hg38 = readRDS('data/repeatmasker.hg38.LINE.bed')
SINE_dt_hg38 = readRDS('data/repeatmasker.hg38.SINE.bed')
LINE_dt_hg19 = readRDS('data/repeatmasker.hg19.LINE.bed')
SINE_dt_hg19 = readRDS('data/repeatmasker.hg19.SINE.bed')
hg19_genes = readRDS('data/gencode.genes.hg19.rds')
hg19_exons=readRDS('data/gencode.exons.hg19.rds')
hg38_genes=readRDS('data/gencode.genes.hg38.rds')
hg38_exons=readRDS('data/gencode.exons.hg38.rds')
reptimedata_hg19 = readRDS('data/reptime.hg19.rds')
reptimedata_hg38 = readRDS('data/reptime.hg38.rds')
scaling_mat <- fread("data/scalingmatrix.txt")

#############################################
tmp_env <- new.env()

# Load the file into the temporary environment
load("/data/beroukhim1/oumayma/no_bkpt_GaTSV/GaTS_BpA_031826.rda", envir = tmp_env)

# Extract the first (and likely only) object inside it and name it GaTSV
GaTSV <- tmp_env[[ls(tmp_env)[1]]]

# Optional: Clean up the temporary environment
rm(tmp_env)
#############################################

#running the classifier on example data
metadata <- fread("/data/beroukhim1/oumayma/logs/pcawg_tp53_metadata.txt") #metadata file that contains the sample_ids (same as basename of filepath without '.vcf' extension) and associated tp53_mutation_status

# Read patient prefix passed by Slurm array
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No patient prefix provided by Slurm.")

patient_prefix <- args[1]
sample_id <- basename(patient_prefix) 

# Reconstruct the full paths for both VCFs
germline_vcf <- paste0(patient_prefix, ".germline.sv.vcf")
somatic_vcf <- paste0(patient_prefix, ".somatic.sv.vcf")

run_GaTSV(
  germline_pth = germline_vcf, 
  somatic_pth = somatic_vcf, 
  sample = sample_id, 
  n_cores = 1, 
  genome = 'hg19', 
  output_path = '/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/'
)



