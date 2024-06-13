
################################################################################
### Set up
################################################################################

library(argparser)
library(rmarkdown)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("Generate trial specific project level SNiPDx reports")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (choose from 'RP3500-01', 'RP3500-03' ,'Non-Clinical', 'RP6306-01', 'RP6306-02')")
p <- add_argument(p, arg = "--processingDir",
                  help = "specify path to working directory to export results and data folder ie) SNiPDx/[trial]/processing/")
p <- add_argument(p, arg = "--inInventory",
                  help = "specify path to latest genosity inventory file")
p <- add_argument(p, arg = "--inSamplesFolder",
                  help = "specify path to local trial specific Genosity facets output folder")
p <- add_argument(p, arg = "--inVendorData",
                  help = "specify path to store local trial specific Genosity stats and variants folder ie.) in [SEQRUN]/RESULTS/")
p <- add_argument(p, arg = "--in_nf_snipdx",
                  help = "specify path nf-snipdx repo folder (ensure current branch is set to master)")
# p <- add_argument(p, arg = "--panalAnnotations",
#                   help = "specify path to gene panel annotaion file ie.) nf-snipdx/refs/panel_genes_annotated-v2.tsv")
# p <- add_argument(p, arg = "--inExtensions",
#                   help = "specify path to file containing snipdx report specific file")

argv <- parse_args(p)

trial <- argv$trial
inInventory <- argv$inInventory
inSamplesFolder <- argv$inSamplesFolder
inVendorData <- argv$inVendorData
processingDir <- argv$processingDir
# panalAnnotations <- argv$panalAnnotations
# inExtensions <- argv$inExtensions
in_nf_snipdx <- as.character(argv$in_nf_snipdx)

################################################################################
### Set testing arguments 
################################################################################
if(FALSE){
  trial <- "Non-Clinical"
  
  inInventory <- ("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/INVENTORY/rp_nonclin_genosity_inv_annotated_2022-12-06.csv")
  
  # inSamplesFolder <- ("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/ARCHER/")
  # 
  # inVendorData <- ("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/ARCHER/")
  
  inSamplesFolder <- ("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/")
  
  inVendorData <- ("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/")
  
  panalAnnotations <- ("/ClinBio/SP-ClinicalBioinformatics/ntran/analysis/SNIPDx/nf-snipdx-dev/nf-snipdx/refs/panel_genes_annotated-v2.tsv")
  
  processingDir <- "/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/processing/"
  
  inExtensions <- "/ClinBio/SP-ClinicalBioinformatics/ntran/analysis/SNIPDx/nf-snipdx-dev/nf-snipdx/src/utils/snipdx_reports.R"
  
  in_nf_snipdx <- "/ClinBio/SP-ClinicalBioinformatics/ntran/analysis/SNIPDx/nf-snipdx-dev/nf-snipdx/"
  
}

################################################################################
### Match trials with gcp naming conventions
################################################################################
trial_to_gcp <- list("RP3500-01" = "350001",
                     "RP3500-03" = "350003",
                     "RP6306-01" = "rp630601",
                     "RP6306-02" = "rp630602",
                     "Non-Clinical" = "nonclin")

################################################################################
### Sync GCP watch folder - regex needs work
################################################################################
### Sync folder structures from GCP AWS watchfolder to rserver, read files to build project level SNiPDx reports

# example:
# gsutil rsync -mr -x ".*\.vcf.gz$|.*\.fastq.gz$|.*\.bam$|.*\.bai$" gs://quartz-bio/prd/repare_nonclin_bmdb/data/source/GENOSITY/ /ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/

# -m: activates multithreading processes for faster syncing
# -r: syncs folders recursively 
# -x: specifies regex expression to exclude files matching pattern
rsync <- paste0("gsutil -m rsync -r -x")
# exclude vcfs, fastqs, bams, and bai and duplicate project level reports folder (REPORTS/)
regex_exclude <- paste("'REPORTS/.*.$|.*./FACETS/GENOSITY/.*.$|WGS/.*.$|.*./METADATA/.*.$|.*\\.vcf.gz$|.*\\.fastq.gz$|.*\\.bam$|.*\\.bai$'")

gcp <- paste0("gs://quartz-bio/prd/repare_",trial_to_gcp[[trial]], "_bmdb/data/source/GENOSITY/")
# local <- paste0("/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/", trial , "/vendor-data/gcp/")
command <- paste(rsync,regex_exclude, gcp, inVendorData)
system(command)

################################################################################
### Sync GCP watch folder - regex needs work
################################################################################

### load in inventory file to get list of all possible patients 
Inventory <- readr::read_csv(inInventory)
#new_batch <- tail(str_sort(unique(Inventory$SEQRUN)),n = 1)

# knit_root_dir <- getwd()

# set wd to nf-snipdx bitbucket folder
setwd(in_nf_snipdx)
print(getwd())

if(trial == "Non-Clinical"){
  projects <- unique(Inventory$`Project`)
  
  # individual project test
  #projects <- "MDACC_LMS_RNASEH2"
  
  results_dir <- paste0(Sys.Date(), "_",trial, "_", "SNiPDx_by_Projects")
  output_dir <- paste0(processingDir, 
                       "/results/",results_dir)
  
  
  for(currProject in projects){
    rmarkdown::render("./src/generate_SNiPDx_project_level_reports_clinical_nonclinical_sgz_segments.Rmd", 
                      knit_root_dir = in_nf_snipdx,
                      output_dir = output_dir,
                      output_file = paste0(output_dir, "/",
                                           Sys.Date(),"_", 
                                           trial,"_", 
                                           'SNiPDx',"_", 
                                           currProject, 
                                           '.html'),
                      params=list(
                        inExtensions = "./src/utils/snipdx_reports.R",
                        inInventory = inInventory,
                        inSamplesFolder = inSamplesFolder,
                        inVendorData = inVendorData,
                        currTrial = trial,
                        currPatient = currProject,
                        panalAnnotations = "./refs/panel_genes_annotated-v2.tsv",
                        processingDir = processingDir)
    )
  }
  ### clinical trials
}else{
  patients <- unique(Inventory$SUBJID)
  
  #single patient 
  #patients <- c("1006-0332")
  
  results_dir <- paste0(Sys.Date(), "_",trial, "_", "SNiPDx_by_Patients")
  output_dir <- paste0(processingDir, 
                       "/results/",results_dir)
  
  for(currPatient in patients){
    rmarkdown::render("./src/generate_SNiPDx_project_level_reports_clinical_nonclinical_sgz_segments.Rmd", 
                      knit_root_dir = in_nf_snipdx,
                      output_dir = output_dir,
                      output_file = paste0(output_dir, "/",
                                           Sys.Date(),"_", 
                                           trial,"_", 
                                           'SNiPDx',"_", 
                                           currPatient, 
                                           '.html'),
                      params=list(
                        inExtensions = "./src/utils/snipdx_reports.R",
                        inInventory = inInventory,
                        inSamplesFolder = inSamplesFolder,
                        inVendorData = inVendorData,
                        currTrial = trial,
                        currPatient = currPatient,
                        panalAnnotations = "./refs/panel_genes_annotated-v2.tsv",
                        processingDir = processingDir)
    )
  }
  
}

################################################################################
### Sync back to GCP watch bucket
################################################################################
gcp_destination <- paste0(gcp, "REPORTS/", results_dir)
system(paste("gsutil -m rsync -r", output_dir, gcp_destination ))
