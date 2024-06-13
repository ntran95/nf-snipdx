# nf-snipdx Analysis
User Manual For Repare SNiPDX analysis pipeline. This is a version refactored in NextFlow from the
original version available at https://bitbucket.org/reparecompbio/osloh/src/master/

Version: 2.5.0

Copyright 2022, Repare Therapeutics
______________________

## Table of Contents:

- [Updates](CHANGELOG.md)
- [Overview](#overview)
- [Installation](#installation)
- [Command line parameters](#command-line-parameters)
- [Running the pipeline](#running-the-pipeline)
- [wiki](docs/wiki.md)
- [Panel of Normal Versions](#panel-of-normal-versions)

-------------------------------------------------------------------------------------------------
## Overview

The SNiPDX pipeline requires performing the following steps. The nf-SNiPDX nextflow pipeline performs steps 3-7, 
and the `refs` folder contains the files for the 'identified germline SNPs' (from step 1) and 'PON' (from step 2),
references for PureCN are located in GCP  
(`gs://repare-reference-data/purecn/for_snipdx/`), if running the pipeline
locally these files will need to be downloaded locally

1) Identify germline SNPs (from Gnomad database, not covered by this repository - see [the Wiki page](](docs/wiki.md)))

2) Analyze set of normal (control) samples. Perform pileup at SNP positions, using .bam files. Select SNPs that are 
   appropriately covered (possibly stratify by sample QC). Only keep clean SNPs. Eliminate SNPs that are not 
   heterozygous in normal samples, as this may be indicative of mapping problems. (not covered by this repository - 
   see [the Wiki page](docs/wiki.md))

3) Perform a pileup (at SNPs from 1) and gather insert size metrics. Simulate a tumor/normal input file from Facets, 
   where normal coverage is derived from the coverage across normal samples.

4) Perform Facets analysis. For the tumor-only analysis, use Facet's tumor-only settings.

5) After facets, save a table specifying copy number (major, minor) in genomic segments, specified by chromosome, 
   start and end positions. Make standard Facets plot and plots with different layout so they can be shown 
   side-to-side with Ascat plots.

6) Run PureCN and prerequisite freebayes variant calling annotated with dbSNP to get germline variants required for
PureCN modeling

7) Get diplogr from PureCN (mean log2 ratio of segments at copy number 2 in PureCN model) and rerun facets with
   dipLogR parameter set to this value
   
8) Compile sample level facets results and build project level reports, partitioned by patient (for clinical samples) or project id (for non-clinical samples) for each trial

__________________________________________________________________________________________________
## Installation

This pipeline requires [Docker](https://www.docker.com/) and [Nextflow](https://www.nextflow.io/):

You can clone the repository locally using the following command:

```
git clone https://[username]@bitbucket.org/reparecompbio/nf-snipdx.git
```

This pipeline uses a custom container. To build the custom containers within docker
navigate to the `docker/snipdx_tools/` directory and build using the following commands:

```
docker build -t snipdx_tools:prod_2.2 .
```

To build a version of snipdx_tools that allows for creating chromosome plots
navigate to the `docker/snipdx_tools/` directory and build using the following commands:

```
docker build -t snipdx_tools:prod_2.3 .
```

To build the reporting script docker navigate to the `docker/RMarkDown/` directory and run the following:

```
docker build -t rmarkdown:prod_1.0 .
```

To build docker containers required for PureCN steps of the pipeline:  
   
Navigate to the `docker/python_tools` directory and run the following:

```
docker build -t python_tools:1.0 .
```
Navigate to the `docker/snpeff` directory and run the following:

```
docker build -t snpeff:1.0 .
```

##### Notes on Requirements:

The Docker container includes [snp-pileup](https://bioconda.github.io/recipes/snp-pileup/README.html) and 
[r-facets](https://bioconda.github.io/recipes/r-facets/README.html), both of which have public repositories. 
These are brought included in a single container for this pipeline because the plotting R scripts also require package 
'RcppRoll'. Instead of creating a separate custom container for each dependency we pull all tools into a single
container. To run the pipeline locally without Docker all the aforementioned dependencies will be required.

The previous repo is available at Bitbucket:
https://bitbucket.org/reparecompbio/osloh/src/master/

-------------------------------------------------------------------------------------------------
## Command line parameters

When using the provided config there are several profiles set up for basic use cases

|Profile Name|Description                                                                                                                   |
|-----------|------------------------------------------------------------------------------------------------------------------------------|
|gcp        |outlines google life science connection parameters for the repare-database project                                            |
|gcpClin    |outlines google life science connection parameters for the repare-clinical project                                            |
|ref        |uses reference files in the local filesystem (files included in repo)                                                         |
|local      |use this profile for local compute. This will lower the system requirements and set gcp region parameters for local execution.|
|gcpRef     |uses reference files on GCP for PureCN

Profiles are used as a shortcut for specifying multiple common command line parameters simultaneously, 
and can be specified in a comma separated list with the following syntax:

```
nextflow -c snipdx.config run snipdx.nf -profile gcp,ref 
```

Individual command line parameters:

| Parameter                | Description                                                                                                                                             | Default                   |
|--------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|
| output_folder            | path to output the results to (defaults to present working directory)                                                                                   | current working directory | 
| vcf_fp                   | path to ref vcf file with targetted SNPs                                                                                                                |                           |
| panel_ga                 | path to panel genes annotation file (tsv)                                                                                                               |                           |
| pon_version              | string of pon version to use. One of 'pon.012222' or 'pon.092722'                                                                                       |                           |
| pon_fp_012222            | path to pon.012222 (sequenced pre SETD2 artifact, see below for sequencing runs affected)                                                               |                           |
| pon_fp_092722            | path to pon.092722 (sequenced post SETD2 artifact)                                                                                                      |                           |
| metadata_file            | path to the metadata file. Use this to specify multiple bams, including parameters sample IDs and trial in inventory if generating markdown reports     |                           |
| experiment_name          | string of the experiment name for the batch, used to name the folder containing the cohort of sample results, required when --output_to_quartz is false |                           |
| clust_info               | path to cluster info file (RData file)                                                                                                                  |                           |
| local_insertsize_norm    | boolean, whether to normalize around local SNP insert sizes                                                                                             | false                     |
| pon_target_location      | path to the pon target location file for use in local insert size normalization                                                                         |                           |
| rmarkdown_report         | boolean, whether to generate R Markdown report                                                                                                          | false                     |
| report_only              | boolean, whether to only generate the R Markdown report. Requires metadata file and experiment name                                                     | false                     |
| purecn_normal_db         | path to PureCN PON file                                                                                                                                 |                           |
| purecn_interval_file     | path to txt file with PureCN intervals                                                                                                                  |                           |
| purecn_interval_bed      | path to bed file with PureCN intervals                                                                                                                  |                           |
| purecn_mapping_bias_file | path to file with PON mapping bias for PureCN                                                                                                           |                           |
| dbsnp                    | path to vcf file with dbsnp annotations                                                                                                                 |                           |
| dbsnp_index              | path to dbsnp vcf .tbi index file                                                                                                                       |                           |
| ref_genome               | path to reference genome fasta                                                                                                                          |                           |
| ref_genome_index         | path to index file for reference genome fasta                                                                                                           |                           |
| purecn_min_af            | float, minimum variant allele frequency for PureCN to consider a variant in copy number modeling                                                        | 0.05                      |
| purecn_min_count         | integer, minimum number of reads in tumor + normal to keep an interval                                                                                  | 100                       |
| purecn_seg_algorithm     | string, segmentation algorithm for PureCN. One of 'CBS', 'PSCBS', 'Hclust'                                                                              | PSCBS                     |
| purecn_alpha             | float, significance of break points                                                                                                                     | 0.005                     |
| purecn_min_purity        | float, minimum considered purity                                                                                                                        | 0.10                      |
| purecn_max_purity        | float, maximum considered purity                                                                                                                        | 0.95                      |
| purecn_min_ploidy        | float, minimum considered ploidy                                                                                                                        | 1.4                       |
| purecn_max_ploidy        | float, maximum considered ploidy                                                                                                                        | 6.0                       |
| purecn_max_ascn          | int, maximum allele specific copy                                                                                                                       | 6                         |
| purecn_model             | string, model used to fit variants. One of 'beta' or 'betabin'                                                                                          | betabin                   |
| purecn_max_nonclonal     | float, maximum genomic fraction assigned to a subclonal copy number state                                                                               | 0.30                      |
| output_to_quartz         | boolean, whether to write outputs to quartz-bio GCP bucket under individual sample directories                                                          | false                     |


The following Trials are acceptable and recognized in nf-SNiPDX:

| Trial Name  | Description                                                  |
|-------------|--------------------------------------------------------------|
| nonclinical | trial name and vendor data path based on nonclinical samples |
| rp350001    | trial name and vendor data path based on rp3500-01 trial     |
| rp350003    | trial name and vendor data path based on rp3500-03 trial     |
| rp630601    | trial name and vendor data path based on rp6306-01 trial     |
| rp630602    | trial name and vendor data path based on rp6306-02 trial     |


------------------------------------------------------------------------------------------------
## Running the pipeline

#### Using GCP profiles


First make sure you have [gcloud CLI](https://cloud.google.com/sdk/gcloud) installed, and have `gcloud` configured to 
the appropriate project based on where the pipeline config points. If you currently use `gcloud`/`gsutils` for 
interacting with GCP, then these first commands may not need to be run the pipeline, but the next likely still will:

```
gcloud auth login
```

Following the above login you should set the appropriate google project for working in.

```
gcloud config set project repare-database
```

You need the environment you run the pipeline in to have access to the GCP backend. The most common way to do this
will be to pick up your end-user Google credentials from your workstation. You can create these by running the command:

```
gcloud auth application-default login
```

then running through the authentication flow. This will write a credential file to your gcloud configuration directory 
that will be used for any tool you run on your workstation that picks up default credentials. From here you should
be able to use GCP related profiles for reference files and compute.


#### Process Samples

Previous versions of the pipeline allowed for running samples individually without using a metadata file, in the 
current version of the pipeline individual samples should be run using a metadata file with one entry for the sample

Example command line for running samples:

```
nextflow -c snipdx.config run snipdx.nf -profile ref --metadata_file "path/to/metadata.tsv" --experiment_name "test_experiment" --pon_version pon.092722
```

Pipeline results can also be written to the GCP quartz-bio bucket using the  
`--output_to_quartz` flag, these will be written to a location under the individual samples genosity results, example
location:  
```
gs://quartz-bio/prd/repare_350001_bmdb/data/source/GENOSITY/ARCHER/SEQ_2111010396/5455505/FACETS/REPARE/
```

Example command line for running samples and writing results to quartz-bio:

```
nextflow -c snipdx.config run snipdx.nf -profile ref --metadata_file "path/to/metadata.tsv" --pon_version pon.092722 --output_to_quartz
```

The metadata file should have the following format (tab separated):

```
Sample_ID         Sample_Type      Bam                                   Bai
some_sample_1     FFPE             gs://bucket/some_sample_1.bam         gs://bucket/some_sample_1.bai
some_sample_2     blood            gs://bucket/some_sample_2.bam         gs://bucket/some_sample_2.bai
```

The metadata file can additionally include columns for `Trial` and `Inventory_File` with the command line flag 
`--rmarkdown_report` to automatically generate reports for all samples in the batch.

#### Report Processing Only

In cases where the samples need to be run through facets ahead of time before the data required for the reports is 
available, you can also run just the reporting section using the following metadata format:

```
Sample_ID         Facets_Folder                  Trial          Inventory_File
some_sample_1     gs://bucket/facets_out_1/      rp350001       gs://bucket/inventory.csv
some_sample_2     gs://bucket/facets_out_2/      nonclinical    gs://bucket/inventory.csv
```

to run the pipeline to only generate the report use the `--report_only` option, like below:

```
nextflow -c snipdx.config run snipdx.nf -profile ref --metadata_file "path/to/metadata.tsv" --experiment_name "test_experiment" --report_only 
```

#### Project level reports

This final stage of the pipeline should be executed after sample level facets have been ouputted for a new sequencing batch and stored to their respective trial and sample id in the quartz-bio google storage bucket.

This step compiles all sample level facets results and renders markdown reports partitioned by unique patients (for clinical samples) or project id (for non-clinical samples) specific to each trial.

For reproducibility purposes, a conda environment used to generate the project level rmarkdown reports can be replicated by importing the [rmarkdown.yml](rmarkdown.yml) file.

*Example command line for generating project level reports:*

```
## run if conda env does not exist
#conda env create --name R-4.0.5 -f environment.yml

conda activate R-4.0.5

trial="Non-Clinical"
inInventory="/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/INVENTORY/rp_nonclin_genosity_inv_annotated_2022-12-06.csv"
inSamplesFolder="/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/"
inVendorData="/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/vendor-data/gcp/"
processingDir="/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Non-Clinical/processing/"
in_nf_snipdx="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/nf-snipdx/"

Rscript /ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/nf-snipdx/generate_projects_report.R \
  --trial "$trial" \
  --inInventory "$inInventory" \
  --inSamplesFolder "$inSamplesFolder" \
  --inVendorData "$inVendorData" \
  --processingDir "$processingDir" \
  --in_nf_snipdx "$in_nf_snipdx"
```

*Individual command line parameters:*

| Parameters      | Description                                                                                                                                                                                                |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| trial           | trial number (choose from 'RP3500-01', 'RP3500-03' ,'Non-Clinical', 'RP6306-01', 'RP6306-02')                                                                                                              |
| processingDir   | specify path to working directory to export results and data folder ie) SNiPDx/[trial]/processing/                                                                                                         |
| inInventory     | specify path to latest genosity inventory file                                                                                                                                                             |
| inSamplesFolder | specify path to local trial specific Genosity facets output folder. If files were extracted from the `quartz-bio` bucket, then inSamplesFolder and inVendorData should point to the same folder            |
| inVendorDataw   | specify path to store local trial specific Genosity stats and variants folder. If files were extracted from the `quartz-bio` bucket, then inSamplesFolder and inVendorData should point ot the same folder |
| in_nf_snipdx    | specify path nf-snipdx repo folder (ensure current branch is set to master)                                                                                                                                |

------------------------------------------------------------------------------------------------

## Panel of Normal Versions


The panel of normal used for facets was updated due to a change in primers in the sequencing lab. After this change copy 
number artifacts appeared, primarily identified by the elevation of SETD2 targets. Normals were resequenced with the 
new primers and a new panel of normals was created. If older samples that were sequenced prior to the primer change need 
to be rerun they should be run with the old PON (pon.012222). New samples and samples from after the change should be 
run with the later PON (pon.092722).

### Sequencing runs to be run with pon.012222:


SEQ_2012080335  
SEQ_2101080007  
SEQ_2103150109  
SEQ_2105120174  
SEQ_2106010196  
SEQ_2107090249  

------------------------------------------------------------------------------------------------



