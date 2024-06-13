# Archer Panel Wiki Doc

This Wiki is copied from the Teams wiki and is left unmodified except for places that would cause confusion otherwise.
It is the original Archer wiki by Domonik and contains some explanations and guides on Archer processing and how to 
build the appropriate reference files for processing. 

# ArcherPanel Copy number analysis

Last edited: 4/2

_______________________________________________________________________________________________________________________
## Presentation summarizing the design of the panel
 
[Slides - concise introduction:](https://www.dropbox.com/s/o0rxm0gsdgb60o2/Panel.Overview.Jan2021.pptx?dl=0)
 
[December 2020](https://www.dropbox.com/s/q8i55hm5zrzxwgv/Panel.Overview.Dec2020.pptx?dl=0)
 
_______________________________________________________________________________________________________________________
## Overview of the analysis process

1) Identify germline SNPs (from Gnomad database) 
 
2) Analyze set of normal (control) samples. Perform pileup at SNP positions, using .bam files. Select SNPs that are 
appropriately covered (possibly stratify by sample QC). Only keep clean SNPs. Eliminate SNPs that are not heterozygous 
in normal samples, as this may be indicative of mapping problems.
 
3) On a new tumor sample, perform a pileup (at SNPs from 1). Simulate a tumor/normal input file from Facets, where normal 
coverage is derived from the coverage across normal samples.
 
4) Perform Facets analysis. For the tumor-only analysis, use Facet's tumor-only settings.
 
5) After facets, save a table specifying copy number (major, minor) in genomic segments, specified by chromosome, start 
and end positions. Make standard Facets plot and plots with different layout so they can be shown side-to-side with 
Ascat plots.

_______________________________________________________________________________________________________________________
## Previous Analysis of a new panel sample by Facets

### This section refers exclusively to the un-refactored Archer pipeline

Repository:
https://bitbucket.org/reparecompbio/osloh/src/master/
 
Bash script prepares a pileup at selected SNPs.
https://bitbucket.org/reparecompbio/osloh/src/master/src/script_run_single_sniprx_facets.sh
 
Estimate of purity, ploidy and copy number states by Facets (R script)
https://bitbucket.org/reparecompbio/osloh/src/master/src/run_sniprx_facets.R
 
Reference files:
list of SNPs condsidered: https://bitbucket.org/reparecompbio/osloh/src/master/data/external/gnomad2.genomes.and.exomes.400bp.dedup.vcf.gz
reference file (panel of normals): https://bitbucket.org/reparecompbio/osloh/src/master/data/external/Genosity.pon.19thJan.RData
both available in the Git repository and the relative paths should just work.
 
Dependencies
snp-pileup (I had installed with anaconda https://anaconda.org/bioconda/snp-pileup)
R ( I run 3.6.3)
facets R package (https://github.com/mskcc/facets)

_______________________________________________________________________________________________________________________
## Generation of SNP .vcf file
 
The file is used by 'snp-pileup' to specify which SNP locations to interrogate. Used for the pileup stage, and for 
everything else. The SNPs are taken from the Gnomad database, and only selected those in proximity of panel primers. 
The conservative threshold is 400bp.
 
https://bitbucket.org/reparecompbio/athenscn/src/master/src/facets/script_intersct_SNPs.sh

_______________________________________________________________________________________________________________________
## Generation of pileup files
 
"snp-pileup" in a loop over available .bam files. The format of bam files is as we get them from Genosity.
https://bitbucket.org/reparecompbio/athenscn/src/master/src/facets/script_run_facets_v2.sh
 
Visualization of a pileup file from version 2 of the panel
https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/facets/facetsExploreRaw_v2.ipynb?viewer=nbviewer

_______________________________________________________________________________________________________________________
## Panel of normals generation
 
[This notebook](https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/facets/facets_dominik_PON_v2_normals.ipynb?viewer=nbviewer) 
processes pileup files from normal samples, to generate a set of reliable SNPs and to estimate expected 
sequencing coverage at each. This reference data is used as reference when analyzing tumor samples.

_______________________________________________________________________________________________________________________
## Matched tumor-normal analysis

Tumor-normal analysis is what Facets had been designed for.
[See an example notebook here.](https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/facets/allFacetsPaired.ipynb?viewer=nbviewer)

_______________________________________________________________________________________________________________________
## Evaluation of the panel (version 2) calls
 
The accuracy of panel CNV calls is estimated by comparison to whole genome calls, 
 
[This notebook](https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/summaries/prepareSummary.ipynb?viewer=nbviewer) 
generates a [table](https://www.dropbox.com/s/u9gc3zufulst6bo/panel.v2.cnv.muts.v2.28Jan.csv?dl=0), (sample, gene) x (different copy number estimates). 
 
[This notebook](https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/summaries/prepareSummaryP2.ipynb?viewer=nbviewer) 
generates visual summaries of the tables, used in presentations.
 
[This notebook](https://bitbucket.org/reparecompbio/athenscn/src/master/notebooks/displayResults/displayWithPanel.ipynb?viewer=nbviewer) 
shows genome and panel visualizations - side to side.

_______________________________________________________________________________________________________________________

## Most recent dataset

 
Please see attached a table summarizing results from panel v2 with the following columns:

[panel.v2.cnv.muts.v2.April2nd.simplified.csv](https://www.dropbox.com/s/wtjvenp1lpih5i7/panel.v2.cnv.muts.v2.April2nd.simplified.csv?dl=0)

gene      

External.Sample.ID          

ger.id – Genosity ID

facets.purity – purity estimate from Facets

facets.ploidy – ploidy estimate from Facets

facets.V2.tcn.em – total copy number at locus, by Facets

facets.V2.lcn.em – minor copy number at locus, by Facets

Genosity.CNV.call – CNV call by Bob         

Chr. – chromosome of mutation (if any) 

Pos. – chromosomal position of mutations (if any)

hGVS.cDNA – cDNA position of a mutation            

hGVS.Protein – protein change   

FAF – allele fraction of mutation 

FDP – read depth at mutation position

has.WGS – whether the same sample also has WGS

validationUtility – LOH/amp/HomDel/SNV

Full table, with more columns that may turn out to be useful:

https://www.dropbox.com/s/knz93wh9dy519wp/panel.v2.cnv.muts.v2.April2nd.csv?dl=0

 

Please see a summary of WGS results

[call.df.cnvkit.genosity.ihc.wgs.muts.26Feb.simplified.csv](https://www.dropbox.com/s/cd54zc779glp785/call.df.cnvkit.genosity.ihc.wgs.muts.26Feb.simplified.csv?dl=0)

Gene

Sample, id.underscores – both are sample IDs

Ascat.Curation – whether the sample passed curation. ‘Exclude’ either because of FFPE damage to DNA or purity too low for meaningful analysis. The spreadsheet contains all samples, so do remove columns with ‘Exclude’ for more confident calls.

Purity – sample purity estimate from ASCAT

ascat.ploidy – cancer ploidy estimate from ASCAT

nTot – CNV call from ASCAT: total number of copies of gene         

nMinor – CNV call from ASCAT: minor copy number of gene (0 if LOH)

ascat.homdel,   

ascat.loh – conservative CNV calls derived from ASCAT

Genosity.call.wgs – CNV call made by Bob based on WGS data

panel.v2.genosity.curation.call   - CNV call from panel V2 (if sequenced)

wgs.CHROM - chromosome of a mutation in WGS(if any)

wgs.POS - chromosomal position of a mutation in WGS(if any)

wgs.X6 - consequence of a mutation (if any)

wgs.AF_tumor - allele fraction of a mutation (if any)

wgs.AD_tumor_alt - number of reads reporting a mutation (if any)

panel.V2Chr. - chromsome of a mutation identified by panel (if any)

panel.V2Pos. - position of a mutation identified by panel (if any)

panel.V2Ref. - reference alllele of a mutation

panel.V2Alt. - alternative allele of a mutation

panel.V2FAF - allele fraction of a mutation

panel.V2hGVS.cDNA - cDNA change caused by mutation

validationUtility - amp/HomDel/LOH/SNV

Full table with extra columns:

https://www.dropbox.com/s/7d3o8s418a1jnc5/call.df.cnvkit.genosity.ihc.wgs.muts.26Feb.csv?dl=0

Accessing the data from Genosity's cloud

 _______________________________________________________________________________________________________________________
## Genosity imports

s3://rpr-20-0001-0003/RPR-20-0004/

Panel v2 in the folder - transfer files with prefix   
RPR-20-0004/

Secret ID: AKIAQMOBNXCY45EDREE7

Passcode: (email Dominik for passcode)

For calculations, I copied the .bam files locally to my disk on the cloud:

/home/dglodzik/workdisk_link/VariantPlex/V2/genosityBams/


Copying the .bam files

cd /home/dglodzik/mounts/repare_disk_images

cp -r RPR-20-0004/ Projects/Athens/Genosity/athens_genosity_rpr-20-0001-0003/

/home/dglodzik/mounts/repare_disk_images/Projects/Athens/Genosity/athens_genosity_rpr-20-0001-0003/RPR-20-0004

cp -r RPR-20-0004 /home/dglodzik/workdisk_link/VariantPlex/V2/genosityBams/