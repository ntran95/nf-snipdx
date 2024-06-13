#!/usr/bin/env nextflow

/*
 * SNiPDX pipeline,
 * SNiPDX pipeline, takes snps from a bam and runs r-facets
 * Repare, copyright 2022
 * version 2.5
 */

nextflow.enable.dsl=1

def helpMessage() {
    log.info"""
    ========================
    SNiPDX NextFlow Pipeline
    ========================
    Description: Pipeline for facets analysis. When using a metadata file as input for batch mode, tab separated, in the following format:

    Sample_ID         Sample_Type    Bam                               Bai                               Trial          Inventory_File
    some_sample_1     blood          gs://bucket/some_sample_1.bam     gs://bucket/some_sample_1.bai     rp350001       gs://bucket/inventory1.csv
    some_sample_2     FFPE           gs://bucket/some_sample_2.bam     gs://bucket/some_sample_2.bai     nonclinical    gs://bucket/inventory2.csv

    Note: 'Trial' and 'Inventory_File' columns are only required if running with the pipeline with --rmarkdown_report, otherwise they should be excluded from the metadata file.
    To run just the report, use the flag --report_only and the following tab-separated format for the metadata file:

    Sample_ID         Facets_Folder                  Trial          Inventory_File
    some_sample_1     gs://bucket/facets_out_1/      rp350001       gs://bucket/inventory.csv
    some_sample_2     gs://bucket/facets_out_2/      nonclinical    gs://bucket/inventory.csv

    General usage: nextflow -c snipdx.config run snipdx.nf --metadata_file /path/to/snipdx_meta.tsv --experiment_name "RP3500-01_test_exp" -profile ref,gcp

    Required:
        --vcf_fp                       path to ref vcf file with targetted SNPs
        --panel_ga                     path to panel genes annotation file (tsv)
        --clust_info                   path to cluster info file (RData file)
        --pon_version                  version of panel of normals to be used ('pon.012222' or 'pon.092722')
        --metadata_file                path to the metadata file. Use this to specify multiple bams, including varying sample IDs and types

    Required when not outputting results to quartz-bio GCP bucket:
        --experiment_name              string of the experiment name for the batch, used to name the folder containing the cohort of sample results

    Required for PureCN:
        --purecn_normal_db             path to PureCN pon .rds file
        --purecn_interval_file         path to .txt file with PureCN intervals
        --purecn_interval_bed          path to .bed file with PureCN intervals
        --purecn_mapping_bias_file     path to mapping bias rds file
        --dbsnp                        path to dbsnp vcf.gz file
        --dbsnp_index                  path to dbsnp index .tbi file
        --ref_genome                   path to the reference genome
        --ref_genome_index             path to the reference genome index

    Required for local insert size normalization mode:
        --local_insertsize_norm        boolean, whether to normalize around local SNP insert sizes (default: false)
        --pon_target_location          path to the pon target location file for use in local insert size normalization

    Required for R Markdown report:
        --rmarkdown_report             boolean, whether to generate R Markdown report (default: false)
        --report_only                  boolean, whether to only generate the R Markdown report. Requires metadata file and experiment name (default: false)

    Optional:
        --output_folder                path to output the results to if not outputting to quartz-bio GCP bucket (defaults to present working directory)
        --output_to_quartz             boolean, whether to write the results to the quartz-bio GCP bucket (default: false)

    Optional PureCN Parameters:
        --purecn_min_af                float, minimum allele frequency to consider a variant (default: 0.05)
        --purecn_min_count             int, minimum number of reads in tumor and normal to keep an interval (default: 100)
        --purecn_seg_algorithm         string, segmentation algorithm for PureCN. One of 'CBS', 'PSCBS', 'Hclust' (default: 'PSCBS')
        --purecn_alpha                 float, significance of break points (default: 0.005)
        --purecn_min_purity            float, minimum considered purity (default: 0.1)
        --purecn_max_purity            float, maximum considered purity (default: 0.95)
        --purecn_min_ploidy            float, minimum considered ploidy (default: 1.4)
        --purecn_max_ploidy            float, maximum considered ploidy (default: 6)
        --purecn_max_ascn              int, maximum allele specific copy number (default: 6)
        --purecn_model                 string, model used to fit variants. One of 'beta' or 'betabin' (default: 'betabin')
        --purecn_max_nonclonal         float, maximum genomic fraction assigned to a subclonal copy number state (default: 0.3)

    Profiles:
        ref                            uses the path to local reference files included in this repo
        gcp                            use gcp containers and resources for computation for the repare-database project
        gcpClin                        outlines google life science connection parameters for the repare-clinical project
        gcpRef                         use reference files for purecn from their gcp location
        local                          use this profile for local compute. This will lower the system requirements and set gcp region parameters for local execution

    Trial Names:
        rp350001                       trial name and vendor data path based on rp3500-01 trial
        rp350003                       trial name and vendor data path based on rp3500-03 trial
        nonclinical                    trial name and vendor data path based on nonclinical samples
        rp630601                       trial name and vendor data path based on rp6306-01 trial
        rp630602                       trial name and vendor data path based on rp6306-02 trial
    """.stripIndent()
}

// command line parameters
params.output_folder                  = false
params.help                           = false
params.metadata_file                  = false
params.experiment_name                = false
params.vcf_fp                         = false
params.pon_fp_012222                  = false
params.pon_fp_092722                  = false
params.panel_ga                       = false
params.clust_info                     = false
params.pon_version                    = false
params.local_insertsize_norm          = false
params.pon_target_location            = false
params.rmarkdown_report               = false
params.report_only                    = false
params.snipdx_gene_coords             = false
params.ref_genome                     = false
params.ref_genome_index               = false
params.dbsnp                          = false
params.dbsnp_index                    = false
params.purecn_normal_db               = false
params.purecn_interval_file           = false
params.purecn_interval_bed            = false
params.purecn_mapping_bias_file       = false
params.purecn_min_af                  = 0.05
params.purecn_min_count               = 100
params.purecn_seg_algorithm           = 'PSCBS'
params.purecn_alpha                   = 0.005
params.purecn_min_purity              = 0.1
params.purecn_max_purity              = 0.95
params.purecn_min_ploidy              = 1.4
params.purecn_max_ploidy              = 5
params.purecn_max_ascn                = 6
params.purecn_model                   = 'betabin'
params.purecn_max_nonclonal           = 0.3
params.output_to_quartz               = false

// display help message if requested then exit
if(params.help) {
    helpMessage()
    exit 0
}

// Show help message if min required components are not present
batch_mode = !( (!params.metadata_file || !params.experiment_name) && params.output_to_quartz)
if(!batch_mode && !params.output_to_quartz) {
    helpMessage()
    println('Error! Missing one or more required command parameters! Exiting.')
    exit 1
}

// check for parameters required for PureCN
if(!params.ref_genome || !params.ref_genome_index || !params.dbsnp || !params.dbsnp_index || !params.purecn_normal_db || !params.purecn_interval_file || !params.purecn_mapping_bias_file) {
    helpMessage()
    println('Error! Missing 1 or more required reference files for PureCN analysis!')
    exit 1
}

// report only mode check
if(params.report_only && (!params.metadata_file || !params.experiment_name)){
    helpMessage()
    println('Error! Running report only mode but missing either the metadata or experiment name parameters! Exiting.')
    exit 1
}

// set globals val for whether we are running a report and whether we are running facets
run_facets      = !params.report_only
generate_report = params.rmarkdown_report || params.report_only

// trial information contained here. This is all hardcoded here, sorry
def acceptable_trial_name_vals = ['rp350001', 'rp350003' ,'nonclinical', 'rp630601', 'rp630602']
def trial_name_long = [rp350001: 'rp3500-01', rp350003: 'rp3500-03',nonclinical: 'non-clinical', rp630601: 'rp6306-01', rp630602: 'rp6306-02']
def trial_vendor_basepaths = [rp350001: 'gs://quartz-bio/prd/repare_350001_bmdb/data/source/GENOSITY/ARCHER/',
                             rp350003: 'gs://quartz-bio/prd/repare_350003_bmdb/data/source/GENOSITY/ARCHER/',
                             nonclinical: 'gs://quartz-bio/prd/repare_nonclin_bmdb/data/source/GENOSITY/ARCHER/',
                             rp630601: 'gs://quartz-bio/prd/repare_rp630601_bmdb/data/source/GENOSITY/ARCHER/',
                             rp630602: 'gs://quartz-bio/prd/repare_rp630602_bmdb/data/source/GENOSITY/ARCHER/']

// set the acceptable sample type vals:
def acceptable_sample_type_vals = ['blood', 'FFPE']


// pull in utils scripts and SNiPDX facets script
run_sniprx_facets_script_file          = file("${baseDir}/src/run_sniprx_facets.R")
draw_raw_sample_script_file            = file("${baseDir}/src/utils/drawRawSample.R")
draw_raw_sample_gene_script_file       = file("${baseDir}/src/utils/drawRawSampleGenes.R")
plot_facets_script_file                = file("${baseDir}/src/utils/plot.facets.R")
plot_IS_script_file                    = file("${baseDir}/src/utils/plotIS.R")
classify_sample_script_file            = file("${baseDir}/src/utils/classifySample.R")
get_sample_primer_distance_script_file = file("${baseDir}/src/utils/getSamplePrimerDistance.R")
normalize_wrt_primer_script_file       = file("${baseDir}/src/utils/normalizeWrtPrimer.R")
plot_sample_modified_script_file       = file("${baseDir}/src/utils/plotSampleModified.R")
find_closest_samples_script_file       = file("${baseDir}/src/utils/findClosestSamples.R")
plot_cover_scatterplot_script_file     = file("${baseDir}/src/utils/plotCoverageScatterplot.R")
logRlog_or_spider_script_file          = file("${baseDir}/src/utils/logRlogORspider.mod.R")
write_bw_script_file                   = file("${baseDir}/src/utils/writeBw.R")
plot_gc_script_file                    = file("${baseDir}/src/utils/plotGC.R")
extract_ei_bias_script_file            = file("${baseDir}/src/utils/extract.ei.bias.R")
plot_IS_analysis_file                  = file("${baseDir}/src/utils/plotISanalysis.R")
normalize_by_is_file                   = file("${baseDir}/src/utils/normalize_by_is.R")
get_chrom_plots_and_segments_file      = file("${baseDir}/src/utils/get_chrom_plots_and_segments.py")
get_diplogr_from_purecn_file           = file("${baseDir}/src/utils/get_diplogr_from_purecn.py")

// read in reference genome fa file
ref_genome_ch = load_and_check_file(params.ref_genome)

// read in reference genome index fai file
ref_genome_index_ch = load_and_check_file( params.ref_genome_index)

// purecn params
purecn_normal_db_ch               = file(params.purecn_normal_db)
purecn_interval_file_ch           = file(params.purecn_interval_file)
purecn_interval_bed_ch            = file(params.purecn_interval_bed)
purecn_mapping_bias_file_ch       = file(params.purecn_mapping_bias_file)
dbsnp_ch                          = file(params.dbsnp)
dbsnp_index_ch                    = file(params.dbsnp_index)
purecn_min_af_ch                  = params.purecn_min_af
purecn_min_count_ch               = params.purecn_min_count
purecn_seg_algorithm_ch           = params.purecn_seg_algorithm
purecn_alpha_ch                   = params.purecn_alpha
purecn_min_purity_ch              = params.purecn_min_purity
purecn_max_purity_ch              = params.purecn_max_purity
purecn_min_ploidy_ch              = params.purecn_min_ploidy
purecn_max_ploidy_ch              = params.purecn_max_ploidy
purecn_max_ascn_ch                = params.purecn_max_ascn
purecn_model_ch                   = params.purecn_model
purecn_max_nonclonal_ch           = params.purecn_max_nonclonal


// check the vcf file exists
vcf_fp = load_and_check_file(params.vcf_fp)

// check the cluster info file exists
clust_info = load_and_check_file(params.clust_info)

// check the gene annotation file exists
panel_ga = load_and_check_file(params.panel_ga)

// check snipdx gene coords file exists
snipdx_gene_coords = load_and_check_file(params.snipdx_gene_coords)

// if running local insert size normalization load these files
if(params.local_insertsize_norm){
    // check the insert targets file exists
    pon_target_locations_file = load_and_check_file(params.pon_target_location)
    calc_IS_at_targets_script_file = file("${baseDir}/src/calcIsAtTargets.py")
}else{
    // if we aren't running insert size normalization we need to load a dummy file for
    // input into the facets process
    // for more details see the old NextFlow patterns documentation:
    // https://nextflow-io.github.io/patterns/index.html#_optional_input
    Channel.fromPath("${baseDir}/src/utils/NO_FILE").set{ dummy_file_ch }
}

// get correct pon version
if(params.pon_version == 'pon.012222') {
    pon_ch = file(params.pon_fp_012222)
} else if(params.pon_version == 'pon.092722') {
    pon_ch = file(params.pon_fp_092722)
} else {
    println('Error! pon_version parameter value not recognized. Please specify one of "pon.012222" or "pon.092722"')
    exit 1
}

// ensure that no output_folder parameter is given when option to write to quartz-bio is used
if(params.output_folder && params.output_to_quartz) {
    println('No output_folder argument should be given when writing to quartz-bio bucket')
    exit 1
}

// read in metadata file
if ( run_facets ) {
    //  read in the metadata file
    // format for metadata validation lovingly stolen from nf-core at
    // https://github.com/nf-core/eager/blob/master/main.nf
    Channel
        .fromPath(params.metadata_file, checkIfExists: true)
        .splitCsv(header: true, sep:'\t')
        .map{ row->
            // expected columns dependant on if we are generating a report
            def expected_keys = ['Sample_ID', 'Sample_Type', 'Bam', 'Bai']
            if (params.rmarkdown_report){
                expected_keys = ['Sample_ID', 'Sample_Type', 'Bam', 'Bai', 'Trial', 'Inventory_File']
            }

            // check all headers are present
            if ( !row.keySet().containsAll(expected_keys) ) exit 1, "Error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

            checkNumberOfItem(row, expected_keys.size())

            // make sure no elements are empty
            if ( row.Sample_ID == null || row.Sample_ID.isEmpty() ) exit 1, "Error: the Sample_ID column is empty. Check row:\n ${row}"
            if ( row.Sample_Type == null || row.Sample_Type.isEmpty() ) exit 1, "Error: the Sample_Type column is empty. Check row:\n ${row}"
            if ( row.Bam == null || row.Bam.isEmpty() ) exit 1, "Error: the Bam column is empty. Check row:\n ${row}"
            if ( row.Bai == null || row.Bai.isEmpty() ) exit 1, "Error: the Bai column is empty. Check row:\n ${row}"

            // set row vals to variables
            def sample_name_val = row.Sample_ID
            def sample_type_val = row.Sample_Type
            def bam_file_val = row.Bam
            def bai_file_val = row.Bai
            def trial_name = params.rmarkdown_report ? row.Trial : null

            // parse bam file name to get path to output to quartz-bio GCP bucket
            quartz_path = (bam_file_val =~ /gs:\/\/quartz-bio\/(.+)BAM.+/)[0][1] + "FACETS/REPARE"

            // check no empty metadata fields
            if (sample_name_val == '' || sample_type_val == '' || bam_file_val == '' || bai_file_val == '' ) exit 1, "Error: a field/column does not contain any information. Ensure all cells are filled. Check row:\n ${row}"

            // Check sample type is valid
            if (!acceptable_sample_type_vals.contains(sample_type_val)) exit 1, "Error: Sample Type in TSV must be one of ${acceptable_sample_type_vals.join(", ")}. Check row:\n ${row}"

            // So we don't accept existing files that are wrong format
            if (!has_extension(bam_file_val, "bam")) exit 1, "Error: A specified Bam file either has a non-recognizable BAM extension. Check row:\n ${row}"
            if (!has_extension(bai_file_val, "bai")) exit 1, "Error: A specified Bai file either has a non-recognizable BAI extension. Check row:\n ${row}"

            // row passes basic checks, output it into the channel
            tuple(sample_name_val, sample_type_val, tuple(file(bam_file_val), file(bai_file_val)), quartz_path) }
        .into { bam_channel_1; bam_channel_2; bam_channel_3; bam_channel_4; bam_channel_5; bam_channel_6 }
}

if(generate_report && !run_facets){
    // here we can assume we are not running facets and doing a report only mode
    // run again building another channel based on the metadata input for the r markdown if we are
    // generating a report. We assume all other columns were tested reading it in the previous time
    if (generate_report){
        Channel
        .fromPath(params.metadata_file, checkIfExists: true)
        .splitCsv(header: true, sep:'\t')
        .map{ row->
            // expected columns dependant on if we are generating a report
            def expected_keys = ['Sample_ID', 'Facets_Folder', 'Trial', 'Inventory_File']

            // check all headers are present
            if ( !row.keySet().containsAll(expected_keys) ) exit 1, "Error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

            checkNumberOfItem(row, expected_keys.size())

            if (row.Sample_ID == null || row.Sample_ID.isEmpty()) exit 1, "Error: The Sample_ID column is empty. Check row:\n ${row}"
            if (row.Facets_Folder == null || row.Facets_Folder.isEmpty()) exit 1, "Error: The Facets_Folder column is empty. Check row:\n ${row}"

            // set row vals to variables
            def sample_name_val = row.Sample_ID
            def facets_folder = row.Facets_Folder

            if (sample_name_val == '' || facets_folder == '') exit 1, "Error: a field/column does not contain any information. Ensure all cells are filled. Check row:\n ${row}"

            // read in the files from the facets output folder
            if (facets_folder[-1] != "/"){
                facets_folder_regex = facets_folder + '/*'
            } else{
                facets_folder_regex = facets_folder + '*'
            }
            facets_files = file(facets_folder_regex)

            // row passes basic checks, output it into the channel
            tuple(sample_name_val, facets_files) }
        .set { facet_script_output_channel }
    }
}

// run again building another channel based on the metadata input for the r markdown if we are
// generating a report. We assume all other columns were tested reading it in the previous time
if ( generate_report && run_facets ){
    Channel
    .fromPath(params.metadata_file, checkIfExists: true)
    .splitCsv(header: true, sep:'\t')
    .map{ row->
        if (row.Trial == null || row.Trial.isEmpty()) exit 1, "Error: The Trial column is empty. Check row:\n ${row}"
        if (row.Inventory_File == null || row.Inventory_File.isEmpty()) exit 1, "Error: The Inventory File column is empty. Check row:\n ${row}"

        // set row vals to variables
        def sample_name_val = row.Sample_ID
        def trial_name = row.Trial
        def inventory_file = row.Inventory_File

        if (trial_name == '' || inventory_file == '' || sample_name_val == '') exit 1, "Error: a field/column does not contain any information. Ensure all cells are filled. Check row:\n ${row}"

        // make sure trial type is recognized if we are doing trials
        if (!acceptable_trial_name_vals.contains(trial_name)) exit 1, "Error: Trial Name in TSV must be one of ${acceptable_trial_name_vals.join(", ")}. Check row:\n ${row}"

        // get batch_id from sample_name if applicable
        // this will re-evaluate every loop but in cases where we use it
        // to determine the vendor data path for the r markdown report
        // (which should be its only use), we should
        // assume it is the same for all samples within the cohort
        batch_id = get_batch_id_from_sample_name(sample_name_val)

        vendor_data_base_path   = trial_vendor_basepaths[trial_name]
        trial_name_long_val     = trial_name_long[trial_name]
        // vendor data path is inferred here
        vendor_data_path        = "${vendor_data_base_path}${batch_id}/RESULTS/"

        vendor_cnv_file         = "${vendor_data_path}${batch_id}.CNV.tsv"
        vendor_stats_file       = "${vendor_data_path}${batch_id}.Stats_Summary.tsv"
        vendor_variant_file     = "${vendor_data_path}${batch_id}.variant_summary.tsv"

        // row passes basic checks, output it into the channel
        tuple(sample_name_val, trial_name_long_val, file(inventory_file), file(vendor_cnv_file), file(vendor_stats_file), file(vendor_variant_file)) }.view{"$it"}
    .set { rmarkdown_files_ch }
}


// this section will deal with R markdown variables
// for cases where we are running the R markdown report script
if(generate_report){
    // pull in R markdown report script and dependant files
    rmarkdown_report_script_ch = file("${baseDir}/src/generate_SNiPDx_clinical_by_sample.Rmd")
}

// create channel to output metadata in run_info
if (batch_mode){
   run_info_metadata_file_ch = file(params.metadata_file)
}

// check if we have an output dir, otherwise we set it to the current working directory or to prefix of quartz-bio
// bucket if --output_to_quartz true
if (params.output_to_quartz) {
    output_dir = "gs://quartz-bio/"
} else if ( params.output_folder ) {
    output_dir = "${params.output_folder}/${params.experiment_name}/"
} else {
    output_dir = "$PWD/${params.experiment_name}/"
}

// get GIT commit sha
try {
   git_commit_sha = new File("${projectDir}/.git/refs/heads/master").text.replaceAll("[\n\r]", "")
} catch(Exception e) {
   println("Error retrieving the git commit SHA!")
   git_commit_sha = "N/A"
}

/*
 *
 * Auxiliary Functions
 *
 */

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Error:  Invalid TSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag for more information"
    return true
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// check file exists and return it if it does, otherwise exit with error
def load_and_check_file(file_path) {
    file_to_check = file(file_path)
    if (!file_to_check.exists())  exit 1, "Error: Input file does not exist, check path at $file_path"
    return file_to_check
}


// remove the batch_id from the sample_name. We can assume the
// batch_id is the last section of the sample_name, and consists of the following format:
// {sample_name}_{prefix}_{batch_number}, meaning we want to return the '{prefix}_{batch_number}' portion.
def get_batch_id_from_sample_name(sample_name) {
    sample_name_split = "${sample_name}".split('_')
    batch_id = sample_name_split[-2] + "_" + sample_name_split[-1]
    return batch_id
}

/*
 *
 * Runtime Parameters Log
 *
 */

// print run params to stdout
run_cmd_line_info = """Run Parameters:
==================================================
Metadata File                : ${params.metadata_file}
Experiment Name              : ${params.experiment_name}
BAM File                     : ${params.bam_file}
BAM Index File               : ${params.bai_file}
VCF File                     : ${params.vcf_fp}
PON Version                  : ${params.pon_version}
Panel Gene Annotation File   : ${params.panel_ga}
Clust Info File              : ${clust_info}
Insert Size Normalization    : ${params.local_insertsize_norm}
Pon Target Location File     : ${params.pon_target_location}
Generate R MarkDown Report   : ${generate_report}
Run Facets                   : ${run_facets}
Reference Genome             : ${params.ref_genome}
Reference Genome Index       : ${params.ref_genome_index}
PureCN PON                   : ${params.purecn_normal_db}
PureCN Interval Text File    : ${params.purecn_interval_file}
PureCN Interval Bed File     : ${params.purecn_interval_bed}
PureCN Mapping Bias File     : ${params.purecn_mapping_bias_file}
dbSNP Vcf                    : ${params.dbsnp}
dbSNP Vcf Index              : ${params.dbsnp_index}
PureCN Min VAF               : ${params.purecn_min_af}
PureCN Min Variant Count     : ${params.purecn_min_count}
PureCN Segmentation Method   : ${params.purecn_seg_algorithm}
PureCN Segmentation Alpha    : ${params.purecn_alpha}
PureCN Min Purity            : ${params.purecn_min_purity}
PureCN Max Purity            : ${params.purecn_max_purity}
PureCN Min Ploidy            : ${params.purecn_min_ploidy}
PureCN Max Ploidy            : ${params.purecn_max_ploidy}
PureCN Max ASCN              : ${params.purecn_max_ascn}
PureCN Variant Model         : ${params.purecn_model}
PureCN Max Nonclonal         : ${params.purecn_max_nonclonal}
Output Folder                : ${output_dir}
Output to Quartz             : ${params.output_to_quartz}
Pipeline GIT Commit          : ${git_commit_sha}
Command                      : $workflow.commandLine
Manifest Pipeline Ver        : $workflow.manifest.version
Nextflow Ver                 : $nextflow.version
==================================================
"""
log.info("${run_cmd_line_info}")


/*
 *
 * Finished reading in inputs, run the pipeline
 *
 */

if (run_facets){

    // run freebayes to get variant calls used in PureCN modeling
    process freebayes {

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/freebayes/$filename"
            else "${sample_name}/freebayes/$filename"
        }

        input:
        set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_4
        file(ref_genome) from ref_genome_ch
        file(ref_genome_index) from ref_genome_index_ch

        output:
        tuple sample_name, file("${sample_name}.freebayes.vcf"), quartz_path into freebayes_out_ch
        file(".command.log") into freebayes_log_out_ch

        script:
        """
        touch *.bai
        touch *.fai
        freebayes -b ${bam_files[0]} -f $ref_genome -v ${sample_name}.freebayes.vcf \
            -r chr1 \
            -r chr2 \
            -r chr3 \
            -r chr4 \
            -r chr5 \
            -r chr6 \
            -r chr7 \
            -r chr8 \
            -r chr9 \
            -r chr10 \
            -r chr11 \
            -r chr12 \
            -r chr13 \
            -r chr14 \
            -r chr15 \
            -r chr16 \
            -r chr17 \
            -r chr18 \
            -r chr19 \
            -r chr20 \
            -r chr21 \
            -r chr22 \
            -r chrX \
            -r chrY \
            --min-mapping-quality 19 \
            --min-alternate-count 5 \
            --min-alternate-fraction 0.05 \
            --min-coverage 50 \
            --min-base-quality 20
        """
    }

    // annotate variants with dbsnp ID to infer germline variants for PureCN modeling
    process dbsnp_annotate_variants {

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/freebayes/$filename"
            else "${sample_name}/freebayes/$filename"
        }

        input:
        set sample_name, file(vcf), quartz_path from freebayes_out_ch
        file(dbsnp) from dbsnp_ch
        file(dbsnp_index) from dbsnp_index_ch

        output:
        tuple sample_name, file("${sample_name}.freebayes.dbsnp_annotated.vcf") into freebayes_annotate_out_ch
        file(".command.log") into dbsnp_log_out_ch

        script:
        """
        touch *.tbi
        java -jar -Xmx64g /snpEff/SnpSift.jar annotate \
            $dbsnp \
            $vcf > "${sample_name}.freebayes.dbsnp_annotated.vcf"
        """
    }

    // get sample coverage per PureCN intervals
    process purecn_coverage {

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/purecn/$filename"
            else "${sample_name}/purecn/$filename"
        }

        input:
        set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_5
        file(intervals) from purecn_interval_file_ch

        output:
        tuple sample_name, file("coverage/${sample_name}.final_coverage.txt.gz"), quartz_path into coverage_out_ch
        file("coverage/.command.log") into coverage_log_out_ch

        script:
        """
        touch *.bai
        mkdir coverage
        Rscript /opt/PureCN/Coverage.R \
            --out-dir coverage/ \
            --bam ${bam_files[0]} \
            --intervals $intervals \
            --remove-mapq0 \
            --cores $task.cpus \
            --skip-gc-norm
        cp .command.log coverage/
        """
    }

    purecn_input_ch = freebayes_annotate_out_ch.join(coverage_out_ch)

    // run PureCN
    process purecn {

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/purecn/$filename"
            else "${sample_name}/purecn/$filename"
        }

        input:
        set sample_name, file(vcf), file(coverage), quartz_path from purecn_input_ch
        file(normal_db) from purecn_normal_db_ch
        file(mapping_bias) from purecn_mapping_bias_file_ch
        file(intervals) from purecn_interval_file_ch
        val(seg_algorithm) from purecn_seg_algorithm_ch
        val(model) from purecn_model_ch
        val(max_nonclonal) from purecn_max_nonclonal_ch
        val(min_af) from purecn_min_af_ch
        val(min_counts) from purecn_min_count_ch
        val(alpha) from purecn_alpha_ch
        val(min_purity) from purecn_min_purity_ch
        val(max_purity) from purecn_max_purity_ch
        val(min_ploidy) from purecn_min_ploidy_ch
        val(max_ploidy) from purecn_max_ploidy_ch
        val(max_ascn) from purecn_max_ascn_ch

        output:
        file("results/*") into purecn_out_ch
        tuple sample_name, file("results/*.rds") into purecn_rds_out_ch
        tuple sample_name, quartz_path, file("results/${sample_name}_loh.csv") into purecn_seg_out_ch
        file("results/.command.log") into purecn_log_out_ch

        script:
        """
        mkdir results
        Rscript /opt/PureCN/PureCN.R \
            --sampleid $sample_name \
            --out results \
            --tumor $coverage \
            --vcf $vcf \
            --normaldb $normal_db \
            --mapping-bias-file $mapping_bias \
            --intervals $intervals \
            --genome hg19 \
            --fun-segmentation $seg_algorithm \
            --model $model \
            --max-non-clonal $max_nonclonal \
            --min-af $min_af \
            --min-total-counts $min_counts \
            --alpha $alpha \
            --min-purity $min_purity \
            --max-purity $max_purity \
            --min-ploidy $min_ploidy \
            --max-ploidy $max_ploidy \
            --max-copy-number $max_ascn \
            --cores $task.cpus \
            --post-optimize --seed 123
        cp .command.log results/
        """
    }

    // get the mean log2 ratio of segments where C=2 from PureCN output; will be used as argument to facets
    process diplogr_calculation {

        label 'teeny_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/purecn/$filename"
            else "${sample_name}/purecn/$filename"
        }

        input:
        set sample_name, quartz_path, file(purecn_seg_file) from purecn_seg_out_ch
        file(get_diplogr_from_purecn_script) from get_diplogr_from_purecn_file

        output:
        tuple sample_name, file("${sample_name}_purecn_diplogr.txt") into diplogr_out_ch

        script:
        """
        python3 $get_diplogr_from_purecn_script $purecn_seg_file $sample_name
        """
    }

    // get callable regions from PureCN intervals based on quality + coverage thresholds, used for PureCN biomarkers
    process callstate {

        label 'small_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/purecn/callstate/$filename"
            else "${sample_name}/purecn/callstate/$filename"
        }

        input:
        set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_6
        file(intervals) from purecn_interval_bed_ch

        output:
        tuple sample_name, file("${sample_name}_callable_edited.bed"), quartz_path into callable_out_ch
        file(".command.log") into callstate_log_out_ch

        script:
        """
        touch *.bai
        callstate $intervals ${bam_files[0]} -o "${sample_name}_callable_status.bed" --min-mapq 19 --min-base-qual 20 --min-depth 100
        grep CALLABLE "${sample_name}_callable_status.bed" > "${sample_name}_callable.bed"
        awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4}' "${sample_name}_callable.bed" > tmp && mv tmp "${sample_name}_callable_edited.bed"
        """
    }

    biomarker_input_ch = purecn_rds_out_ch.join(callable_out_ch)

    // run PureCN biomarkders script to output CIN,TMB,mutational signatures; TMB + signatures may not be available
    // or accurate due to small target region
    process purecn_biomarkers {

        label "reg_task"

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/purecn/biomarkers/$filename"
            else "${sample_name}/purecn/biomarkers/$filename"
        }

        input:
        set val(sample_name), file(purecn_rds), file(callable_bed), quartz_path from biomarker_input_ch

        output:
        file("*.csv") into biomarker_out_ch
        file(".command.log") into biomarkers_log_out_ch

        script:
        """
        Rscript /opt/PureCN/Dx.R --out $sample_name --rds $purecn_rds --callable $callable_bed --signatures
        """
    }

    // create pileup from bam file
    process snp_pileup{

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/pileup/$filename"
            else "${sample_name}/pileup/$filename"
        }

        input:
        set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_1
        file vcf_fp from vcf_fp

        output:
        tuple val("${sample_name}"), val("${sample_type}"), file("${sample_name}.snp_pileup*") into pileup_output_channel
        file(".command.log") into pileup_log_ch

        script:
        """
        snp-pileup -A --min-map-quality=15 --min-base-quality=15 --gzip --max-depth=15000 \
        "${vcf_fp}" \
        "${sample_name}.snp_pileup" \
        "${bam_files[0]}"
        """
    }

    // get insert metrics
    process insert_metrics{

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/insert_metrics/$filename"
            else "${sample_name}/insert_metrics/$filename"
        }

        input:
        set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_2

        output:
        tuple val("${sample_name}"), file ("${sample_name}.final.bam.stats"), file ("${sample_name}.final.bam.insert.stats"), quartz_path into insert_stats_file_ch
        file(".command.log") into insert_metrics_log_ch

        script:
        """
        samtools stats --remove-dups ${bam_files[0]} > "${sample_name}.final.bam.stats"
        cat "${sample_name}.final.bam.stats" | grep ^IS | cut -f 2- > "${sample_name}.final.bam.insert.stats"
        """
    }

    if(params.local_insertsize_norm){
        // calculate local insert size metrics
        process local_insert_metrics{

            label 'reg_task'

            publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
                if (params.output_to_quartz) "${quartz_path}/local_insert_metrics/$filename"
                else "${sample_name}/local_insert_metrics/$filename"
            }

            input:
            set sample_name, sample_type, file(bam_files), quartz_path from bam_channel_3
            file pon_target_locations from pon_target_locations_file
            file calc_IS_at_targets_script from calc_IS_at_targets_script_file

            output:
            tuple val("${sample_name}"), file("${sample_name}.insert.at.targets") into local_insert_metric_out_ch
            file(".command.log") into local_insert_metrics_log_ch

            script:
            """
            python ${calc_IS_at_targets_script} \
                        --targets_csv ${pon_target_locations} \
                        --bam "${bam_files[0]}" \
                        --output_fp "${sample_name}.insert.at.targets"
            """
        }

        // merge pileup and insert metrics and local insert metrics output into 1 channel for gather step
       pileup_output_channel.combine(insert_stats_file_ch, by:0).combine(local_insert_metric_out_ch, by:0).into { facets_input_tuple_channel_1; facets_input_tuple_channel_2 }
    }else{
        // merge pileup and insert metrics output into 1 channel for gather step
        // We also need to add a dummy file to keep tuple ordinality going into the facets process
        // for more details see the old NextFlow patterns documentation:
        // https://nextflow-io.github.io/patterns/index.html#_optional_input
        pileup_output_channel.combine(insert_stats_file_ch, by:0).combine( dummy_file_ch ).into { facets_input_tuple_channel_1; facets_input_tuple_channel_2 }
    }

    // run main facets script
    process run_facet_script{

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/facets/$filename"
            else "${sample_name}/facets/$filename"
        }

        input:
        // sample files/info
        set sample_name, sample_type, file(snp_pileup), file(samtools_stat), file(insert_metrics), quartz_path, file(local_insert_metrics) from facets_input_tuple_channel_1
        // ref files
        file pon_fp from pon_ch
        file panel_ga from panel_ga
        file clust_info from clust_info
        // import sniprx main script
        file run_sniprx_facets_script from run_sniprx_facets_script_file
        // import util scripts
        file draw_raw_sample_script from draw_raw_sample_script_file
        file draw_raw_sample_gene_script from draw_raw_sample_gene_script_file
        file plot_facets_script from plot_facets_script_file
        file plot_IS_script from plot_IS_script_file
        file classify_sample_script from classify_sample_script_file
        file get_sample_primer_distance_script from get_sample_primer_distance_script_file
        file normalize_wrt_primer_script from normalize_wrt_primer_script_file
        file plot_sample_modified_script from plot_sample_modified_script_file
        file find_closest_samples_script from find_closest_samples_script_file
        file plot_cover_scatterplot_script from plot_cover_scatterplot_script_file
        file logRlog_or_spider_script from logRlog_or_spider_script_file
        file write_bw_script from write_bw_script_file
        file plot_gc_script from plot_gc_script_file
        file extract_ei_bias_script from extract_ei_bias_script_file
        file plot_IS_analysis from plot_IS_analysis_file
        file normalize_by_is from normalize_by_is_file

        output:
        tuple val("${sample_name}"), quartz_path, file("*") into facet_script_output_channel_1, facet_script_output_channel_2
        file(".command.log") into facets_log_ch

        script:
        // start by setting the local insert metric normalization params if we are using them:
        def insert_size_norm = local_insert_metrics.name != 'NO_FILE' ? "--normalize_is TRUE --sample_is_fn ${local_insert_metrics}" : "--normalize_is FALSE"
        // not run the facets R script
        """
        Rscript ${run_sniprx_facets_script} \
            --pileup_fn "${snp_pileup}" \
            --sample_id "${sample_name}" \
            --sample_output_folder "${sample_name}/" \
            --sample_insert_distribution "${insert_metrics}" \
            --min_cov 500 \
            --het_threshold 0.05\
            --secondary_normalization TRUE \
            --cval 50 \
            --deltaCN 0.2 \
            --sample_type "${sample_type}" \
            --snp.nbhd 250 \
            --only_hets FALSE \
            --reference_PON "${pon_fp}" \
            --panel_genes_file "${panel_ga}" \
            --cluster_info "${clust_info}" \
            --utils_path "" ${insert_size_norm} \
            --diplogr_file "NA"
        # move all output files to top directory then delete the directory
        mv ${sample_name}/* .
        rm -rf ${sample_name}/
        """
    }

    facets_w_diplogr_input_ch = facets_input_tuple_channel_2.join(diplogr_out_ch)

    // run main facets script with diplogr from purecn
    process run_facet_script_with_purecn_diplogr {

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/facets_with_purecn_diplogr/$filename"
            else "${sample_name}/facets_with_purecn_diplogr/$filename"
        }

        input:
        // sample files/info
        set sample_name, sample_type, file(snp_pileup), file(samtools_stat), file(insert_metrics), quartz_path, file(local_insert_metrics), file(diplogr) from facets_w_diplogr_input_ch
        // ref files
        file pon_fp from pon_ch
        file panel_ga from panel_ga
        file clust_info from clust_info
        // import sniprx main script
        file run_sniprx_facets_script from run_sniprx_facets_script_file
        // import util scripts
        file draw_raw_sample_script from draw_raw_sample_script_file
        file draw_raw_sample_gene_script from draw_raw_sample_gene_script_file
        file plot_facets_script from plot_facets_script_file
        file plot_IS_script from plot_IS_script_file
        file classify_sample_script from classify_sample_script_file
        file get_sample_primer_distance_script from get_sample_primer_distance_script_file
        file normalize_wrt_primer_script from normalize_wrt_primer_script_file
        file plot_sample_modified_script from plot_sample_modified_script_file
        file find_closest_samples_script from find_closest_samples_script_file
        file plot_cover_scatterplot_script from plot_cover_scatterplot_script_file
        file logRlog_or_spider_script from logRlog_or_spider_script_file
        file write_bw_script from write_bw_script_file
        file plot_gc_script from plot_gc_script_file
        file extract_ei_bias_script from extract_ei_bias_script_file
        file plot_IS_analysis from plot_IS_analysis_file
        file normalize_by_is from normalize_by_is_file

        output:
        tuple val("${sample_name}"), quartz_path, file("*") into facet_w_purecn_diplogr_output_channel
        file(".command.log") into facets_w_purecn_diplogr_log_ch

        script:
        // start by setting the local insert metric normalization params if we are using them:
        def insert_size_norm = local_insert_metrics.name != 'NO_FILE' ? "--normalize_is TRUE --sample_is_fn ${local_insert_metrics}" : "--normalize_is FALSE"
        // now run the facets R script
        """
        Rscript ${run_sniprx_facets_script} \
            --pileup_fn "${snp_pileup}" \
            --sample_id "${sample_name}" \
            --sample_output_folder "${sample_name}/" \
            --sample_insert_distribution "${insert_metrics}" \
            --min_cov 500 \
            --het_threshold 0.05\
            --secondary_normalization TRUE \
            --cval 50 \
            --deltaCN 0.2 \
            --sample_type "${sample_type}" \
            --snp.nbhd 250 \
            --only_hets FALSE \
            --reference_PON "${pon_fp}" \
            --panel_genes_file "${panel_ga}" \
            --cluster_info "${clust_info}" \
            --utils_path "" ${insert_size_norm} \
            --diplogr_file "${diplogr}"
        # move all output files to top directory then delete the directory
        mv ${sample_name}/* .
        rm -rf ${sample_name}/
        """
    }

    // get chromosome plot png + pdf file and segment csv file
    process get_chrom_plots_and_segments{

        label 'teeny_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (filename ==~ /google.*/) null
            else if (params.output_to_quartz) "${quartz_path}/facets/chromosome_plots/$filename"
            else "${sample_name}/facets/chromosome_plots/$filename"
        }

        input:
        set sample_name, quartz_path, file("*") from facet_script_output_channel_1
        file get_chrom_plots_and_segments_script from get_chrom_plots_and_segments_file
        file snipdx_gene_coords from snipdx_gene_coords

        output:
        tuple val("${sample_name}"), file("*") into chrom_plots_out_ch
        file(".command.log") into chrom_plots_log_out_ch

        script:
        """
        python get_chrom_plots_and_segments.py \
            "${sample_name}".cnv.table.txt \
            "${sample_name}".jointseg.csv  \
            "${snipdx_gene_coords}"
        """
    }

    // get chromosome plot png + pdf file and segment csv file
    process get_chrom_plots_and_segments_w_purecn_diplogr {

        label 'teeny_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (filename ==~ /google.*/) null
            else if (params.output_to_quartz) "${quartz_path}/facets_with_purecn_diplogr/chromosome_plots/$filename"
            else "${sample_name}/facets_with_purecn_diplogr/chromosome_plots/$filename"
        }

        input:
        set sample_name, quartz_path, file("*") from facet_w_purecn_diplogr_output_channel
        file get_chrom_plots_and_segments_script from get_chrom_plots_and_segments_file
        file snipdx_gene_coords from snipdx_gene_coords

        output:
        tuple val("${sample_name}"), file("*") into chrom_plots_w_purecn_diplogr_out_ch
        file(".command.log") into chrom_plots_w_purecn_diplogr_log_out_ch

        script:
        """
        python get_chrom_plots_and_segments.py \
            "${sample_name}".cnv.table.txt \
            "${sample_name}".jointseg.csv  \
            "${snipdx_gene_coords}"
        """
    }
}


// process for generating an R Markdown Report
// based on the results form the sample facets data
if (generate_report){

    rmarkdown_files_ch.combine(facet_script_output_channel_2, by:0).combine(chrom_plots_out_ch, by:0).set { rmarkdown_report_all_inputs_ch }

    process rmarkdown_report_process{

        label 'reg_task'

        publishDir "${output_dir}/", mode: 'copy',  saveAs: { filename ->
            if (params.output_to_quartz) "${quartz_path}/report/$filename"
            else "${sample_name}/report/$filename"
        }

        input:
        file rmarkdown_report_script from rmarkdown_report_script_ch
        file panel_ga from panel_ga
        set sample_name, trial_name, file(inventory_file), file(vendor_cnv_file), file(vendor_stats_file), file(vendor_variant_file), quartz_path, file('*'), file('*') from rmarkdown_report_all_inputs_ch

        output:
        file("*SNiPDx_by_samples.html") into rmarkdown_report_html_ch
        file("*SNiPDx_by_samples_Inventory.xlsx") into sample_inventory_ch
        file("${sample_name}/*.png") into sample_input_png_out_ch
        file(".command.log") into rmarkdown_report_log_ch

        script:
        """
        ls
        mkdir ${sample_name}
        mv *.png  ${sample_name}/
        mv *.cnv.table.txt  ${sample_name}/
        mv *segments.csv ${sample_name}/
        mkdir VendorDat/
        mv ${vendor_cnv_file} VendorDat/
        mv ${vendor_stats_file} VendorDat/
        mv ${vendor_variant_file} VendorDat/
        echo "rmarkdown::render('${rmarkdown_report_script}', \
            knit_root_dir = '\$PWD/', \
            output_dir = '\$PWD/', \
            output_file = '${sample_name}.${trial_name}_SNiPDx_by_samples.html', \
            params=list( \
                inInventory = '${inventory_file}', \
                inSamplesFolder = '${sample_name}/', \
                inVendorData = 'VendorDat/', \
                currTrial = '${trial_name}', \
                currSample = '${sample_name}', \
                panalAnnotations = '${panel_ga}' \
                ) \
            )" > generate_report.R
        Rscript generate_report.R
        """
    }
}


// collect the command line parameters for historical record keeping
if (!params.output_to_quartz){
    process batch_run_info{

        label 'teeny_task'

        publishDir "${output_dir}/run_info/", mode: 'copy', saveAs: { filename ->
            if (filename ==~ /google.*/) null
            else "$filename"
        }

        input:
        val report_run_cmd_line_info from run_cmd_line_info
        file metadata_file from run_info_metadata_file_ch

        output:
        file("*") into run_info_out_ch

        script:
        """
        TIMESTAMP=`date --rfc-3339=seconds`
        echo "${report_run_cmd_line_info}" > run_info.txt
        echo "Current Date + Time          : \$TIMESTAMP" >> run_info.txt
        if [ -d input_files ]; then
            rm -R input_files
        fi
        mkdir input_files; mv $metadata_file input_files/
        """
    }
}
