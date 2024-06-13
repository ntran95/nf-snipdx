## nf-snipdx: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),

_____________________________________________________________________

## TODO:

- change name to SNiP-Dx
- add production profile that outputs to qwuartz-bio bucket in GCP at the appropriate folder, with pipeline version 
  subfolder
- when doing prod pipelines check that basefolder exists, make sure you aren't overwriting data

_____________________________________________________________________________________

## [2.5.0] - 10-02-2022

### Added

- gcpClin profile for running the pipeline in `repare-clinical`

### Changed

- Moved some configs to seperate sub-configs in the `conf/` directory for reusability

## [2.4.1] - 2022-06-27

### Fixed

- Changed the `trial_name_long` variable assigned in the loop reading the metadata file. This was overwriting the
    global dict value causing errors when there was more than one sample in the metadata file.

## [2.4] - 2022-06-14

### Added

- 2 entry modes now - one for normal snipDX that includes optional reporting, one for just reporting

### Changed

- removed trial profiles, trial information will be hardcoded and paths inferred from it. Reporting will
    require this structure to be accurate to work.
- moved updates to `CHANGELOG.md` file
- changed `gls` profiles to `gcp` since its more understandable
- changed machine configurations on GCP to use preemptibles (to lower cost) and specify machines types

## [2.3]

- added R Markdown report process for automatically generating reports

## [2.2]

- changed name of repo to SNiPDX

## [2.1]

- updated sniprx_tools to include additional python packages, updated pipeline to allow for local insert size 
  normalization around individual snps
- added new PON RData file, and included a reference pon target location file in the ref profile

## [2.0]

- updated codebase to SnipRX (ver 2.0)

## [1.0]

- added batch run mode
- added run_info process to record run parameters, and stderr/stdout logs for each process
- standardized pipeline stdout display

## [0.1] (Beta)

- Converted pipeline from bash scripts to NextFlow pipeline
- Created Docker file that contains all dependencies to run the pipeline
- Created a test data set for installation testing and basic end-to-end validation
