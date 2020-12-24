## Common Changelog

This is a common changelog for the following repositories:

- [ubiquity PKPD tools](https://github.com/john-harrold/ubiquity-pkpd)
- [ubiquity R Package](https://github.com/john-harrold/ubiquity)


## 2020-12-23

### R Workflow

uploaded to github for package:
- Added system_fetch_report_format
- Renamed system_report_fetch to system_fetch_report
- Integrated the markdown (md) header format option with flextable outputs
- Added default for cfg$reporting$ebabled (FALSE)
- Added the support for markdown in flextable headers in NCA summary tables
- Added the ability to use markdown for the NCA summary tables generated with
  system_nca_summary
- Fixed the placeholder text, now the delimiters are === on either side of the text

## 2020-12-20

### R Workflow
- Added as_paragraph output to the code generation in md_to_officer
- Added md_to_oo allow running md_to_officer with small snippets of markdown to be used with flextable

## 2020-12-14

### R Workflow
- Added system_set_option general group and output_directory option
- Updated the default report.docx
- Updated custom shiny app server.R file 

## 2020-12-13

### R Workflow
- For markdown in Word reporting, default font properties have been defined
  and the ability to specify them in the org_functions.R template has been
  added. 
- Moved annotated layout generation of PowerPoint files in
  system_report_view_layout() over to the annotate_base() command in officer
- Updated documentation (functions and vignettes) 
- Added reporting scripts to the unit tests

## 2020-12-05

### R Workflow
- Added import of officer functions starting with body_end_
- Changed how table and figure captions were being numbered in word reporting
- Added optional key fields for tables and figures in word reporting 
- Updated some more elements where list keys are being referenced with $
- Created testthat scripts to run through workshop functions
- Removed aberrant gdata require calls in templates
- Updated documentation

## 2020-11-15

### R Workflow
- Fixed bug in system_nca_run when DOSE is a factor 

## 2020-11-10

### R Workflow
- Fixed bug in system_nca_run where back extrapolation was being done for doses that should have been skipped due to insufficient points. This was causing an error.
- Removed coercion warnings: In tmpsum$halflife = NCA.res$result[NCA.res$result$PPTESTCD == "half.life",  : Coercing LHS to a list
- Added the ability to pass PKNCA.options through to system_nca_run
- Added verbose option to the system_view command
- Allowing pass through of dataset columns to summary NCA output
- Removed digits input from system_nca_run (this is now handled with system_nca_summary below)
- Updated system_view to include nca results
- Updated system_report_ph_content and system_report_doc_add_content to allow for inclusion of flextable objects
- Updated NCA and reporting examples
- Added the following functions: 
  - system_fetch_nca - function to fetch NCA results
  - system_fetch_nca_columns - function to fetch column descriptors of a specific analysis
  - system_nca_parameters_meta - list of standard NCA parameters
  - system_nca_summary - creates summary tables for NCA results
- Updated documentation/vignettes

## 2020-09-27

### R Workflow
- Updated the RData files for the titration vignette 

## 2020-09-27

### R Workflow
- Adding shortcuts to officer exports
- Adding checks in simulate_subjects to ensure required columns are present

## 2020-09-20

### R Workflow
- Fixed coercion warnings for covariates when building the system
- Using explicit declaration of officer functions and specifying those in importFrom this is to prevent namespace issues with readxl to allow the function readxl::read_xlsx
- Added the ability to read xlsx data sets with system_load_data

## 2020-09-12

### R Package V 1.0.3 sent to CRAN

### R Workflow
- Removed URL redirects from documentation to resolve CRAN submission warnings

## 2020-09-11

### R Workflow
- Removed the gdata dependency
- Removed URL redirects from documentation to resolve CRAN submission warnings
- Fixed bug in system_nca_run() that would fail when the dataset was a tibble
- Fixed NCA template
- Converted system_nca_run() from using $ to mostly using [[""]]
- Fixed bug in ShinyApp template where iiv tab was being displayed even when
  the system had no iiv elements


## 2020-07-21
- Updated example system files
- Fixed broken reporting examples from the workshop

### R Workflow
- Updated testthat scripts

## 2020-07-09

### R Workflow
- Bug fixes to stochastic simulations reading from files with covariates

## 2020-07-08

### R Workflow
- Added the option in simulate subjects to specify secondary parameters to be
  saved 

## 2020-07-05

### R Package V 1.0.2 sent to CRAN

### R Workflow
- Uploaded 1.0.2 to CRAN
- Updated the components to fix issues encountered with R 4.0
