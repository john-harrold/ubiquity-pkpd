## Common Changelog

This is a common changelog for the following repositories:

- [ubiquity PKPD tools](https://github.com/john-harrold/ubiquity-pkpd)
- [ubiquity R Package](https://github.com/john-harrold/ubiquity)

## 2023-12-09
- Removed the NONMEM and Monolix outputs
- Added translation of rxode2 output into NONMEM and Monolix outputs 

## 2023-12-09
- Fixed bug in nlmxir2 output where comments in IIV broke things
- Fixed bug in R build where spaces in file path broke compilation

## 2022-03-14
- updated custom ubiquity apps (`ubiquity_server_custom.R`)

## 2021-12-31

### R Workflow
- Removed some cruft (extra files that were no longer used and should have
  been removed before)
- Added documentation for creating user defined observation functions when running parameter estimation
- Updated the reporting template to work with the development version of `onbrand`
- Changed terminal messaging to use `cli`
- Added termination criteria in estimation output
- Added template creation for NONMEM, Monolix, and nlmixr

### ubiquity 
- Fixed wrapping issue in Fortran output where unwrappable strings resulted in infinte loop
- Added the ability to specify `<EST:LT>` parameters to be estimated in the log space
- Added the ability to specify `<EST:P>` parameters to be estimated leaving the others fixed
- Deprecaited the `<GUIDE>` descriptor
- Added nlmxir output target
- Depreciated the CVTYP now using CVINTERP
- Added system-testing.txt example 

## 2021-05-05

### Matlab Workflow
- Fixed bug in `system_plot_cohorts` where titles were not being placed on
  figures

### R Workflow

- Lots of small updates to function documentation and the vignettes

- Completely replaced the reporting functionality. Now using the `onbrand`
  package for templated reporting.  **Note:** This will break all reporting in 
  version 1 of the R package. Old Word and PowerPoint templates will work, but
  you will need to create a yaml mapping file.

- Reporting changes added the following functions
  - `system_rpt_add_slide()`
  - `system_rpt_add_doc_content()`
  - `system_rpt_read_template()`
  - `system_rpt_save_report()`
  - `system_rpt_template_details()`
  - `system_fetch_rpt_officer_object()`
  - `system_fetch_rpt_onbrand_object()`
  - `system_set_rpt_officer_object()`
  - `system_set_rpt_onbrand_object()`

- Reporting changes removed the following functions:
  - `system_report_doc_add_content()`
  - `system_report_doc_format_section()`
  - `system_report_doc_set_ph()`
  - `system_report_estimation()`
  - `system_report_glp()`
  - `system_report_init()`
  - `system_report_nca()`
  - `system_report_ph_content()`
  - `system_report_save()`
  - `system_report_set()`
  - `system_report_slide_content()`
  - `system_report_slide_section()`
  - `system_report_slide_title()`
  - `system_report_slide_two_col()`
  - `system_report_view_layout()`
  - `system_fetch_report()`
  - `system_fetch_report_format()`
  - `md_to_officer()`
  - `md_to_oo()`

- Reporting changes updated vignettes, example scripts and function templates


## 2021-04-25

### ubiquity 
- Added Fortran wrapping at "." character to account for conditional
  statements that are too long
- Removed debug printing of rate information  

### R Workflow
- Fixed error in server.R for Shiny app when the tic() and toc() commands 
  were removed 
- Added tic() and toc() functions

## 2021-04-18

### R Package V 1.0.4 sent to CRAN

### R Workflow
- CRAN erroring out because of tic() function. Removed tic() and toc()
  functions
- Other small updates to get things ready for CRAN submission

## 2021-03-15

### R Workflow
- Fixed a bug where the dataset given to system_run_nca does not have enough valid data to actually run NCA

## 2021-03-07

### R Workflow
- Fixed "Coordinate system already present..."  warning in gg_axis
- Added shaded region for observed AUC in NCA reporting
- Fixed location of table generation in example script analysis_nca_md.R so that table styles would be picked up properly
- Added scales package requirement
- Fixed missing rptname inputs to system_report_ph_content

## 2021-02-15

### R Workflow
- Added cohort-specific output times option to the estimation workflow
- Bug in system view when there are cohort-specific parameters defined
- Fixed figure generation errors in the estimation workflow
- Updated references in the documentation and templates to point to r.ubiquity.tools


## 2020-12-31

### R Workflow
- Updated Reporting vignette
- Updated NCA vignette
- Added checks in estimation routines to check for reasonable bounds for
  global optimizer and to notify users when estimates are near bounds

## 2020-12-27

### R Workflow
- Added sessionInfo() to estimation Word reporting
- Cleaned up some of the reporting code and fixed some random bugs

## 2020-12-26

### R Workflow
- Fixed a bug in system_report_save()

## 2020-12-23

### R Workflow
- Added system_fetch_report_format
- Renamed system_report_fetch to system_fetch_report
- Integrated the markdown (md) header format option with flextable outputs
- Added default for cfg$reporting$enabled (FALSE)
- Added the ability to use markdown for the NCA summary tables generated with system_nca_summary
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
