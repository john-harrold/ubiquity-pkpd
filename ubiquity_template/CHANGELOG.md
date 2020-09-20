## Common Changelog

This is a common changelog for the following repositories:

- [ubiquity PKPD tools](https://github.com/john-harrold/ubiquity-pkpd)
- [ubiquity R Package](https://github.com/john-harrold/ubiquity)


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
