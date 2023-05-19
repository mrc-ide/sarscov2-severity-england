source("global_util.R")

version_check("sircovid", "0.14.11")
version_check("spimalot", "0.8.23")

## Define date at which the data is capped for analysis
date <- "2022-02-24"

## Five epochs after starting with a single strain model (without vaccination)
## * mid August 2020: Alpha appears, expand strains
## * early December 2020: vaccination starts, expand vaccine classes 
##                        but without boosters
## * early March 2021: delta appears, expand strains
## * mid September 2021: booster programme starts, expand vaccine classes
## * early November 2021: omicron appears, rotate strains
epoch_dates <- c("2020-09-17", "2020-12-07", "2021-03-08", "2021-09-14", "2021-11-01")

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(assumptions, deterministic)

## The baselines are always region-specific
regions <- sircovid::regions("england")

baseline <- lapply(regions, create_baseline,
                   date, NULL, # setting restart_date to NULL
                   epoch_dates, pars$info, assumptions)
names(baseline) <- regions

message("Writing parameters_info.csv")
write_csv(pars$info, "parameters_info.csv")
message("Writing parameters_proposal.csv")
write_csv(pars$proposal, "parameters_proposal.csv")
message("Writing parameters_prior.csv")
write_csv(pars$prior, "parameters_prior.csv")

message("Writing parameters_base.rds")
saveRDS(baseline, "parameters_base.rds")

message("Writing parameters_transform.R")
fs::file_copy("R/transform.R",
              "parameters_transform.R", overwrite = TRUE)

message("Printing supplementary figures")
png("fig_sup_vacc_age.png", units = "in", width = 6, height = 6, res = 300)
supl_fig_vac_age()
dev.off()
