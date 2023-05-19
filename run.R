deterministic <- TRUE
short_run <- TRUE

# assumptions can also be any of:
# crim_infect_high, crim_infect_low
# crim_hospi_high, crim_hospi_low
# crim_death_high, crim_death_low
# booster_ve_high, booster_ve_low
# alpha_ve_high, alpha_ve_low
# delta_ve_high, delta_ve_low
# mu_d_winter, mu_d_summer
# fixed_si_high, fixed_si_low
assumptions <- "central"

## 1. severity_parsed_data
orderly::orderly_run("severity_parsed_data",
                     use_draft = "newer")

## 2. severity_parameters 
orderly::orderly_run("severity_parameters", parameters = list(
  deterministic = deterministic, assumptions = assumptions),
  use_draft = "newer")

## 3. severity_fits
for (r in regions) {
  orderly::orderly_run("severity_fits",
                       parameters = list(region = r,
                                         short_run = short_run,
                                         deterministic = deterministic,
                                         assumptions = assumptions),
                       use_draft = "newer")
}


## 4. severity_fits_combined
orderly::orderly_run("severity_fits_combined",
                     parameters = list(short_run = short_run,
                                       deterministic = deterministic,
                                       assumptions = assumptions),
                     use_draft = "newer")

## 5. severity_sensitivity_analysis
orderly::orderly_run("severity_sensitivity_analysis",
                     parameters = list(short_run = short_run,
                                       deterministic = deterministic),
                     use_draft = "newer")
