# sarscov2-severity-england

This is an [orderly](https://www.vaccineimpact.org/orderly/) repository which contains the analysis to our preprint

> Epidemiological drivers of transmissibility and severity of SARS-CoV-2 in England

## Update

This repository was updated on 6 November 2023, in line with a correction submitted to the [published manuscript](https://doi.org/10.1038/s41467-023-39661-5) in Nature Communications.

Changes were made to model parameters of cross-immunity against hospitalisation and death given severe disease with the Omicron variant, following prior infection with a previous variant.

For the hospitalisation parameter, we describe it in supplementary table S11 as a 0.55 probability, as informed by a 2022 study by [Nyberg et al.](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00462-7/fulltext). In our original submission’s code, however, we had erroneously coded it as 0.45 (i.e. 1 – 0.55), so we corrected this.

For the probability of death, we updated this to 0.33 (i.e. 0.18 / 0.55, in the latest version of this code). As described for our model, disease severity progression is a process conditional upon transition to the immediately preceding severity pathway compartment. In the same study by Nyberg et al., authors report an unconditional hazard ratio of 0.18 (95%CI 0.06 – 0.57) given infection with Omicron amongst unvaccinated individuals with a prior history of infection, compared to those without. We believe the updated parameter is thus a more accurate interpretation of the source literature within the context of our model’s dynamics.

All analyses in the related publication have been rerun, and this version of the code includes refitted parameters for all modelled scenarios. We only found slight variations in inferred variant-specific transmissibility and severity, primarily due to stochasticity in the model’s processes. The interpretation of our results remains unchanged and an erratum was submitted for this publication to account for the aforementioned changes.

We would like to express our sincerest thanks to [Dr Lloyd Chapman](https://www.research.lancs.ac.uk/portal/en/people/lloyd-chapman(83f5fa20-c802-49e2-82e7-43e111f41f70).html) from Lancaster University, who advised on the need for this correction. His support for our work has been paramount to ensure its accuracy and reproducibility.


## Running

A sequence of tasks needs to be run with a set of parameters to generate the final results.  This is sketched out in the [`run.R`](run.R) script, though this is provided only as a form of documentation. In practice these were run over several days on a HPC.

* regions: `c("north_west", "north_east_and_yorkshire", "midlands", "east_of_england", "london", "south_west", "south_east")`

1. Run the `severity_parsed_data` task to prepare the raw data for fitting.
2. Run the `severity_parameters` task.
3. Run the `severity_fits` task for each of the seven regions. 
4. Run the `severity_fits_combined` task.
5. Run the `severity_sensitivity_analysis` task.

## Requirements

The core requirement is our [sircovid](https://mrc-ide.github.io/sircovid/) package and its dependencies. Because that package is in constant development you will probably want to pin your versions of the software to the versions we used for preparation:

```r
remotes::install_github(c(
  "mrc-ide/dust@v0.13.2",
  "mrc-ide/mcstate@v0.9.14",
  "mrc-ide/sircovid@v0.14.11",
  "mrc-ide/spimalot@v0.8.24"))
```

However, you can always install the versions that we are using with

```r
drat:::add("ncov-ic")
install.packages(c("sircovid", "spimalot"))
```

You will also need a recent [orderly](https://www.vaccineimpact.org/orderly/) which can be installed with

```r
drat:::add("vimc")
install.packages("orderly")
```

## License

MIT © Imperial College of Science, Technology and Medicine
