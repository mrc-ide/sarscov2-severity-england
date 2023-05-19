# sarscov2-severity-england

This is an [orderly](https://www.vaccineimpact.org/orderly/) repository which contains the analysis to our preprint

> Epidemiological drivers of transmissibility and severity of SARS-CoV-2 in England

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
  "mrc-ide/sircovid@0.14.11",
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
