# simulation to quantify the prevalence of anomalous 4-taxon networks

## how to run the pipeline

open file #0: edit as desired to set parameters.  
run file #0 from command line, like this to run all first 3 steps:

```shell
Rscript 0_control_params.R
```

or like this to run step 1 only (network simulation with SiPhyNetwork).
Replace 1 by 2 or by 3 to run step 2 only (network classification and
gene tree simulation with PhyloCoalSimulations),
or step 3 only (summarize all replicates from each scenario):

```shell
Rscript 0_control_params.R 1
```

check the number of networks simulated after step 1:
```shell
ls job[0-9]*/SiPhyNetwork_output/sim*.tree | wc -l # overall: should be #jobs * #nets
for jobdir in job[0-9]* # if more detail needed
do
  echo -n "$jobdir: "
  ls $jobdir/SiPhyNetwork_output/sim*.tree | wc -l
done
```

Step 4 files (starting with `4_*`) are for sanity checks, and
to compare the age of the simulated networks with the genetic distance
threshold above which reticulation fails.

Step 5, to plot results (after sanity checks) and make figures using R:
see `5_plot_results.Rmd`.

Non-numbered `.jl` julia files: helper functions used in steps 2 and 4.

## dependencies

[SiPhyNetwork](https://github.com/jjustison/SiPhyNetwork): now on CRAN,
and other packages.
Install within R with:
```r
install.packages("SiPhyNetwork")
install.packages("tictoc")
```

[PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl) and
[PhyloCoalSimulations](https://github.com/cecileane/PhyloCoalSimulations.jl):
and others: all registered julia packages. They are listed, with specific
version numbers, in `Manifest.toml`. Install them within Julia:
first launch julia with `julia --project` to use `Project.toml` then
install all packages with specific versions with this in package mode:
```julia
pkg> instantiate
```
