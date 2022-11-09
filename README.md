# six-fingered-hand

2022-06-13: moving files from john's home machine to this github repository.  (they still only work on john's home machine.)

2022-08-30: works on franklin00, running many scenarios in parallel.

pipeline: put everything in one folder.  open file #0.
set parameters as you wish.
run file #0 from command line, like this to run all steps:

```shell
Rscript 0_control_params.R
```

or like this to run step 1 only (network simulation with SiPhyNetwork)
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
