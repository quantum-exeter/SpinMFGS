# SpinMFGS
Spin Mean Force Gibbs State

## Repo structure

* **SMFGS**: contains the code for computing the MFGS expectation values of a spin; structured as a Julia module
* **run**: contains example code using the module; in particular contains the scripts used to generate the plots for the paper
* **data**: contains data generated using this library; in particular contains data used for the plots of the paper
* **plot**: contains a Jupyter notebook with code to plot the results; in particular contains the code used for the plots of the paper

## Usage of the **run** scripts
To use the contained run scripts, the user should change into the folder where the scripts are located and then run them as
```
julia -t n ./script-name.jl
```
where *script-name* is the name of the script that one wants to run, and *n* is the number of threads for parallelization of the computations (for example "-t 3" would spawn 3 threads).
If run this way, the output of the scripts will be saved into the **data** folder. If one wants to change the output directory or filename, all that is required is modifying the *filename* variable defined inside of each script.
