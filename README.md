# Clonal dynamics code - Lineage Tracing Model (LTM)

## introduction

This codebase is used for simulations of lineage tracing and subclone labelling in the manuscript Bailey_Bhargava_Fu_etal ([Preprint](https://doi.org/10.1101/2025.09.15.675999); Manuscript under review).

## overview of codes 

* *initTumour.hpp* and *initTumour.cpp* include functions for initialising the data structures for the tumour and writing outputs. Key model parameters are stored in *initTumour.hpp*. 
* *evolveTumour.hpp* and *evolveTumour.cpp* include functions for simulating the spatio-temporal evolution of tumour growth and clonal dynamics. 
* *particleCell_main.cpp* include codes to call key functions for running a simulation.
* *Makefile* instructs code complilation. 
* *virtualTumour* is an executible created after the code is successfully compiled. 
* *slurm-submit-particleCell.sh* configures cluster jobs on HPC using the SLURM system. 
* *visualise_subclone_patterns.ipynb* is a Jupyter Notebook providing example codes to visualise the subclone patterns in the simulated tumours. 

## demo simulations 

The demo simulations are designed to show the effect of the cell's intrinsic random motility on subclone patterns. Specifically, these two parameter settings are employed:
- **BROWNIAN_DIFF_SCALE = 0**, which reflects Low Intrinsic Motility. Codes are in the subfolder "demo_subclone_labelling_d5_low_motility"
- **BROWNIAN_DIFF_SCALE = 6e5**, which reflects High Intrinsic Motility. Codes are in the subfolder "demo_subclone_labelling_d5_high_motility"

To make the computational cost manageable (under 1 hour of walltime using 4 CPU cores per simulation), these parameters are selected:
- **T = 240**, which reflects model time of 10 days.
- **NUM_HOUR_INDUCE_SUBCL = 168**, which reflects subclone labelling at day 5. 

**seed = 123** has been set for the random number generator for reproducible results in this demo. 

## how to run the demo 

Note that the following applies to Linux environment on a High-Performance-Computing cluster at the Francis Crick Institute. 
Information about `uname --kernel-name --kernel-release --machine` returns `Linux 4.18.0-553.126.1.el8_10.x86_64 x86_64` (as of August 2026).

1. Compile the code via `make -j` in a command line terminal. If `virtualTumour` already exists before compilation, run `make clean` first. Successful compilation should generate an executible called `virtualTumour`. 

2. Run a simulation. If on a local device, execute `./virtualTumour`; If on an HPC using the SLURM job submission system, execute `sbatch slurm-sumbit-particleCell.sh` (MAKE SURE to modfiy the email address for receiving HPC job information! Also update other job configuration as the user sees fit.)

3. Collect output files. `*nodeDynamics.txt` records the information about individual cells over time, including the lineage ids; `*tumour_size.txt` records the tumour sizes over time.

4. Explore the script to visualise subclone patterns. Follow the script in `visualise_subclone_patterns.ipynb`

## feedback and contact

Any feedback on the agent-based model, codes and demo simulations would be very welcome and should be sent to Xiao Fu `xiao.fu@crick.ac.uk`. For wider discussions on the clonal dynamics project and manuscript, please also contact Erik Sahai `erik.sahai@crick.ac.uk`.
