# HYBRID10.2
To run HYBRID10.2:

1. Put the source code and makefile in the same directory

2. "make" will compile the code

3. Edit the file slurm_submit.peta4-skylake, changing, if necessary:
  3.1 Project to be charged.
  3.2 Nodes and ntasks
  3.3 Wallclock time
  3.4 mail_user
  3.5 Path to application executable
  
You should not need to change anything else.

4. Then to run in batch mode: "sbatch slurm_submit.peta4-skylake"
