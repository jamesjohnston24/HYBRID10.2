# HYBRID10.2

To run HYBRID10.2:

1. Put the source code and makefile in the same directory

2. "make" will compile the code

To run in batch mode:

3. Edit the file slurm_submit.peta4-skylake, changing, if necessary:
    Project to be charged.
    Nodes and ntasks
    Wallclock time
    mail_user
    Path to application executable
  Edit 'shared.f90' ntasks to equal the ntasks in the slurm file
Then to run submit: "sbatch slurm_submit.peta4-skylake"

To run interactively:

4. Edit ntasks in 'shared.f90' to = 4.
5. Make sure the run will only last up to about 1 min. by limiting the no. years in 'shared.f90' (I just did a 1 yr run, which took 50 s).
6. "time nice -19 mpirun -n 4 ./dev.exe" will run it. Restart files are produced that can be used for a new run by changing the RSF logical to .TRUE. in 'dev.f90', and making sure the path in 'RSF_In.f90' is the one you want.
