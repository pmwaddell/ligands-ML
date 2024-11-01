# https://sites.google.com/site/orcainputlibrary/geometry-optimizations
# Use GGA DFT functionals if they are accurate enough (depends on your system), with the RI-J approximation (default)
# as that is often the fastest useful optimization one can do. Use of the RI-J approximation leads to minimal
# geometrical errors. Often the slightly higher accuracy from hybrid functionals is not worth the effort.
import subprocess
import time


start = time.time()

# Construct the ORCA command
orca_command = "orca pph3_manual3.inp > ph3_manual3.out"

# Execute the command
subprocess.run(orca_command, shell=True)

end = time.time()
print(end-start)
# ok so you're going to have to get ORCA working for yourself in order to do all this.


# make directory based on CID
# get xyz file for lowest E conformer
# use that to make input file for ORCA
# then run calculation
# use our script to pull the results out into a DataFrame I guess