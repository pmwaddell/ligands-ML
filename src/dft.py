# https://sites.google.com/site/orcainputlibrary/geometry-optimizations
# Use GGA DFT functionals if they are accurate enough (depends on your system), with the RI-J approximation (default)
# as that is often the fastest useful optimization one can do. Use of the RI-J approximation leads to minimal
# geometrical errors. Often the slightly higher accuracy from hybrid functionals is not worth the effort.
import os
import subprocess
import datetime
from utils import mkdir


def make_geom_opt_inp_from_xyz(xyz_filename: str, inp_destination_path: str,
                               functional: str="B3LYP", basis_set: str="def2-SVP") -> None:
    """Produces an ORCA input file for geom. opt. from xyz file."""
    header = f"! {functional} {basis_set} Opt\n\n* xyz 0 1\n"
    with open(xyz_filename, 'r') as xyz_file:
        # Remove the initial lines of the xyz file, leaving only the atoms and their coordinates:
        inp_contents = "\n".join(xyz_file.read().splitlines()[2:])

    with open(inp_destination_path, "w") as inp_file:
        inp_file.write(header + inp_contents + "\n*")

# TODO: change log messages so they all start with the CID?
def geom_opt(functional: str="B3LYP", basis_set: str="def2-SVP", redo_all=False) -> None:
    """Performs geometry optimizations based on the xyz files in the conf_search directory."""
    log = f"Geometry optimization started {datetime.datetime.now()}\nCreating directories:\n"
    for filename in os.listdir("data/cf_tst"):
        cid = filename[:-4]
        # Make directories and input files from the xyz files from the conformer search:
        try:
            mkdir(f"data/geom_opt/{cid}")
        except Exception as e:
            log += f"{cid}: Directory creation in geom_opt unsuccessful with exception {e}.\n"
            continue

        make_geom_opt_inp_from_xyz(
            functional=functional,
            basis_set=basis_set,
            xyz_filename=f"data/conf_search/{cid}.xyz",
            inp_destination_path=f"data/geom_opt/{cid}/{cid}.inp")

    log += "\nGeometry optimzation part:\n"
    print()
    for filename in os.listdir("data/cf_tst"):
        cid = filename[:-4]

        if not redo_all:
            # If an existing .out file is found in the directory, check if it seems that it was from a successful calc.
            # If so, skip doing the geometry optimization for that CID.
            if os.path.exists(f"data/geom_opt/{cid}/{cid}.out"):
                with open(f"data/geom_opt/{cid}/{cid}.out", 'r') as out_file:
                    if out_file.read().splitlines()[-2].strip() == "****ORCA TERMINATED NORMALLY****":
                        msg = (f"Existing .out file from completed calculation found in data/geom_opt/{cid}/, "
                               f"skipping calculation.")
                        print(msg)
                        log += msg + "\n"
                        continue
                    else:
                        msg = (f"Existing .out file from apparently unsuccessful calculation found "
                               f"in data/geom_opt/{cid}/, redoing calculation.")
                        print(msg)
                        log += msg + "\n"

        print(f"Performing geometry optimization on {cid}: ", end="")
        orca_command = f"orca data/geom_opt/{cid}/{cid}.inp > data/geom_opt/{cid}/{cid}.out"
        subprocess.run(orca_command, shell=True)
        print("complete.")

    print("\nGeometry optimizations complete.")

    with open("logs/geom_opt_log.txt", "w") as log_file:
        log_file.write(log)


if __name__ == "__main__":
    geom_opt()
