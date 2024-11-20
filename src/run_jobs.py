from src.orca_jobs import orca_batch_job


if __name__ == '__main__':
    print("\n\nNi ethylene complex single point calcs: ")
    orca_batch_job(path_to_xyz_files="data/diimine_data/geom_opt_NiE",
                   destination_path="data/diimine_data/single_pt_NiE",
                   functional="B3LYP",
                   basis_set="def2-TZVP",
                   RI='RIJCOSX',
                   dispersion_correction="D3BJ",
                   charge=1,
                   NMR=True,
                   freq=True,
                   job_type="Single Point Calculation")

    print("\n\nNiCO2 complex single point calcs: ")
    orca_batch_job(path_to_xyz_files="data/diimine_data/geom_opt_NiCO2",
                   destination_path="data/diimine_data/single_pt_NiCO2",
                   functional="B3LYP",
                   basis_set="def2-TZVP",
                   RI='RIJCOSX',
                   dispersion_correction="D3BJ",
                   NMR=True,
                   freq=True,
                   job_type="Single Point Calculation")

# TODO: run these on the second laptop? I guess we'll wait for it to finish the current computation
# we'll have to send the geom_opt_ligand_BP86 folder over as well I believe...
# TODO:
#     print("\n\nFree ligand geom opts at higher level of theory or whatever: ")
#     orca_batch_job(path_to_xyz_files="data/diimine_data/geom_opt_ligand_BP86",
#                    destination_path="data/diimine_data/geom_opt_ligand_B3LYP",
#                    functional="B3LYP",
#                    basis_set="def2-TZVP",
#                    RI='RIJCOSX',
#                    dispersion_correction="D3BJ",
#                    job_type="Geometry Optimization")
# then, single point calculations for those too w freq and NMR

# TODO:
#     print("\n\nFree ligand single pt calcs: ")
#     orca_batch_job(path_to_xyz_files="data/diimine_data/geom_opt_ligand_B3LYO",
#                    destination_path="data/diimine_data/single_pt_ligand",
#                    functional="B3LYP",
#                    basis_set="def2-TZVP",
#                    RI='RIJCOSX',
#                    dispersion_correction="D3BJ",
#                    NMR=True,
#                    freq=True,
#                    job_type="Single Point Calculation")
# Also calcs for protonated ligand for proton affinity eventually I guess