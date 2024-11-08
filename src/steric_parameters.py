# # Note: dbstep requires the use of a separate Anaconda environment, at least for me.
# import dbstep.Dbstep as db
# from dbstep.sterics import buried_vol
#
# # make sure that atom1 and atom2 orient things the right way, so we're measuring the phosphine instead of the CO3...
# mol = db.dbstep("phosphine_data/geom_opt_PNiCO3_B3LYP/11776/11776.xyz", atom1=7, atom2=8, sterimol=True, volume=True)
#
# # Sterimol parameters
# print([mol.L, mol.Bmin, mol.Bmax, mol.bur_vol])


# Plan: use the optimized geometries for the NiCO3 complexes, measure along the P-Ni bond:
# we are always guaranteed to have only one P and one Ni so you can just find which atom is which number and use them.

# BUT we may have to clip the COs off, which is doable
# Alternative: I'm pretty sure the O always follows the P, at which point you could use the phosphine oxides as well

# maybe calc both and see which is "better"?





# ALSO use morfeus: https://github.com/digital-chemistry-laboratory/morfeus
from morfeus import read_xyz, ConeAngle, Sterimol


elements, coordinates = read_xyz("data/phosphine_data/geom_opt_PO_BP86/11776/11776.xyz")
# cone_angle = ConeAngle(elements, coordinates, 7)
# print(cone_angle.cone_angle)

sterimol = Sterimol(elements, coordinates, 8, 7)
sterimol.print_report()
