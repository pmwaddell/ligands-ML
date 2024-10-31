from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs, UFFOptimizeMoleculeConfs


def lowest_conf_search_from_smiles(smiles: str, num_confs: int=100, use_seed: bool=False, force_field: str="MMFF") \
        -> Chem.Mol:
    """Performs a search for the lowest energy conformer of the molecule corresponding to the input SMILES string."""
    p = Chem.MolFromSmiles(smiles)
    p = Chem.AddHs(p)  # Necessary for reasonable conformers.
    if use_seed:
        rdDistGeom.EmbedMultipleConfs(p, numConfs=num_confs, randomSeed=0xf00d)
    else:
        rdDistGeom.EmbedMultipleConfs(p, numConfs=num_confs)

    if force_field == "MMFF":
        results = MMFFOptimizeMoleculeConfs(p)
    elif force_field == "UFF":
        results = UFFOptimizeMoleculeConfs(p)
    else:
        raise Exception("Unrecognized force_field.")

    # Find the lowest energy conformer:
    min_energy = results[0][1]
    min_energy_index = 0
    for i, conf in enumerate(results):
        if results[i][0] == 0 and results[i][1] < min_energy:
            min_energy_index = i
            min_energy = results[i][1]

    return Chem.Mol(p, confId=min_energy_index)


# r = lowest_conf_search_from_smiles("C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3")
# print(r.GetConformers())
