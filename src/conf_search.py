import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs, UFFOptimizeMoleculeConfs


def lowest_conf_search_from_smiles(smiles: str, CID: str, destination_directory, num_confs: int=100,
                                   use_seed: bool=False, force_field: str="MMFF") -> None:
    """Performs a search for the lowest energy conformer of the mol from input SMILES string, saves as xyz file."""
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

    Chem.MolToXYZFile(mol=p, filename=f"data/{destination_directory}/{CID}.xyz", confId=min_energy_index)


def conf_search_xyzs_from_phos_set(phos_set_filename: str="data/phosphine_set_redux.csv",
                                   destination_directory="conf_search") -> None:
    """Performs conformer search for phosphines from a DataFrame from a PubChem request."""
    phos_set = pd.read_csv(phos_set_filename)
    phos_set = phos_set.iloc[:, 2:]

    for index, row in phos_set.iterrows():
        print(f"Conformer search for CID {row["CID"]}: ", end="")
        lowest_conf_search_from_smiles(smiles=row["CanonicalSMILES"], CID=row["CID"],
                                       destination_directory=destination_directory)
        print("complete.")


if __name__ == "__main__":
    conf_search_xyzs_from_phos_set(phos_set_filename="data/phosphine_set_redux.csv",
                                   destination_directory="conf_search_P=O")
