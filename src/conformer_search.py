import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs, UFFOptimizeMoleculeConfs


def lowest_conf_search_from_smiles(smiles: str, label: str, destination_path: str, force_field: str="MMFF",
                                   num_confs: int=100, use_seed: bool=False, ) -> None:
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

    Chem.MolToXYZFile(mol=p, filename=f"{destination_path}/{label}.xyz", confId=min_energy_index)


def conformer_search(mols_filename: str, smiles_col: str, label_col: str,
                     destination_path: str, force_field: str) -> None:
    """Performs conformer search for phosphines from a DataFrame from a PubChem request."""
    for index, row in pd.read_csv(mols_filename).iterrows():
        print(f"Conformer search for label {row[f"{label_col}"]}: ", end="")
        lowest_conf_search_from_smiles(smiles=row[f"{smiles_col}"], label=row[f"{label_col}"],
                                       destination_path=destination_path, force_field=force_field)
        print("complete.")


if __name__ == "__main__":
    conformer_search(mols_filename="diimine_data/dis.csv",
                     smiles_col="ligand",
                     label_col="ligand_key",
                     destination_path="diimine_data/conf_search_ligand",
                     force_field="MMFF")
    print()
    conformer_search(mols_filename="diimine_data/dis.csv",
                     smiles_col="NiCO2",
                     label_col="ligand_key",
                     destination_path="diimine_data/conf_search_NiCO2",
                     force_field="MMFF")
    print()
    conformer_search(mols_filename="diimine_data/dis.csv",
                     smiles_col="NiACN",
                     label_col="ligand_key",
                     destination_path="diimine_data/conf_search_NiACN",
                     force_field="MMFF")
    print()
    conformer_search(mols_filename="diimine_data/dis.csv",
                     smiles_col="NiE",
                     label_col="ligand_key",
                     destination_path="diimine_data/conf_search_NiE",
                     force_field="MMFF")