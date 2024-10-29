import requests
import re
from time import sleep

from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd


def prettyprint_phos_df(phos_df: pd.DataFrame) -> None:
    """Prints a phosphine DataFrame in a nice and legible way."""
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        print(phos_df.drop("CanonicalSMILES", axis=1))  # note that dropping does not occur in-place


def request_monophosphines(properties: tuple=
                           ("MolecularFormula","MolecularWeight","CanonicalSMILES","Complexity","Charge"),
                           max_records: int=-1,
                           smarts: str="[CX4,c][PX3]([CX4,c])[CX4,c]") -> pd.DataFrame:
    """Makes an API request to PubChem for phosphines (i.e. PR3) and prunes the results (see prun_phos_df)."""
    api = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"

    props = "property/"
    for p in properties:
        props += p + ","
    props += "/"

    # If max_rec < 0, there is no maximum:
    if max_records > 0:
        max_records = f"?MaxRecords={max_records}"
    else:
        max_records = ""

    request = requests.get(api + "fastsubstructure/smarts/" + smarts + "/" + props + "JSON" + max_records).json()
    # Prevent too many requests from taking place in too short a timeframe:
    sleep(1)

    phos_df = pd.DataFrame(request["PropertyTable"]["Properties"])
    phos_df["MolecularWeight"] = phos_df["MolecularWeight"].astype(float)
    phos_df["Complexity"] = phos_df["Complexity"].astype(int)
    phos_df["Charge"] = phos_df["Charge"].astype(int)
    return prune_phos_df(phos_df)


def prune_phos_df(phos_df: pd.DataFrame) -> pd.DataFrame:
    """Removes undesired results from a PubChem API request for general phosphines."""
    # Require Charge = 0:
    phos_df = phos_df[phos_df.Charge == 0]
    # Charge column is no longer needed, so it can be dropped:
    phos_df = phos_df.drop("Charge", axis=1)

    # Require MW <= 400:
    phos_df = phos_df[phos_df.MolecularWeight <= 400]

    # Require SMILES does not have "+", "-" or "."
    # This removes cases which include multiple molecules ("disconnected structures") and salts.
    phos_df = phos_df[~phos_df.CanonicalSMILES.str.contains(r"\.|\+|\-")]

    # Require one P atom only and only certain main group elements:
    phos_df = phos_df[prune_by_elements(phos_df.MolecularFormula)]

    return phos_df.reset_index()


def prune_by_elements(mol_formulas: pd.Series) -> pd.Series:
    result = []
    for mol_formula in mol_formulas:
        result.append(prune_by_elements_helper(mol_formula))
    return pd.Series(result, index=mol_formulas.index)


def prune_by_elements_helper(mol_formula: str) -> bool:
    # Need to add an explicit 1 to the formula after elements that occur only once:
    indices_to_add_1 = []
    for i, c in enumerate(mol_formula):
        if c.isalpha():
            if i == len(mol_formula) - 1 or mol_formula[i + 1].isupper():
                # The index must increase for each 1 that appears earlier, since the str gets longer each time:
                indices_to_add_1.append((i + 1) + len(indices_to_add_1))
    for i in indices_to_add_1:
        mol_formula = mol_formula[:i] + "1" + mol_formula[i:]

    # Split the string into each element and the number of times it appears:
    s = re.split(r"(\d+)", mol_formula)
    if s[-1] == "":
        s.pop()

    elems_to_count = {}
    for i in range(len(s)):
        if i % 2 == 0:
            elems_to_count[s[i]] = int(s[i + 1])

    desired_elems = {
        "H",
        "B", "Ga", "In", "Tl",
        "C", "Si", "Ge", "Sn", "Pb",
        "N", "P", "As", "Sb", "Bi",
        "O", "S", "Se", "Te",
        "F", "Cl", "Br", "I"
    }

    # We only want monophosphines, so limit P count to 1:
    if elems_to_count["P"] > 1:
        return False
    # Limit ourselves to compounds that contain only the desired elements:
    for elem in elems_to_count.keys():
        if elem not in desired_elems:
            return False
    return True


def draw_from_phos_df(phos_df: pd.DataFrame, filename: str="phosphine_set.png", legend: str="CID") -> None:
    """Draws the set of phosphines from a phosphine DataFrame."""
    if legend:
        Draw.MolsToGridImage(
            [Chem.MolFromSmiles(smiles) for smiles in phos_df["CanonicalSMILES"]],
            molsPerRow=20,
            subImgSize=(400,400),
            legends=phos_df["legend"].astype(str).to_list()
        ).save("images/" + filename)
    else:
        Draw.MolsToGridImage(
            [Chem.MolFromSmiles(smiles) for smiles in phos_df["CanonicalSMILES"]],
            molsPerRow=20,
            subImgSize=(400,400)
        ).save("images/" + filename)


if __name__ == "__main__":
    request_monophosphines(max_records=1500).to_csv('data/phosphine_set_2.csv')
