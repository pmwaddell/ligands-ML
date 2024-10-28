import requests
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

    # If max_rec >= 0, there is no maximum:
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
    phos_df.drop("Charge", axis=1, inplace=True)

    # Require MW <= 400:
    phos_df = phos_df[phos_df.MolecularWeight <= 400]

    # Require SMILES does not have "+", "-" or "."
    # This removes cases which include multiple molecules ("disconnected structures") and salts.
    phos_df = phos_df[~phos_df.CanonicalSMILES.str.contains(r"\.")]
    phos_df = phos_df[~phos_df.CanonicalSMILES.str.contains(r"\+")]
    phos_df = phos_df[~phos_df.CanonicalSMILES.str.contains(r"\-")]

    # Require one P atom only and only main group elements:
    # parse molecular formula, etc.

    # could we also enforce no salts and restrict to cases that represent just one molecule?

    return phos_df.reset_index()


def draw_from_phos_df(phos_df: pd.DataFrame, filename: str="phos.png") -> None:
    """Draws the set of phosphines from a phosphine DataFrame."""
    Draw.MolsToGridImage(
        [Chem.MolFromSmiles(smiles) for smiles in phos_df["CanonicalSMILES"]],
        molsPerRow=20,
        subImgSize=(400,400),
        legends=phos_df['CID'].astype(str).to_list()
    ).save("images/" + filename)


# prettyprint_phos_df(
#     request_monophosphines(max_records=100)
# )

draw_from_phos_df(request_monophosphines(max_records=1000))
