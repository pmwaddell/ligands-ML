from typing import Callable

import pandas as pd


def mod_phos_smiles(phos_smiles, mod: str) -> str:
    """Modifies a phosphine and modifies it according to the input string."""
    for i, c in enumerate(phos_smiles):
        if c == "P":
            return phos_smiles[:i] + mod + phos_smiles[i+1:]


def phos_smiles_to_nico3_complex_smiles(phos_smiles: str) -> str:
    """Takes SMILES string for phosphine and returns SMILES string for corresponding R3P-Ni(CO)3 complex."""
    return mod_phos_smiles(phos_smiles, "[P+]([Ni]([C-]#[O+])([C-]#[O+])([C-]#[O+]))")


def phos_smiles_to_phos_oxide_smiles(phos_smiles: str) -> str:
    """Takes SMILES string for phosphine and returns SMILES string for corresponding phosphine oxide."""
    return mod_phos_smiles(phos_smiles, "P(=O)")


def make_mod_phos_df(phos_df: pd.DataFrame, mod: Callable) -> pd.DataFrame:
    """Returns a DataFrame of modified phosphines from a phosphine dataframe."""
    # TODO: fill out the df with additional properties i.e. MW that we can get from the Mol or whatever?
    result = pd.DataFrame({"CanonicalSMILES": [mod(smiles) for smiles in phos_df.CanonicalSMILES]})
    result["CID"] = phos_df.CID
    return result


# from request_monophosphines import draw_from_phos_df
# df = pd.read_csv('data/phosphine_set.csv')
# comp_df = make_mod_phos_df(df, phos_smiles_to_nico3_complex_smiles)
# draw_from_phos_df(comp_df, filename='P-NiCO3_set.png')
