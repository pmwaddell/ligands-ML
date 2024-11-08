import pandas as pd
from rdkit import Chem


def format_R_smiles(R):
    if R == "H":
        return ""
    else:
        return f"({R})"


def smiles_to_canonical_smiles(smiles):
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles)
    )


def make_diimine_complex_smiles(R1="H", R2="H", backbone_sub="H", Ni_L1="C", Ni_L2="C"):
    R1, R2, backbone_sub = format_R_smiles(R1), format_R_smiles(R2), format_R_smiles(backbone_sub)
    Ar1 = f"c2c{R1}cc{R2}cc{R1}2"
    Ar2 = f"c3c{R1}cc{R2}cc{R1}3"
    if Ni_L1 == "E" and Ni_L2 == "E":
        return f"C{backbone_sub}1=N({Ar1})[Ni]45(N({Ar2})=C{backbone_sub}1)(CC4)CC5"
    elif Ni_L1 == "E":
        return f"C{backbone_sub}1=N({Ar1})[Ni]4(N({Ar2})=C{backbone_sub}1)(CC4){Ni_L2}"
    elif Ni_L2 == "E":
        return f"C{backbone_sub}1=N({Ar1})[Ni]4(N({Ar2})=C{backbone_sub}1)({Ni_L1})CC4"
    return f"C{backbone_sub}1=N({Ar1})[Ni](N({Ar2})=C{backbone_sub}1)({Ni_L1}){Ni_L2}"


def make_diimine_ligand_smiles(R1="H", R2="H", backbone_sub="H"):
    R1, R2, backbone_sub = format_R_smiles(R1), format_R_smiles(R2), format_R_smiles(backbone_sub)
    Ar1 = f"c1c{R1}cc{R2}cc{R1}1"
    Ar2 = f"c2c{R1}cc{R2}cc{R1}2"
    return f"N({Ar1})=C{backbone_sub}C{backbone_sub}=N({Ar2})"


if __name__ == "__main__":
    subs = [
        "H",
        "C",
        "F",
        "Br",
        "C#N",
        "OC",
        "C(F)(F)F",
        "N(C)C",
        "C(C)C",
        "N(=O)=O",
        "S(C)(=O)=O",
        "C=O",
        "OC(C)C",
        "N(CC)CC"
    ]

    ligand_key, R2, backbone_sub, ligand, NiCO2, NiACN, NiE = [], [], [], [], [], [], []
    for sub in subs:
        R2.append("H")
        backbone_sub.append("H")
        ligand.append(smiles_to_canonical_smiles(make_diimine_ligand_smiles(R1=sub)))
        NiCO2.append(
            smiles_to_canonical_smiles(make_diimine_complex_smiles(R1=sub, Ni_L1="[C-]#[O+]", Ni_L2="[C-]#[O+]")))
        NiACN.append(smiles_to_canonical_smiles(make_diimine_complex_smiles(R1=sub, Ni_L1="C", Ni_L2="[N+]#CC")))
        NiE.append(smiles_to_canonical_smiles(make_diimine_complex_smiles(R1=sub, Ni_L1="C", Ni_L2="E")))

    dis = pd.DataFrame.from_dict(
        {"R1": subs, "R2": R2, "backbone_sub": backbone_sub, "ligand": ligand,
         "NiCO2": NiCO2, "NiACN": NiACN, "NiE": NiE})

    # Make a primary key based on R1, R2, and backbone_sub strings (this will be used to distinguish each ligand):
    ligand_keys = []
    for row in dis.iterrows():
        ligand_keys.append(
            row[1].loc["R1"] + "_" + row[1].loc["R2"] + "_" + row[1].loc["backbone_sub"]
        )

    dis.insert(0, "ligand_key", ligand_keys)
    dis.set_index("ligand_key", inplace=True)
    dis.to_csv("diimine_data/diimine_smiles.csv")
