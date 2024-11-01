import os

import pandas as pd


def prettyprint_phos_df(phos_df: pd.DataFrame) -> None:
    """Prints a phosphine DataFrame in a nice and legible way."""
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        print(phos_df.drop("CanonicalSMILES", axis=1))  # note that dropping does not occur in-place


def mkdir(dir_name: str) -> None:
    """Creates a directory."""
    try:
        os.mkdir(dir_name)
        print(f"Directory '{dir_name}' created successfully.")
    except FileExistsError:
        print(f"Warning: Directory '{dir_name}' already exists.")
    except PermissionError:
        print(f"Warning: Permission denied: Unable to create '{dir_name}'.")
    except Exception as e:
        print(f"Warning: An error occurred: {e}")