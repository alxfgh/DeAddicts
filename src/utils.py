from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def featurize_smiles(smiles):
    """Featurize a SMILES string into a list of molecular descriptors.

    Args:
        smiles (str): A SMILES string.

    Returns:
        list: A list of molecular descriptors.
            - MolWt: Molecular weight.
            - MolLogP: LogP.
            - MaxAbsPartialCharge: Maximum absolute partial charge.
            - MinAbsPartialCharge: Minimum absolute partial charge.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0, 0, 0, 0]
    else:
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.MaxAbsPartialCharge(mol),
            Descriptors.MinAbsPartialCharge(mol),
        ]
