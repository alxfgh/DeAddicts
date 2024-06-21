from rdkit import Chem
from rdkit.Chem import Descriptors

def featurize_smiles(smiles):
    """
    Featurize a SMILES string into a list of molecular descriptors.

    Args:
        smiles (str): A SMILES string.

    Returns:
        list: A list of molecular descriptors.
            - MolWt: Molecular weight.
            - ExactMolWt: Exact molecular weight.
            - MolLogP: LogP.
            - MaxAbsPartialCharge: Maximum absolute partial charge.
            - MinAbsPartialCharge: Minimum absolute partial charge.
            - MaxPartialCharge: Maximum partial charge.
            - MinPartialCharge: Minimum partial charge.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0, 0, 0, 0, 0, 0, 0]
    else:
        try:
            max_abs_partial_charge = Descriptors.MaxAbsPartialCharge(mol)
            min_abs_partial_charge = Descriptors.MinAbsPartialCharge(mol)
            max_partial_charge = Descriptors.MaxPartialCharge(mol)
            min_partial_charge = Descriptors.MinPartialCharge(mol)
        except Exception as e:
            # Handle cases where partial charge calculation fails
            max_abs_partial_charge = 0
            min_abs_partial_charge = 0
            max_partial_charge = 0
            min_partial_charge = 0
        
        return [
            Descriptors.MolWt(mol),
            Descriptors.ExactMolWt(mol),
            Descriptors.MolLogP(mol),
            max_abs_partial_charge,
            min_abs_partial_charge,
            max_partial_charge,
            min_partial_charge,
        ]

# Example usage:
smiles = "CCO"
print(featurize_smiles(smiles))
