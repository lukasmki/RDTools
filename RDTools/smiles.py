from typing import Self
from rdkit.Chem import Mol, MolFromSmiles, MolFromSmarts, MolToSmiles
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

from rdkit.Chem import DeleteSubstructs, ReplaceSubstructs, CombineMols


class SMILES:
    def __init__(self, smiles: str | Mol):
        self.mol = None
        if isinstance(smiles, Mol):
            self.mol: Mol = smiles
        elif isinstance(smiles, str):
            self.mol: Mol = MolFromSmiles(smiles, sanitize=True, replacements={})
        else:
            print(type(smiles), smiles)

    def __repr__(self):
        return MolToSmiles(self.mol)

    def __str__(self):
        return MolToSmiles(self.mol)

    def remove(self, target_smarts: str) -> Self:
        target = MolFromSmarts(target_smarts)
        new = DeleteSubstructs(self.mol, target)
        return SMILES(new)

    def replace(self, target_smarts: str, fragment_smarts: str) -> Self:
        target = MolFromSmarts(target_smarts)
        fragment = MolFromSmarts(fragment_smarts)
        new = ReplaceSubstructs(self.mol, target, fragment)[0]
        return SMILES(new)

    def add(self, target_smarts: str, fragment_smarts: str) -> Self:
        pass


if __name__ == "__main__":

    print(SMILES("CC(=O)(O)").remove("[OX1]"))
    print(SMILES("CC(=O)(O)").replace("[OX1]", "O"))
