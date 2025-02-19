from typing import Self
from rdkit.Chem import Mol, MolFromSmiles, MolFromSmarts, MolToSmiles
from rdkit.Chem import DeleteSubstructs, ReplaceSubstructs
from rdkit.Chem.Descriptors import CalcMolDescriptors
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts


class SMILES:
    def __init__(self, smiles: str | Mol | list[str | Mol]):
        if not isinstance(smiles, list):
            smiles = smiles.split(".")
        self.mol = [self._convert_to_mol(s) for s in smiles]

    def _convert_to_mol(self, smiles):
        if isinstance(smiles, Mol):
            return smiles
        elif isinstance(smiles, str):
            return MolFromSmiles(smiles, sanitize=True, replacements={})
        else:
            raise ValueError(f"Input type {type(smiles)} not Mol or str {smiles}")

    def __repr__(self):
        return ".".join([MolToSmiles(m) for m in self.mol])

    def __str__(self):
        return ".".join([MolToSmiles(m) for m in self.mol])

    @property
    def description(self) -> list[dict]:
        return [CalcMolDescriptors(mol) for mol in self.mol]

    def remove(self, target_smarts: str) -> Self:
        target = MolFromSmarts(target_smarts)
        self.mol = [DeleteSubstructs(mol, target) for mol in self.mol]
        return self

    def replace(self, target_smarts: str, fragment_smarts: str) -> Self:
        target = MolFromSmarts(target_smarts)
        fragment = MolFromSmarts(fragment_smarts)
        self.mol = [ReplaceSubstructs(mol, target, fragment)[0] for mol in self.mol]
        return self

    def add(self, smiles: str | Mol | list[str | Mol]) -> Self:
        if not isinstance(smiles, list):
            smiles = smiles.split(".")
        self.mol += [self._convert_to_mol(s) for s in smiles]
        return self

    def react(self, reaction_smarts: str) -> Self:
        reaction: ChemicalReaction = ReactionFromSmarts(reaction_smarts)
        products = reaction.RunReactants(self.mol)
        if not products:
            return self
        self.mol = list({MolToSmiles(m) for p in products for m in p})
        self.mol = [self._convert_to_mol(m) for m in self.mol]
        return self


if __name__ == "__main__":
    print(SMILES("CC(=O)(O)").remove("[OX1]"))
    print(SMILES("CC(=O)(O)").replace("[OX1]", "[O]"))
    print(SMILES("CC(=O)(O)").add("CNC"))
    print(
        SMILES("CC(=O)O").add("CNC").react("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")
    )
    print(
        SMILES("C(=O)OC(=O)O")
        .add("CNC")
        .react("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")
    )

    a = SMILES("C(=O)(O)C(=O)(O)")
    print(a)
    a.add("CNC")
    print(a)
    a.react("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")
    print(a)

    print(SMILES("CC(=O)(O)").description)
