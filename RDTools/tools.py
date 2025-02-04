from rdkit import Chem
from rdkit.Chem import AllChem, Mol

from langgraph.types import Command
from langchain_core.tools import tool
from langgraph.prebuilt import ToolNode
from langchain_core.tools.base import InjectedToolCallId
from langgraph.prebuilt import InjectedState
from langchain_core.messages import ToolMessage

from typing_extensions import Annotated


@tool
def load_molecule_from_smiles(
    tool_call_id: Annotated[str, InjectedToolCallId],
    molecule: Annotated[str, "Molecule in SMILES format"],
) -> Command:
    """
    Loads SMILES string for manipulation
    """
    return Command(
        update={
            "molecule": molecule,
            "messages": [
                ToolMessage(
                    f"Loaded molecule {molecule}",
                    tool_call_id=tool_call_id,
                )
            ],
        }
    )


@tool
def add_smarts_to_smiles(
    tool_call_id: Annotated[str, InjectedToolCallId],
    state: Annotated[dict, InjectedState],
    attachment_target_smarts: Annotated[
        str, "Fragment to be attached to in SMARTS format"
    ],
    attachment_smarts: Annotated[str, "Fragment to attach in SMARTS format"],
) -> Command:
    """
    Adds a fragment to a target molecule at a specified attachment point.
    """
    input_smiles = state["molecule"]

    input_mol = Chem.MolFromSmiles(input_smiles)
    fragment_mol = Chem.MolFromSmiles(attachment_target_smarts)
    attachment = Chem.MolFromSmarts(attachment_smarts)

    if not input_mol or not fragment_mol or not attachment:
        return Command(update={"messages": [ToolMessage("Input parsing failed")]})

    # Find attachment point
    match = input_mol.GetSubstructMatch(attachment)
    if not match:
        return Command(
            update={
                "messages": [
                    ToolMessage(
                        f"No valid attachment site found for target {attachment_target_smarts}"
                    )
                ]
            }
        )

    # Merge molecules
    combined = Chem.CombineMols(input_mol, fragment_mol)
    combined_rw = Chem.RWMol(combined)

    # Add bond between attachment point and fragment
    target_idx = match[0]
    fragment_idx = input_mol.GetNumAtoms()  # First atom in the fragment
    combined_rw.AddBond(target_idx, fragment_idx, Chem.BondType.SINGLE)

    return Command(
        update={
            "molecule": Chem.MolToSmiles(combined_rw.GetMol()),
            "messages": [
                ToolMessage(
                    f"Added {attachment_smarts} to {input_smiles}",
                    tool_call_id=tool_call_id,
                )
            ],
        }
    )


@tool
def remove_smarts_from_smiles(
    tool_call_id: Annotated[str, InjectedToolCallId],
    state: Annotated[dict, InjectedState],
    removal_target_smarts: Annotated[str, "Fragment to remove in SMARTS format"],
) -> Command:
    """
    Removes a substructure defined by a SMARTS pattern.
    """
    input_smiles = state["molecule"]

    mol = Chem.MolFromSmiles(input_smiles)
    pattern = Chem.MolFromSmarts(removal_target_smarts)

    if not mol or not pattern:
        return Command(update={"messages": [ToolMessage("Input parsing failed")]})

    # Remove substructure
    modified_mol = AllChem.DeleteSubstructs(mol, pattern)
    return Command(
        update={
            "molecule": Chem.MolToSmiles(modified_mol)
            if modified_mol
            else input_smiles,
            "messages": [
                ToolMessage(
                    f"Removed {removal_target_smarts} from {input_smiles}",
                    tool_call_id=tool_call_id,
                )
            ],
        }
    )


@tool
def replace_smarts_in_smiles(
    tool_call_id: Annotated[str, InjectedToolCallId],
    state: Annotated[dict, InjectedState],
    replacement_target_smarts: Annotated[
        str, "Fragment to be replaced in SMARTS format"
    ],
    replacement_smarts: Annotated[str, "Fragment that will replace in SMARTS format"],
) -> Command:
    """
    Replaces a substructure defined by a SMARTS pattern with another fragment.
    """
    input_smiles = state["molecule"]

    mol = Chem.MolFromSmiles(input_smiles)
    pattern = Chem.MolFromSmarts(replacement_target_smarts)
    replacement = Chem.MolFromSmiles(replacement_smarts)

    if not mol or not pattern or not replacement:
        return Command(update={"messages": [ToolMessage("Input parsing failed")]})

    rxn = AllChem.ReactionFromSmarts(
        f"[{replacement_target_smarts}:1]>>[{replacement_smarts}:1]"
    )
    products = rxn.RunReactants((mol,))

    return Command(
        update={
            "molecule": Chem.MolToSmiles(products[0][0]) if products else input_smiles,
            "messages": [
                ToolMessage(
                    f"Replaced {replacement_target_smarts} with {replacement_smarts} in {input_smiles}",
                    tool_call_id=tool_call_id,
                )
            ],
        }
    )


RDTools: list = [
    load_molecule_from_smiles,
    add_smarts_to_smiles,
    remove_smarts_from_smiles,
    replace_smarts_in_smiles,
]
RDToolNode: ToolNode = ToolNode(RDTools)
