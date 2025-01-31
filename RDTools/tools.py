from rdkit import Chem
from rdkit.Chem import Mol, Draw

from langgraph.types import Command
from langchain_core.tools import tool
from langgraph.prebuilt import ToolNode
from langchain_core.tools.base import InjectedToolCallId
from langchain_core.messages import ToolMessage

from typing_extensions import Any, Annotated


@tool
def load_from_smiles_string(
    tool_call_id: Annotated[str, InjectedToolCallId], smiles_string: str
) -> Mol:
    """Load RDKit Mol object from SMILES string"""
    return Command(
        update={
            "molecule": Chem.MolFromSmiles(smiles_string),
            "messages": [
                ToolMessage(
                    f"Successfully converted {smiles_string} to RDKit Mol",
                    tool_call_id=tool_call_id,
                )
            ],
        }
    )

RDTools: list = [load_from_smiles_string]
RDToolNode: ToolNode = ToolNode(RDTools)