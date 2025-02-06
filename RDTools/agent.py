from typing import Annotated, TypedDict, Sequence
from langchain_core.messages import AnyMessage

from langgraph.managed import RemainingSteps
from langgraph.graph import StateGraph, START, END, add_messages

from langgraph.prebuilt import ToolNode


class MolState(TypedDict):
    messages: Annotated[Sequence[AnyMessage], add_messages]
    molecule: str


class MolConfig(TypedDict):
    thread_id: str = "thread-placeholder"


class MolAgent:
    def __init__(self, llm, state_schema, config_schema=None):
        # Create graph
        workflow = StateGraph(
            state_schema=state_schema,
            config_schema=config_schema,
        )

        workflow.add_node("route", self.route)
        workflow.add_node("tools", self.tools)
        workflow.add_node("chat", self.chat)

    # Nodes
    def route(self, state: MolState):
        pass

    def tools(self, state: MolState):
        pass

    def chat(self, state: MolState):
        pass
