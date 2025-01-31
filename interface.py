import gradio as gr

from RDTools import RDToolNode, RDTools
from rdkit import Chem
from rdkit.Chem import Draw

from langchain_community.tools import WikipediaQueryRun
from langchain_community.utilities import WikipediaAPIWrapper

from typing import Annotated
from typing_extensions import TypedDict

from langchain_core.messages import HumanMessage, AIMessage, ToolMessage
from langgraph.prebuilt import create_react_agent
from langgraph.managed import RemainingSteps
from langgraph.graph.message import add_messages
from langchain_ollama import ChatOllama


def respond(message, history):
    img = None
    content = ""

    user_message = [{"role": "user", "content": message}]
    for state in graph.stream(
        {"messages": history + user_message},
        {"configurable": {"thread_id": "thread-1"}},
        stream_mode="values",
    ):
        new_message = state["messages"][-1]
        if isinstance(new_message, AIMessage):
            content += new_message.content

        if "molecule" in state:
            try:
                img = Draw.MolToImage(state["molecule"])
            except:
                pass

        yield {"role": "assistant", "content": content}, img


class State(TypedDict):
    messages: Annotated[list, add_messages]
    remaining_steps: RemainingSteps
    query: str
    molecule: Chem.Mol


if __name__ == "__main__":
    llm = ChatOllama(
        model="llama3.1:8b",
        temperature=0.0,
        num_predict=256,
    )
    graph = create_react_agent(llm, RDTools, state_schema=State)

    with gr.Blocks() as app:
        with gr.Row(equal_height=True):
            with gr.Column():
                image = gr.Image(interactive=False)

            with gr.Column():
                chat = gr.ChatInterface(
                    respond,
                    additional_outputs=[image],
                    type="messages",
                )

    app.launch()
