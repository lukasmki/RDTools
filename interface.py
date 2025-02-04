import gradio as gr

from RDTools import RDToolNode, RDTools
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL.PngImagePlugin import PngImageFile

from typing import Annotated
from typing_extensions import TypedDict

from langchain_core.messages import HumanMessage, AIMessage, ToolMessage
from langgraph.prebuilt import create_react_agent
from langgraph.managed import RemainingSteps
from langgraph.graph.message import add_messages
from langchain_ollama import ChatOllama


# Setup agent
class StateSchema(TypedDict):
    messages: Annotated[list, add_messages]
    remaining_steps: RemainingSteps
    molecule: str


graph = create_react_agent(
    ChatOllama(
        model="llama3.1:8b",
        temperature=0.0,
        num_predict=256,
    ),
    RDTools,
    state_schema=StateSchema,
)


# Setup interface
class ChatState:
    history: list = None
    molecule: str = None


def submit_message(chat_state, message):
    user_message = [{"role": "user", "content": message}]

    if chat_state.history is None:
        chat_state.history = user_message
    else:
        chat_state.history += user_message

    return chat_state, chat_state.history, ""


def new_respond(chat_state: ChatState):
    payload = {"messages": chat_state.history}
    if chat_state.molecule is not None:
        print("MOLECULE ADDED TO CONTEXT", chat_state.molecule)
        payload.update({"molecule": chat_state.molecule})

    response_content = ""
    for state in graph.stream(
        payload,
        {"configurable": {"thread_id": "thread-1"}},
        stream_mode="values",
        debug=False,
    ):
        new_message = state["messages"][-1]
        print(new_message)

        if isinstance(new_message, AIMessage):
            response_content += new_message.content

        if "molecule" in state:
            print("MOLECULE ADDED TO STATE", state["molecule"])
            chat_state.molecule = state["molecule"]

        ai_message = {"role": "assistant", "content": response_content}

        yield chat_state, chat_state.history + [ai_message]

    chat_state.history.append(ai_message)


def update_image(chat_state: ChatState) -> PngImageFile:
    if chat_state.molecule is None:
        return None
    else:
        return Draw.MolToImage(Chem.MolFromSmiles(chat_state.molecule))


if __name__ == "__main__":
    with gr.Blocks() as app:
        chat_state = gr.State(ChatState())

        with gr.Row(equal_height=True):
            with gr.Column():
                image = gr.Image(interactive=False)

            with gr.Column():
                chat_box = gr.Chatbot([], type="messages")
                chat_input = gr.Textbox("", show_label=False)

        chat_input.submit(
            submit_message,
            [chat_state, chat_input],
            [chat_state, chat_box, chat_input],
        ).then(
            new_respond,
            [chat_state],
            [chat_state, chat_box],
        ).then(
            update_image,
            chat_state,
            image,
        )
    gr.ChatInterface

    app.launch()
