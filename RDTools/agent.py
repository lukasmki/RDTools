""" Old basic ah agent
def chat(state: State):
    return {"messages": [llm.invoke(state["messages"])]}
def should_continue(state: State):
    print(state)
    last_message = state["messages"][-1]
    if last_message.tool_calls:
        return "tools"
    return END
graph_builder = StateGraph(State)
graph_builder.add_node("chat", chat)
graph_builder.add_node("tools", RDToolNode)
graph_builder.add_edge(START, "chat")
graph_builder.add_conditional_edges("chat", should_continue, ["tools", END])
graph_builder.add_edge("tools", "chat")
graph_builder.add_edge("chat", END)
graph = graph_builder.compile()
"""

