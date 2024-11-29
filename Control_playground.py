import streamlit as st

st.set_page_config(
    page_title="Control playground home",
    page_icon="👋",
)

st.write("# Welcome to the Control playground")

st.sidebar.success("Select a tool above.")

st.markdown(
    """
    This is a web app to observe different concepts of Control theory at a basic level.

    - **Laplace transform viewer**: Observe the inverse time response of a Laplace function.

    - **Feedback control playground**: Observe how different parameters affect a first order model with dead time.

"""
)