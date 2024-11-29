from cpfuncs import get_Laplace_inverse, plot_Laplace_inverse
import sympy as sp
import streamlit as st

st.set_page_config(
    page_title="Laplace transform viewer",
    page_icon="🎚️",
    layout="wide",
) 


with st.expander("🗒️ About", expanded=False):
    st.write("""
    This is a web app to observe the time response of a Laplace transform function. 
    """)

st.markdown("# Laplace domain function")
num = st.text_input("Numerator", value="1")
den = st.text_input("Denominator", value="s+1")

st.markdown("The laplace domain function and its inverse are shown below, where $\\theta(t)$ represents the step function, with a value of 1 at any time $T\\geq t$ and 0 at any time $T<t$ is:")
LI, La_exp, time_exp = get_Laplace_inverse(num, den)
st.markdown(f"### $$F(s) = {sp.latex(La_exp)} \qquad \qquad f(t) =  {sp.latex(time_exp)}$$")
st.pyplot(plot_Laplace_inverse(num, den))