from cpfuncs import pidPlot, FOPDTprocess, calc_FOPDT_closed_loop, calc_SOPDT_closed_loop
import streamlit as st

st.set_page_config(
    page_title="Feedback control playground",
    page_icon="🎚️",
    layout="wide",
) 


with st.expander("🗒️ About", expanded=False):
    st.write("""
    This is a web app to observe how different parameters affect a first
    order model with dead time. The model is simulated using a PID controller. 
    """)

st.markdown("# Process & Controller information")
process_order = st.selectbox("Process order", ["FOPDT", "SOPDT"])
if process_order == "FOPDT":
    st.markdown("The first-order plus dead time model is:")
    st.markdown("#### $$\\tau_p \\frac{dy(t)}{dt} + y(t) = K_p·c(t-\\theta_p) + K_d·d(t-\\theta_p)$$")
    st.markdown('**Parameters of the process**  \nProcess gain $K_p$, process time constant $\\tau_p$ and process dead time $\\theta_p$. $K_d$ refers to the gain of the disturbances  \n  \n**Parameters of the controller**  \nController gain $K_c$, integral time constant $\\tau_I$ and derivative time constant $\\tau_D$.')
elif process_order == "SOPDT":
    st.markdown("The second-order plus dead time model is:")
    st.markdown("#### $$\\tau_s^2 \\frac{d^2y(t)}{dt^2} + 2\\zeta\\tau_s\\frac{dy(t)}{dt} + y(t) = K_p·c(t-\\theta_p) + K_d·d(t-\\theta_p)$$")
    st.markdown('**Parameters of the process**  \nProcess gain $K_p$, second order time constant $\\tau_s$, damping ratio $\\zeta$ and process dead time $\\theta_p$. $K_d$ refers to the gain of the disturbances  \n  \n**Parameters of the controller**  \nController gain $K_c$, integral time constant $\\tau_I$ and derivative time constant $\\tau_D$.')

controller_type = st.selectbox("Controller type", ["P", "I", "D", "PI", "PD", "ID", "PID"])
col1,col2,col3,col4, col5, col6 = st.columns(6)
with col1:
    Kp = st.slider("$K_p$", 0.0, 10.0, value = 3.0)
    Kc = st.slider("$K_c$", 0.0, 10.0, value = 0.0)
with col2:
    if process_order == "FOPDT":
        taup = st.slider("$\\tau_p$", 0.1, 10.0, value = 2.0)
        taus = 0
    elif process_order == "SOPDT":
        taus = st.slider("$\\tau_s$", 0.1, 10.0, value = 2.0)
        taup = 0
    tauI = st.slider("$\\tau_I$", 0.0, 10.0)
with col3:
    thetap = st.slider("$\\theta_p$", 0.0, 10.0)
    tauD = st.slider("$\\tau_D$", 0.0, 2.0)
with col4:
    Kd = st.slider("$K_d$", 0.0, 10.0, value = 1.0)
    if process_order == "SOPDT":
        zeta = st.slider("$\\zeta$", 0.0, 10.0, value = 0.5)
    
DerivativeFix = st.checkbox("Derivative kick fix")

st.markdown("# Simulation data")
st.markdown("Simulation time and input change data. In the case of servo problem, 'input' refers to the set point. In the case of the regulator problem, 'input' refers to the disturbance.")
col1,col2,col3, _, _, _ = st.columns(6)
with col1:
    tf = st.slider("Simulation time", 1, 100, value = 30)
with col2:
    Input_start = st.slider("Time of the input change", 0.0, float(tf))
with col3:
    Input_change = st.slider("Change of the input", 0.0, 10.0, value = 1.0)
col1, col2 = st.columns(2)
with col1:
    Input_change_type = st.selectbox("Input change shape", ["step", "ramp"])
with col2:
    type_of_process = st.selectbox("Type of process", ["servo", "regulator"])    


n = max(1000,tf*100) # Simulation time points

if process_order == "FOPDT":
    dictionary_options = {"Kp":Kp,
                        "taup":taup,
                        "thetap":thetap,
                        "Kc":Kc,
                        "tauI":tauI,
                        "tauD":tauD,
                        "tf":tf,
                        "Input_start":Input_start,
                        "Input_change":Input_change,
                        "n":n,
                        "Input_change_type":Input_change_type,
                        "DerivateKickFix":DerivativeFix,
                        "Kd":Kd,
                        "type_of_process":type_of_process,
                        "controller_type":controller_type}

    t, SP, PV, OP, P, I, D, e, DerivateKickFix = calc_FOPDT_closed_loop(**dictionary_options)
    st.pyplot(pidPlot(t, SP, PV, OP, P, I, D, e, DerivateKickFix)[1])

elif process_order == "SOPDT":
    dictionary_options = {"Kp":Kp,
                        "taus":taus,
                        "thetap":thetap,
                        "Kc":Kc,
                        "tauI":tauI,
                        "tauD":tauD,
                        "tf":tf,
                        "Input_start":Input_start,
                        "Input_change":Input_change,
                        "n":n,
                        "Input_change_type":Input_change_type,
                        "DerivateKickFix":DerivativeFix,
                        "Kd":Kd,
                        "zeta":zeta,
                        "type_of_process":type_of_process,
                        "controller_type":controller_type}
    t, SP, PV, OP, P, I, D, e, DerivateKickFix = calc_SOPDT_closed_loop(**dictionary_options)
    st.pyplot(pidPlot(t, SP, PV[:,0], OP, P, I, D, e, DerivateKickFix)[1])