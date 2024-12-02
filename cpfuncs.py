import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import odeint

def FOPDTprocess(y,t,u, d,Kp, taup, thetap, Kd, Input_start):
    if t<(thetap+Input_start):
        dydt = 0.0  # time delay
    else:
        dydt = (1.0/taup) * (-y + Kp * u + Kd*d)
    return dydt

def SOPDTprocess(y,t, u, d, Kp, zeta, taus, thetap, Kd, Input_start):
    if t<(thetap+Input_start):
        dydt = 0.0  # time delay
        dy2dt2 = 0.0
    else:
        y, dydt = y
        dy2dt2 = (-y + Kp * u + Kd*d - 2*zeta*taus*dydt)/taus**2
    return [dydt, dy2dt2]

def calc_FOPDT_closed_loop(Kp, taup, thetap, Kd, Kc, tauI, tauD, tf, Input_start, Input_change, n, Input_change_type,
                           DerivateKickFix, input_ramp_time = 5, type_of_process = "servo", controller_type = "PID"):
    size_temp_point = n/tf
    t = np.linspace(0,tf,n) # create time vector
    P= np.zeros(n)          # initialize proportional term
    I = np.zeros(n)         # initialize integral term
    D = np.zeros(n)         # initialize derivative term
    e = np.zeros(n)         # initialize error
    OP = np.zeros(n)        # initialize controller output
    PV = np.zeros(n)        # initialize process variable
    SP = np.zeros(n)        # initialize setpoint
    Input_step = int(Input_start/(tf/(n-1))+1) # Input start
    d = np.zeros(n)         # initialize disturbance
    if type_of_process == "servo":
        if Input_change_type == "step":
            SP[0:Input_step] = 0.0     # define setpoint
            SP[Input_step:n] = Input_change     # step up
        elif Input_change_type == "ramp":
            SP[0:Input_step] = 0.0
            SP[Input_step:Input_step+int(input_ramp_time*size_temp_point)] = np.linspace(0,Input_change,int(input_ramp_time*size_temp_point))
            SP[Input_step+int(input_ramp_time*size_temp_point):n] = Input_change
    elif type_of_process == "regulator":
        if Input_change_type == "step":
            d[0:Input_step] = 0.0     # define disturbance
            d[Input_step:n] = Input_change
        elif Input_change_type == "ramp":
            d[0:Input_step] = 0.0
            d[Input_step:Input_step+int(input_ramp_time*size_temp_point)] = np.linspace(0,Input_change,int(input_ramp_time*size_temp_point))
            d[Input_step+int(input_ramp_time*size_temp_point):n] = Input_change
    
    y0 = 0.0                # initial condition
    # loop through all time steps
    for i in range(1,n):
        # simulate process for one time step
        ts = [t[i-1],t[i]]         # time interval
        y = odeint(FOPDTprocess,y0,ts,args=(OP[i-1],d[i-1], Kp, taup, thetap, Kd, Input_start))  # compute next step
        y0 = y[1]                  # record new initial condition
        # calculate new OP with PID
        PV[i] = y[1]               # record PV
        e[i] = SP[i] - PV[i]       # calculate error = SP - PV
        dt = t[i] - t[i-1]         # calculate time step
        P[i] = Kc * e[i]           # calculate proportional term
        if tauI != 0:
            I[i] = I[i-1] + (Kc/tauI) * e[i] * dt  # calculate integral term
        else:
            I[i] = 0
        if DerivateKickFix:
            D[i] = -Kc * tauD * (PV[i]-PV[i-1])/dt # calculate derivative term
        else:
            D[i] = +Kc * tauD * (e[i]-e[i-1])/dt
        
        if "P" not in controller_type:
            P[i] = 0
        if "I" not in controller_type:
            I[i] = 0
        if "D" not in controller_type:
            D[i] = 0

        OP[i] = P[i] + I[i] + D[i] # calculate new controller output
    return t, SP, PV, OP, P, I, D, e, DerivateKickFix


def calc_SOPDT_closed_loop(Kp, taus, thetap, zeta, Kd, Kc, tauI, tauD, tf, Input_start, Input_change, n, Input_change_type,
                           DerivateKickFix, input_ramp_time = 5, type_of_process = "servo", controller_type = "PID"):

    size_temp_point = n/tf
    t = np.linspace(0,tf,n) # create time vector
    P= np.zeros(n)          # initialize proportional term
    I = np.zeros(n)         # initialize integral term
    D = np.zeros(n)         # initialize derivative term
    e = np.zeros(n)         # initialize error
    OP = np.zeros(n)        # initialize controller output
    PV = np.zeros((n,2))    # initialize process variable. There are two because of derivative, so we will consider only the first column vector
    SP = np.zeros(n)        # initialize setpoint
    Input_step = int(Input_start/(tf/(n-1))+1) # Input start
    d = np.zeros(n)         # initialize disturbance

    
    if type_of_process == "servo":
        if Input_change_type == "step":
            SP[0:Input_step] = 0.0     # define setpoint
            SP[Input_step:n] = Input_change     # step up
        elif Input_change_type == "ramp":
            SP[0:Input_step] = 0.0
            SP[Input_step:Input_step+int(input_ramp_time*size_temp_point)] = np.linspace(0,Input_change,int(input_ramp_time*size_temp_point))
            SP[Input_step+int(input_ramp_time*size_temp_point):n] = Input_change
    elif type_of_process == "regulator":
        if Input_change_type == "step":
            d[0:Input_step] = 0.0     # define disturbance
            d[Input_step:n] = Input_change
        elif Input_change_type == "ramp":
            d[0:Input_step] = 0.0
            d[Input_step:Input_step+int(input_ramp_time*size_temp_point)] = np.linspace(0,Input_change,int(input_ramp_time*size_temp_point))
            d[Input_step+int(input_ramp_time*size_temp_point):n] = Input_change

    y0 = [0.0, 0.0]         # initial condition PV, dPV/dt
    # loop through all time steps
    for i in range(1,n):
        # simulate process for one time step
        ts = [t[i-1],t[i]]         # time interval
        y = odeint(SOPDTprocess,y0,ts,args=(OP[i-1],d[i-1], Kp, zeta, 
                                            taus, thetap, Kd, Input_start))  # compute next step
        y0 = y[1,:]                  # record new initial condition
        # calculate new OP with PID
        PV[i,:] = y[1,:]           # record PV
        e[i] = SP[i] - PV[i,0]     # calculate error = SP - PV
        dt = t[i] - t[i-1]         # calculate time step
        P[i] = Kc * e[i]           # calculate proportional term
        if tauI != 0:
            I[i] = I[i-1] + (Kc/tauI) * e[i] * dt  # calculate integral term
        else:
            I[i] = 0
        if DerivateKickFix:
            D[i] = -Kc * tauD * (PV[i,0]-PV[i-1,0])/dt # calculate derivative term
        else:
            D[i] = +Kc * tauD * (e[i]-e[i-1])/dt
        
        if "P" not in controller_type:
            P[i] = 0
        if "I" not in controller_type:
            I[i] = 0
        if "D" not in controller_type:
            D[i] = 0

        OP[i] = P[i] + I[i] + D[i] # calculate new controller output
    return t, SP, PV, OP, P, I, D, e, DerivateKickFix

def pidPlot(t, SP, PV, OP, P, I, D, e, DerivateKickFix):
    # plot PID response
    plt.figure(1,figsize=(15,7))
    plt.subplot(2,2,1)
    plt.plot(t,SP,'k-',linewidth=2,label='Setpoint (SP)')
    plt.plot(t,PV,'r:',linewidth=2,label='Process Variable (PV)')
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.grid(alpha = 0.2)
    plt.subplot(2,2,2)
    plt.plot(t,P,'g.-',linewidth=2,label=r'Proportional = $K_c \; e(t)$')
    plt.plot(t,I,'b-',linewidth=2,label=r'Integral = $\frac{K_c}{\tau_I} \int_{i=0}^{n_t} e(t) \; dt $')
    plt.plot(t,D,'r--',linewidth=2,label=r'Derivative = $-K_c \tau_D \frac{d(PV)}{dt}$' if DerivateKickFix else r'Derivative = $K_c \tau_D \frac{d(e)}{dt}$')    
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.grid(alpha = 0.2)
    plt.subplot(2,2,3)
    plt.plot(t,e,'m--',linewidth=2,label='Error (e=SP-PV)')
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.grid(alpha = 0.2)
    plt.subplot(2,2,4)
    plt.plot(t,OP,'b--',linewidth=2,label='Controller Output (OP)')
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.grid(alpha = 0.2)
    plt.tight_layout()
    ax = plt.gcf()
    return PV.max(), ax


def get_Laplace_inverse(num, den):
    s, t = sp.symbols("s t")
    exp = sp.parse_expr(num)/ sp.parse_expr(den)
    inverse = sp.inverse_laplace_transform(exp, s, t)
    return sp.lambdify(t,inverse),exp, inverse # We lambdify this for usage

def plot_Laplace_inverse(num, den):
    f, exp,_ = get_Laplace_inverse(num, den)
    t = np.linspace(0.00001,20,1000)
    y = f(t)
    plt.figure(dpi = 300, figsize = (6,2))
    plt.plot(t,y, color = "red")
    plt.axhline(0, color="black", ls="--")
    plt.xlabel("t")
    plt.ylabel("f(t)")
    return plt.gcf()


if __name__=="__main__":
    FOPDT_dictionary_options = {"Kp":4,
                          "taup":1,
                          "thetap":3,
                          "Kc":1,
                          "tauI":1,
                          "tauD":00,
                          "tf":30,
                          "Input_start":1,
                          "Input_change":1,
                          "n":10000,
                          "Input_change_type":"ramp",
                          "DerivateKickFix":True,
                          "type_of_process":"regulator",
                          "Kd":1,
                          "controller_type":"PID"}
    
    SOPDT_dictionary_options = {
        "Kp":4,
        "taus":1,
        "thetap":3,
        "zeta":0.5,
        "Kd":1,
        "Kc":1,
        "tauI":2,
        "tauD":0.50,
        "tf":30,
        "Input_start":1,
        "Input_change":1,
        "n":10000,
        "Input_change_type":"step",
        "DerivateKickFix":True,
        "type_of_process":"servo",
        "controller_type":"PID"
    }
    t, SP, PV, OP, P, I, D, e, DerivateKickFix = calc_FOPDT_closed_loop(**FOPDT_dictionary_options)
    PID = pidPlot(t, SP, PV, OP, P, I, D, e, DerivateKickFix)
    plt.show(PID[1])

    t, SP, PV, OP, P, I, D, e, DerivateKickFix = calc_SOPDT_closed_loop(**SOPDT_dictionary_options)
    PID_SO = pidPlot(t, SP, PV[:,0], OP, P, I, D, e, DerivateKickFix)
    plt.show(PID_SO[1])
    print("Maximum overshoot: ", PID)