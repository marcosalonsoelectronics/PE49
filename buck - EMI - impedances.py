# -*- coding: utf-8 -*-
from math import pi, log10, sqrt
from control import tf, bode_plot, nyquist_plot
import matplotlib.pyplot as plt

L= 10e-6; rL=0.05; C=33e-6;rC=0.1; R=1; 
Vi= 10; D=0.53;
# Output voltage and current
Vo=D*Vi/(1 + rL/R)
IL=Vo/R
# PI compensator
R1=1.81e3; R2=1.81e3; C2=10e-9
wz=1/(R2*C2)
s = tf('s')
Cp= (1+s/wz)/(s/wz)
# EMI Filter
Lf=100e-6; Cf=10e-6; Rd=100
# Impedances
Zf= 1/( 1/(s*Lf) + s*Cf + 1/Rd)
Zi= (1/D)*( (1+(R+rC)*C*s)*(rL+L*s) + (Vi*Cp+1)*R*(1+rC*C*s)  )/\
    (  D*(1 + (R+rC)*C*s) -IL*Cp*R*(1+rC*C*s) )
# Characteristic function
F= Zf/Zi
# Plot Plant's Bode
# Note that once Hz is true, omega_limits are in Hz
mag, phase, omega = bode_plot(F, dB=True, Hz=True, omega_limits=(10,1e6),\
                              omega_num=200 )
plt.ylim([-540, -45])
#%%
import numpy as np
wmin=10e3; wmax= 1e6
w = np.arange(wmin,wmax,1, dtype=np.float)
nyquist_plot(F, omega=w, color="blue") 



#%%
