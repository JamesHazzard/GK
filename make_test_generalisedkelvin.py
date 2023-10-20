import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

# Set Yamauchi & Takei 2016 parameters
# Gas constant
R=8.3145
# Reference pressure
Pr=1.5e9
# Reference temperature
TKr=1473.
T0=273.
# Grain size and reference grain size
d=1.e-3
dr=1.e-3
# Reference density
rhor=3291.
# Thermal expansivity
alphaT=3.59e-5
# Bulk modulus
bmod=115.2
# Raw frequency
freq=0.01
A_b=0.664
alpha=0.38
tau_p=6.e-5
Teta=0.94
beta=0.
delphi=0.
gamma=5.
lambdaphi=0.
mu0=74.8
dmudT=-0.0131
dmudP=2.09
eta0=10**22.9
E=476e+03
Va=5.02e-06
solgrad=1.65
sol50=1380

def relaxation_YT16_normtimescale(T_h,tau_n):

    if T_h < 0.91:
        A_p = 0.01
    elif T_h < 0.96:
        A_p = 0.01 + (0.4 * (T_h - 0.91))
    elif T_h < 1:
        A_p = 0.03
    else:
        A_p = 0.03 + beta

    if T_h < 0.92:
        sigma_p = 4
    elif T_h < 1:
        sigma_p = 4 + (37.5 * (T_h - 0.92))
    else:
        sigma_p = 7

    dep=150
    Pg=(dep/30.)
    P=Pg*1.e9
    Tsol=sol50+(solgrad*(dep-50.))
    TK=T_h*(Tsol+273.15)
    Tn=T_h
    
    # Initialise parameters for raw Vs
    if Tn < Teta:
        Aeta=1.
    elif Tn >= Teta and Tn<1.:
        Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
    else:
        Aeta=(1./gamma)*np.exp(-delphi)
    # Work out viscosity given A
    eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)
    # Unrelaxed compliance
    Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*(TK-273))))

    X_b = A_b * (tau_n ** alpha)
    X_p = A_p * np.exp(-((np.log(tau_n / tau_p))**2 / (2 * (sigma_p ** 2))))

    X = X_b + X_p
    
    return X, Ju, eta

def relaxation_YT16(T_h,tau):

    if T_h < 0.91:
        A_p = 0.01
    elif T_h < 0.96:
        A_p = 0.01 + (0.4 * (T_h - 0.91))
    elif T_h < 1:
        A_p = 0.03
    else:
        A_p = 0.03 + beta

    if T_h < 0.92:
        sigma_p = 4
    elif T_h < 1:
        sigma_p = 4 + (37.5 * (T_h - 0.92))
    else:
        sigma_p = 7

    dep=150
    Pg=(dep/30.)
    P=Pg*1.e9
    Tsol=sol50+(solgrad*(dep-50.))
    TK=T_h*(Tsol+273.15)
    Tn=T_h
    
    # Initialise parameters for raw Vs
    if Tn < Teta:
        Aeta=1.
    elif Tn >= Teta and Tn<1.:
        Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
    else:
        Aeta=(1./gamma)*np.exp(-delphi)
    # Work out viscosity given A
    eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)
    # Unrelaxed compliance
    Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*(TK-273))))
    tauM=eta*Ju
    tau_n=tau/tauM

    X_b = A_b * (tau_n ** alpha)
    X_p = A_p * np.exp(-((np.log(tau_n / tau_p))**2 / (2 * (sigma_p ** 2))))

    X = X_b + X_p
    
    return X, Ju, eta

def integrand_compliance_YT16_continuous(T_h,tau,t):

    X, Ju, eta = relaxation_YT16(T_h,tau)
    #integrand = X * (1 - np.exp(-t / tau_n)) * (1 / tau_n)
    integrand = X * (1 - np.exp(-t / tau))

    return integrand

def integrate_compliance_YT16_continuous(T_h,t):

    T_YT16 = np.logspace(start=-256, stop=256, num=3001, base=10.0)
    I_YT16 = cumtrapz(integrand_compliance_YT16_continuous(T_h,T_YT16,t),np.log(T_YT16))[-1]

    return I_YT16

def integrate_compliance_YT16_discreteapprox(T_h,t):

    T_YT16 = np.linspace(-12,3,10)
    delta_ln_tau = np.diff(np.log(10**T_YT16))
    X, Ju, eta = relaxation_YT16_normtimescale(T_h,10**T_YT16)
    tau_i = 10**T_YT16
    J_i = X[1:] * delta_ln_tau
    delta_J = X[1:] * delta_ln_tau * (1 - np.exp(-t/((10**T_YT16[1:])*Ju*eta)))

    return tau_i,J_i,np.sum(delta_J)


theta = np.array([0.9,0.95,1.0])
n_theta = len(theta)
n_t = 1000
_,Ju,eta = relaxation_YT16(0.9,1)
y2s = 365.25*24*60*60
t = y2s*np.linspace(0,100,n_t)
J = np.zeros((n_t,n_theta))
J_GK = J.copy()
t_norm = J.copy()

for i in range(n_t):
    for j in range(n_theta):

        _,Ju,eta = relaxation_YT16(theta[j],1)
        t_norm[:,j] = t / (Ju * eta)
        J[i,j] = integrate_compliance_YT16_continuous(theta[j],t[i])
        tau_i,J_GK_i,J_GK[i,j] = integrate_compliance_YT16_discreteapprox(theta[j],t[i])
        
        J[i,j] += Ju + t[i]/eta
        J_GK[i,j] += Ju + t[i]/eta

fig,ax=plt.subplots(1,1,figsize=(5,5))
colors = ["black","red","green","blue"]
t /= y2s
for j in range(n_theta):
    plot_color = colors[j]
    ax.plot(t,J[:,j],color=plot_color,ls='-',label=f"Th={theta[j]:.2f},YT16")
    ax.plot(t,J_GK[:,j],color=plot_color,ls='--',label=f"Th={theta[j]:.2f},N=10")
ax.set_xlabel(r"$t$ [years]")
ax.set_ylabel(r"$J(t)$")
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig("compliance_GK.jpg",dpi=1200)
