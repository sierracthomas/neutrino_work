import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import math



# Declare symbols 
theta_12, theta_13, theta_23, delta, phi2, phi3 = sp.symbols('theta_12 theta_13 theta_23 delta phi2 phi3')
N_e, N_m, N_t = sp.symbols('N_e N_mu N_tau')
m11, m21, m31 = sp.symbols('m11 m21 m31')

nu_1m, nu_2m, nu_3m = sp.symbols('nu_1m nu_2m nu_3m')

nu_e, nu_mu, nu_tau = sp.symbols('nu_e nu_mu nu_tau')

theta_m, theta = sp.symbols('theta_m theta')

# Components for V_MNS
A = sp.Matrix([[1, 0, 0],
              [0, sp.exp(sp.I * phi2 / 2), 0],
              [0, 0, sp.exp(sp.I * phi3 / 2)]])

B = sp.Matrix([[sp.cos(theta_12), sp.sin(theta_12), 0],
              [-sp.sin(theta_12), sp.cos(theta_12), 0],
              [0, 0, 1]])

C = sp.Matrix([[sp.cos(theta_13), 0, sp.sin(theta_13)* sp.exp(-sp.I * delta)],
              [0, 1, 0],
              [-sp.sin(theta_13) * sp.exp(-sp.I * delta), 0, sp.cos(theta_13)]])


D = sp.Matrix([[1, 0, 0],
              [0, sp.cos(theta_23), sp.sin(theta_23)],
              [0, -sp.sin(theta_23), sp.cos(theta_23)]])


# Components from 3.19

eig_m = sp.Matrix([nu_1m, nu_2m])

eig_f = sp.Matrix([nu_e, nu_mu])

mix_2d = sp.Matrix([[sp.cos(theta), -sp.sin(theta)],
                    [sp.sin(theta), sp.cos(theta)]])


# Equation 3.22
def eff_osc_theta(theta = 45, E_nu = 6, R_sun = False, del_m = 8e-5):
    # eq 3.22
    #A_over_delm = np.linspace(0.001, 100, 10000)
    if R_sun:
        sun_data = np.genfromtxt("electron_density.csv", delimiter = ",")
        N_e = sun_data.T[1][::-1]*10e18
        #N_e = N_e[::-1]
        rad = sun_data.T[0]
    else:
        N_e = np.linspace(0.000001, 400, 10000)
    A_over_delm = (1.54e-4 * E_nu * N_e)/del_m
    osc_prob_pos = 0.5 * np.arcsin(np.sqrt(((np.sin(2 * theta))**2) / ((A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2)))
    osc_prob_neg = 0.5 * np.arcsin(np.sqrt(((np.sin(2 * theta))**2) / ((-A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2)))
    
    if R_sun:
        return rad, osc_prob_pos, osc_prob_neg
    else:
        return N_e, osc_prob_pos, osc_prob_neg
    
# Equation 3.21
def eff_osc_m_m(theta = 45, del_m12 = 8e-5, E_nu = 2, R_sun = False):
    
    if R_sun:
        sun_data = np.genfromtxt("electron_density.csv", delimiter = ",")
        N_e = sun_data.T[1][::-1]*10e18
        rad = sun_data.T[0]
    else:
        N_e = np.linspace(0.000001, 400, 10000)
    
    A_over_delm = (N_e * 1.54e-4 * E_nu)/(del_m12**2)
    osc_prob_pos = del_m12 * np.sqrt(((A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2))
    osc_prob_neg = del_m12 * np.sqrt(((-A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2))
    

    if R_sun:
        return rad, osc_prob_pos, osc_prob_neg
    else:
        return N_e, osc_prob_pos, osc_prob_neg
    
# Units from: https://sites.google.com/view/bentleyphysics/MSWflav

# Fermi coupling constant in GeV^-2
G_f = constants.value('Fermi coupling constant')
# speed of light in m/s
c = constants.value('speed of light in vacuum')
hbar = constants.value('Planck constant over 2 pi in eV s')
hbar_GeVs = hbar / 1e9
c_cm_s = c * 1e2
hbarcGeVcm = hbar_GeVs * c_cm_s
root2 = math.sqrt(2)

def prob(which = "greater", sig_m31_2 = 2.524e-3, sig_m21_2 = 8e-5, length = np.linspace(0.0, 1500, 10000), energy = .4, sig = 0):

    alpha = np.abs(sig_m21_2/sig_m31_2)

    sin_sq_2_theta23 = 1.0 #0.5 # 
    theta23 = np.arcsin(np.sqrt(sin_sq_2_theta23))/2

    sin_sq_2_theta12 = 0.87 #0.306 ##
    theta12 = np.arcsin(np.sqrt(sin_sq_2_theta12))/2

    sin_sq_2_theta13 =  0.1#0.087#
    theta13 = np.arcsin(np.sqrt(sin_sq_2_theta13))/2
    
    x = np.sin(theta23) * np.sin(2 * theta13)
    y = alpha * np.cos(theta23) * np.sin(2*theta12)

    delta_31 =  (sig_m31_2 * length) / (4 * energy) 
    delta_21 =  (sig_m21_2 * length) / (4 * energy)
    
    n_e = 2.5 * 6.022e22
    A = energy * 2.0 * root2 * G_f * n_e * hbarcGeVcm**3 


    a_hat = np.abs(A / sig_m31_2)

    wave_f = (np.sin((1 - a_hat)* delta_31)) / (1 - a_hat)
    wave_f_bar = (np.sin((1 + a_hat)*delta_31)) / (1 + a_hat)
    wave_g = (np.sin(a_hat * delta_31)/ a_hat)

    
    def list_det(a, b):
        if isinstance(a, float) or isinstance(a, int):
            return b 
        else:
            return a
    
    if which == "greater":
        return list_det(energy, length), ((x**2) * (wave_f**2)) + (2*x*y*wave_f*wave_g*(np.cos(sig)*np.cos(delta_31) - np.sin(sig)*np.sin(delta_31))) + ((y**2) * (wave_g**2))
    elif which == "vac":
        return list_det(energy, length), (0.1 * (np.sin(theta23)**2) * (np.sin(delta_31)**2))
    else:
        return list_det(energy, length), ((x**2) * (wave_f_bar**2)) - (2*x*y*wave_f_bar*wave_g*(np.cos(sig)*np.cos(delta_31) + np.sin(sig)*np.sin(delta_31))) + ((y**2) * (wave_g**2))

