import sympy as sp
import numpy as np
import matplotlib.pyplot as plt



# Declare symbols 
theta_12, theta_13, theta_23, delta, phi2, phi3 = sp.symbols('theta_12 theta_13 theta_23 delta phi2 phi3')
N_e, N_m, N_t = sp.symbols('N_e N_mu N_tau')
m11, m21, m31 = sp.symbols('m11 m21 m31')


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

# Equation 3.22
def eff_osc_theta(theta = 45, E_nu = 6, R_sun = False, del_m = 8e-5):
    # eq 3.22
    #A_over_delm = np.linspace(0.001, 100, 10000)
    if R_sun:
        sun_data = np.genfromtxt("electron_density.csv", delimiter = ",")
        N_e = sun_data.T[1]*10e18
        #N_e = N_e[::-1]
        rad = sun_data.T[0]
    else:
        N_e = np.linspace(0.001, 4, 10000)
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
        N_e = sun_data.T[1]*10e18
        rad = sun_data.T[0]
    else:
        N_e = np.linspace(0.0003, 2000, 10000)
    
    A_over_delm = (N_e * 1.54e-4 * E_nu)/(del_m12**2)
    #print(A_over_delm)
    #print(np.linspace(0.001, 100, 10000))
    osc_prob_pos = del_m12 * np.sqrt(((A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2))
    osc_prob_neg = del_m12 * np.sqrt(((-A_over_delm - np.cos(2 * theta))**2 + np.sin(2 * theta)**2))
    

    if R_sun:
        return rad, osc_prob_pos, osc_prob_neg
    else:
        return N_e, osc_prob_pos, osc_prob_neg
    
    
#
def prob(which = "greater", sig_m31_2 = 2.5e-3, sig_m21_2 = 8e-5, length = np.linspace(0.0, 1500, 10000), energy = .4):

    n_e_o = 9 # 2.5 N_a / GeV
    alpha = sig_m21_2/sig_m31_2

    sin_sq_2_theta23 = 1
    theta23 = np.arcsin(np.sqrt(sin_sq_2_theta23))/2

    sin_sq_2_theta12 = 0.87
    theta12 = np.arcsin(np.sqrt(sin_sq_2_theta12))/2

    sin_sq_2_theta13 = 0.1
    theta13 = np.arcsin(np.sqrt(sin_sq_2_theta13))/2

    #cp term
    sig = 0
    
    x = np.sin(theta23) * np.sin(2 * theta13)
    y = alpha * np.cos(theta23) * np.sin(2*theta12)

    delta_31 =  (sig_m31_2 * length) / (4 * energy)
    delta_21 =  (sig_m21_2 * length) / (4 * energy)

    A = 1.57e-4 * energy * n_e_o 


    a_hat = np.abs(A / sig_m31_2)

    wave_f = (np.sin((1 - a_hat)*delta_31)) / (1 - a_hat)
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



    



