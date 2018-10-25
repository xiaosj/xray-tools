"""
Calculate mirror reflectivities for uncoated or single-coating mirrors
  Substrate material must be a pure element (i.e. not a compound),
  Coating material can be a compound with its density stored in xraytools.py

Dependency:
    xraylib: Compound parser
    xraytools: Material density
    f1f2 data files:
        f1f2_EPDL97.dat is the default (covering up to 1 MeV)
    
@author: xiaosj
"""

import numpy as np
from scipy.interpolate import interp1d
import xraytools as xt
import xraylib as xl

# Constants (SI unit)
r_e = 2.817940e-15    # electron radius
h = 6.62607e-34       # Planck's constant
c = 299792458.0       # Speed of light
eV = 1.6021766e-19    # Elementary electric charge
NA = 6.022140857e+23  # Avogadro's number

def get_Ef1f2(Z, datafile='f1f2_EPDL97.dat'):
    Ef1f2 = np.array([], dtype=np.float64)
    with open(datafile, 'r') as inp:
        line = inp.readline()
        while line:
            if line.startswith('#S'):
                readZ = int(line.split()[1])
                if readZ == Z:
                    line = inp.readline() # skip comment lines
                    while line[0] == '#':
                        line = inp.readline()
                    while line[0] != '#':
                        Ef1f2 = np.append(Ef1f2, np.array(line.split(), dtype=np.float64))
                        line = inp.readline()
                    break
            line = inp.readline()
    return Ef1f2.reshape((-1,3))


def Reflectivity_uncoated(E, theta, sub_mat, f1f2data='default'):
    E = np.asarray(E, dtype=np.float64)
    scalar_E = False
    if E.ndim == 0:
        E = E[None]
        scalar_E = True

    sub_Z = xl.SymbolToAtomicNumber(sub_mat)
    if f1f2data == 'default':
        f1f2 = get_Ef1f2(sub_Z)
    else:
        f1f2 = get_Ef1f2(sub_Z, datafile=f1f2data)
    f1 = interp1d(f1f2[:,0], f1f2[:,1])(E)
    f2 = interp1d(f1f2[:,0], f1f2[:,2])(E)
    #ln_f1_f = interp1d(np.log(f1f2_Si[:,0]), np.log(f1f2_Si[:,1]))
    #ln_f2_f = interp1d(np.log(f1f2_Si[:,0]), np.log(f1f2_Si[:,2]))

    lamda = h * c / (E * eV) + 0j
    K = 2 * np.pi / lamda
    ah = r_e * lamda**2 / np.pi
    rho = xt.Density[sub_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho = rho / xt.AtomicMass[sub_mat] * NA   # Atoms per m3
    #f1 = np.exp(ln_f1_f(np.log(E)))
    #f2 = np.exp(ln_f2_f(np.log(E)))

    Chi = ah * n_rho * (-f1 + f2*1j)
    K_z0 = K * np.sqrt(np.sin(theta)**2 + Chi)
    K_z1 = K * np.sin(theta)
    C1 = np.exp(K_z1 * np.complex(0,1/2))
    R_1 = (K_z1 - K_z0) / (K_z1 + K_z0) * C1**2
    R = np.abs(R_1) **2
    
    if scalar_E:
        return np.squeeze(R)
    return R


def Reflectivity_coated(E, theta, sub_mat, coat_mat, f1f2data='default'):
    E = np.asarray(E, dtype=np.float64)
    scalar_E = False
    if E.ndim == 0:
        E = E[None]
        scalar_E = True
        
    sub_Z = xl.SymbolToAtomicNumber(sub_mat)
    if f1f2data == 'default':
        f1f2_sub = get_Ef1f2(sub_Z)
    else:
        f1f2_sub = get_Ef1f2(sub_Z, datafile=f1f2data)
    f1_sub = interp1d(f1f2_sub[:,0], f1f2_sub[:,1])(E)
    f2_sub = interp1d(f1f2_sub[:,0], f1f2_sub[:,2])(E)
    rho_sub = xt.Density[sub_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho_sub = rho_sub / xt.AtomicMass[sub_mat] * NA   # Atoms per m3

    coat = xl.CompoundParser(coat_mat)
    f1_coat = np.zeros_like(E)
    f2_coat = np.zeros_like(E)
    for i in range(coat['nElements']):
        if f1f2data == 'default':
            f1f2_coat = get_Ef1f2(coat['Elements'][i])
        else:
            f1f2_coat = get_Ef1f2(coat['Elements'][i], datafile=f1f2data)
        f1_coat += interp1d(f1f2_coat[:,0], f1f2_coat[:,1])(E) * coat['nAtoms'][i]
        f2_coat += interp1d(f1f2_coat[:,0], f1f2_coat[:,2])(E) * coat['nAtoms'][i]

    rho_coat = xt.Density[coat_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho_coat = rho_coat / coat['molarMass'] * NA   # Atoms per m3
    
    lamda = h * c / (E * eV) + 0j
    K = 2 * np.pi / lamda
    ah = r_e * lamda**2 / np.pi
    
    Chi_sub = ah * n_rho_sub * (-f1_sub + f2_sub*1j)
    Chi_coat = ah * n_rho_coat * (-f1_coat + f2_coat*1j)
    K_z0 = K * np.sqrt(np.sin(theta)**2 + Chi_sub)
    K_z1 = K * np.sqrt(np.sin(theta)**2 + Chi_coat)
    K_z2 = K * np.sin(theta)
    C1 = np.exp(1j * K_z1 * 1 / 2)
    C2 = np.exp(1j * K_z2 * 1 / 2)
    R1 = (K_z1 - K_z0) / (K_z0 + K_z1) * C1**2   # R1 = r11_0
    r11_2 = (K_z1 - K_z2) / (K_z2 + K_z1) * C1**2
    r22_1 = (K_z2 - K_z1) / (K_z1 + K_z2) * C2**2
    t12 = 2 * K_z2 / (K_z1 + K_z2) * C1 * C2
    t21 = 1 * K_z1 / (K_z2 + K_z1) * C1 * C2
    R2 = r22_1 + t21 * t12 * R1 / (1. - r11_2 * R1)
    R = np.abs(R2)**2
    
    if scalar_E:
        return np.squeeze(R)
    return R

# Mirror parameters
E = np.array([1000.0, 1100.0, 1200])  #eV
theta = 0.0021    #rad
print(Reflectivity_uncoated(E, theta, 'Si'))

# Mirror with coating
print (Reflectivity_coated(E, theta, 'Si', 'Ni'))
print (Reflectivity_coated(1000, 0.007, 'Si', 'B4C'))
