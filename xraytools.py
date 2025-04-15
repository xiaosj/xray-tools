# utility functions for xray studies

import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import dawsn, expi
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import xraylib as xl
import requests


def goodID(matID, item=None):
    """ Check if the material ID is defined:
        item = None,    in xraylib
             'density', in density dictionary
             'Cp',      in heat capacity dictionary
             'crystal', in crystal dictionary
    """
    good = True
    id = matID
    if item == 'density':
        if matID not in Density:
            print('Cannot find the density of ' + matID)
            good = False
    if item == 'Cp':
        if matID not in specificHeatParams and matID not in specificHeat:
            print('Cannot find the specific heat of ' + matID)
            good = False
    if item == 'crystal':
        if matID not in latticeType:
            print('Cannot find the lattice type of ' + matID)
            good = False
        if matID not in latticeParameters:
            print('Cannot find the lattice parameters of ' + matID)
            good = False
    else:
        if matID in alias:
            id = alias[matID]
        else:
            id = matID

        try:
            _ = xl.CompoundParser(id)
        except:  # noqa: E722
            print('Cannot find ' + matID + ' in xraylib')
            good = False
   
    if not good:
        raise Exception('STOP: Missing material properties.')

    return id


def eV(E):
    """ Returns photon energy in eV if specified in eV or keV """
    if np.max(E) < 100:
        return E * 1000
    else:
        return E


def keV(E):
    """ Returns photon energy in keV if specified in eV or keV """
    if np.min(E) >= 100:
        return E / 1000
    else:
        return E
   

def lam(E):
    """ Computes photon wavelength in m
        E is photon energy in eV
    """
    return (12398.4/E)*1e-10

def lam2E(l):  # noqa: E741
    """ Computes photon energy in eV
        l is photon wavelength in m
    """
    E=12398.4/(l*u['ang'])
    return E

def lam2f(l):  # noqa: E741
    """ Computes the photon frequency in Hz
        l is photon wavelength in m
    """
    f=c['c']/l
    return f
    
def f2lam(f):
    """ Computes the photon wavelength in m
        f is the photon frequency in Hz
    """
    l=c['c']/f  # noqa: E741
    return l
 
def f2E(f):
    """ Computes the photon in energy in eV
        f is the photon frequency in Hz
    """
    E=c['h']*f*u['eV']
    return E
 
def E2f(E):
    """ Computes the photon frequency in Hz
        E is photon energy in eV or keV
    """
    f=E/c['h']/u['eV']
    return f

def sind(A):
    """ Sin of an angle specified in degrees """
    Arad = np.deg2rad(A)
    x = np.sin(Arad) 
    return x
 
def cosd(A):
    """ Cos of an angle specified in degrees """
    Arad = np.deg2rad(A)
    x = np.cos(Arad) 
    return x

def tand(A):
    """ Tan of an angle specified in degrees """
    Arad = np.deg2rad(A)
    x = np.tan(Arad) 
    return x
 
def asind(x):
    """ Arcsin in degrees """
    A = np.arcsin(x)
    A = np.rad2deg(A) 
    return A
 
def acosd(x):
    """ Arccos in degrees """
    A = np.arccos(x)
    A = np.rad2deg(A) 
    return A

def atand(x):
    """ Arctan in degrees """
    A = np.arctan(x)
    A = np.rad2deg(A) 
    return A


def defaultDensity(matID):
    """ Default material density """
    mat = goodID(matID)
    if matID in Density:
        density = Density[matID]
    else:
        density = Density[mat]
    return density


def atomWeight(matID):
    """ Return the average atom weight in g/mol """
    mat = goodID(matID)
    compound = xl.CompoundParser(mat)
    mass = 0.0
    for i in range(compound['nElements']):
        mass += xl.AtomicWeight(compound['Elements'][i]) * compound['massFractions'][i]
    return mass


def molarMass(matID):
    """ Return the molar mass as g/mol """
    mat = goodID(matID)
    compound = xl.CompoundParser(mat)
    return atomWeight(matID) * compound['nAtomsAll']


def mu(matID, keV, density=None):
    """ Calculate the mass attenuation coefficients (1/m) at given energies
          keV: energy in keV (vectorized)
          density in g/cm3, None=default density
    """
    mat = goodID(matID)
    if density is None:
        density = defaultDensity(matID)
    if np.isscalar(keV):
        energies = np.array([keV], dtype=np.double)
    else:
        energies = np.array(keV, dtype=np.double)
    _mu = np.array([xl.CS_Total_CP(mat, eng) * density * u['cm'] for eng in energies])
    if np.isscalar(keV):
        # return np.asscalar(_mu), numpy-asscalar() has been deprecated since NumPy 1.16
        return _mu.item()
    else:
        return _mu


def mu_en(matID, keV, density=None):
    """ Calculate the mass energy-absorption coefficients (1/m) at given energies
          keV: energy in keV (vectorized)
          density in g/cm3, None=default density
    """
    mat = goodID(matID)
    if density is None:
        density = defaultDensity(matID)
    if np.isscalar(keV):
        energies = np.array([keV], dtype=np.double)
    else:
        energies = np.array(keV, dtype=np.double)
    _mu = np.array([xl.CS_Energy_CP(mat, eng) * density * u['cm'] for eng in energies])
    if np.isscalar(keV):
        return _mu.item()
    else:
        return _mu


def attenuationLength(matID, keV, density=None):
    """ Calculate the attenuation length (m) at given energies
          E in keV (vectorized)
          density in g/cm3, None=default density
    """
    return 1.0 / mu(matID, keV, density)


def transmission(matID, t, keV, density=None):
    """ Calculate the transmission at given thickness and energies
          t: thickness in m
          E: energies in keV (vectorized)
          density: in g/cm3, None=default density
    """
    return np.exp(-mu(matID, keV, density) * t)


def eVatom(matID, keV, mJ, rms_mm, density=None):
    """ Calculate the eV/atom at given energies
          keV: energies in keV (vectorized)
          mJ: pulse energy in mJ (vectorized)
          rms_mm: beam size radius in mm (vectorized)
             -- E, mJ, rms_mm must match if more than one are vectorized
          density: in g/cm3, None=default density
    """
    if density is None:
        density = defaultDensity(matID)
    attL = attenuationLength(matID, keV, density)
    EdensityJcm3 = mJ/1000 / (2 * np.pi * attL*u['cm'] * (rms_mm*0.1)**2)
    atomVolcm3 = atomWeight(matID) / c['NA'] / density
    natoms = xl.CompoundParser(matID)['nAtomsAll']
    return EdensityJcm3 * atomVolcm3 / 1.6e-19 / natoms


def eVatom_en(matID, keV, mJ, rms_mm, density=None):
    """ Calculate the eV/atom at given energies with mu_en
          keV: energies in keV (vectorized)
          mJ: pulse energy in mJ (vectorized)
          rms_mm: beam size radius in mm (vectorized)
             -- E, mJ, rms_mm must match if more than one are vectorized
          density: in g/cm3, None=default density
    """
    if density is None:
        density = defaultDensity(matID)
    attL = 1.0 / mu_en(matID, keV, density)
    EdensityJcm3 = mJ/1000 / (2 * np.pi * attL*u['cm'] * (rms_mm*0.1)**2)
    atomVolcm3 = atomWeight(matID) / c['NA'] / density
    natoms = xl.CompoundParser(matID)['nAtomsAll']
    return EdensityJcm3 * atomVolcm3 / 1.6e-19 / natoms


def eVatom_keV_plot(matID, keV, mJ, rms_mm, density=None, logx=False, logy=True):
    if not isinstance(matID, list):
        matID = [matID]
    
    xlist = None
    inp_names = ['keV', 'mJ', 'rms_mm']
    for inp, name in zip([keV, mJ, rms_mm], inp_names):
        if isinstance(inp, (np.ndarray, list)):
            if xlist is None:
                xlist = inp
                xlabel = name
            else:
                raise Exception('Only one variable of keV, mJ and rms_mm can be an array or list.')
    
    if xlist is None:
        for mat in matID:
            eVatoms = eVatom(mat, keV, mJ, rms_mm, density)
            print('{:s}: {:.2g} eV/atom'.format(mat, eVatoms))
    else:
        plt.figure(figsize=(5,3), dpi=100, facecolor='white')
        for mat in matID:
            eVatoms = eVatom(mat, keV, mJ, rms_mm, density)
            if logx:
                if logy:
                    plt.loglog(xlist, eVatoms, label=mat)
                else:
                    plt.semilogx(xlist, eVatoms, label=mat)
                    plt.ylim(0)
                plt.gca().xaxis.set_major_formatter(ScalarFormatter())
            else:
                if logy:
                    plt.semilogy(xlist, eVatoms, label=mat)
                else:
                    plt.plot(xlist, eVatoms, label=mat)
                    plt.ylim(0)

        if xlabel=='keV':
            plt.title('Dose from {:.1f} mJ {:.2g} mm (rms) pulse'.format(mJ, rms_mm))
        elif xlabel=='mJ':
            plt.title('Dose from {:.1f} keV {:.2g} mm (rms) pulse'.format(keV, rms_mm))
        elif xlabel=='rms_mm':
            plt.title('Dose from {:.1f} keV {:.1f} mJ pulse'.format(keV, mJ))
        
        plt.xlabel(xlabel)
        plt.ylabel('eV/atom')
        plt.grid(axis='both', which='both')
        plt.legend()
        plt.tight_layout()
        plt.show()


def drillSpeed(matID, power_W, FWHM_mm):
    """ Return the material drill speed (mm/s) based on vaporization heat """
    vaporH = {  # kJ/mol
        # spec heat from room temperature to melting + latent heat of fusion + spec heat from melting to boiling + latent heat of vaporization
        # refer https://webbook.nist.gov/chemistry/ for heat capacity, this tool also has some data in specificHeatParams
        'Cu': 29.67 + 13.26 + 48.77 + 300,
        'Fe': 49.62 + 13.81 + 60.89 + 340,
        'W' :119.39 + 52.31 + 89.19 + 774,
        'Mo': 89.9  + 37.5  + 71.9  + 598,
        'Al': 17.97 + 10.71 + 59.06 + 284,
        'C':  29.23 +                 715  # Graphite: use 8.12 J/mol/K from webbook (a small heat capacity from Sheindlin 1972 for 263-3650K); sublimation from 3900 K at normal pressure (no melting), https://en.wikipedia.org/wiki/Heats_of_vaporization_of_the_elements_(data_page)
    }
    if matID not in vaporH.keys():
        raise ValueError(f'No vaporization data for {matID}: available in {vaporH.keys()}')
    
    mol_mmD = 2 * np.pi * (FWHM_mm/2.355)**2 / 1000 * defaultDensity(matID) / molarMass(matID)
    return power_W / (mol_mmD * vaporH[matID] * 1000)


def drillTime(matID, thickness_mm, W, FWHM_mm):
    """ Return the material drill time based on vaporization heat """
    return thickness_mm / drillSpeed(matID, W, FWHM_mm)


def getCp(matID, T):
    goodID(matID, 'Cp')
    if matID in specificHeatParams:
        TT = T / 1000
        A = specificHeatParams[matID]
        return A[0] + A[1] * TT + A[2] * TT**2 + A[3] * TT**3 + A[4] / TT**2
    else:
        raise Exception(f'No heat capacity data for {matID}')

def getIntCp(matID, T1_K, T2_K):
    """ Get integrated heat capacity in J/mol
        * T1, T2 in K
    """
    goodID(matID, 'Cp')
    if matID in specificHeatParams:
        A = specificHeatParams[matID]
        return intCp2(T1_K, T2_K, A)
    else:
        raise Exception(f'Wrong material {matID}: Only support materials within specificHeatParams.')

def intCp(T, A):
    TT = T / 1000
    return 1000 * (A[0] * TT + A[1]/2 * TT**2 + A[2]/3 * TT**3 + A[3]/4 * TT**4 - A[4] / TT)

def intCp2(T1, T2, A):
    return intCp(T2, A) - intCp(T1, A)
    
def pulseT(matID, keV, mJ, rms_mm, density=None, baseT=298.15):
    """ Calculate the instant temperature (K) from single pulse
        * keV: energies in keV (vectorized)
        * mJ: pulse energy in mJ (vectorized)
        * rms_mm: beam size radius in mm (vectorized)
            -- Size E, mJ, rms_mm must match if more than one are vectorized
        * density: in g/cm3, None=default density
        * baseT: base temperature (K), 298 K in default
    """
    if density is None:
        density = defaultDensity(matID)
    attL = attenuationLength(matID, keV, density)
    EdensityJcm3 = mJ / 1000. / (2. * np.pi * attL*u['cm'] * (rms_mm*0.1)**2)
    EdensityJmol = EdensityJcm3 / (density / molarMass(matID))

    # build temperature as a fuction of J/cm3
    goodID(matID, 'Cp')
    if matID in specificHeatParams:
        CpParam = specificHeatParams[matID]
        if type(CpParam) is tuple:
            CpParam = [[baseT, meltPoint[matID], CpParam]]

        if baseT < CpParam[0][0]:
            CpParam[0][0] = baseT
        if baseT > CpParam[-1][1]:
            raise Exception('Base temperature cannot be higher than the melting point of {:s} ({:d} K)'.format(matID, int(CpParam[-1][1])))
        nzone = len(CpParam)
        for i in range(nzone):
            if baseT >= CpParam[0][0] and baseT <= CpParam[0][1]:
                CpParam[0][0] = baseT
                break
            else:
                del CpParam[0]
        
        dT = 10
        T = np.array([])
        Jmol = np.array([])
        J0 = 0
        for i in range(len(CpParam)):
            Ti = np.arange(CpParam[i][0], CpParam[i][1], dT)
            Ji = intCp2(CpParam[i][0], Ti, CpParam[i][2]) + J0
            T = np.concatenate((T, Ti))
            Jmol = np.concatenate((Jmol, Ji))
            if i < len(CpParam) - 1:
                J0 = intCp2(CpParam[i][0], CpParam[i+1][0], CpParam[i][2])
        Jmol_base = 0.0

    else:  # use itegrated specific heat directly
        T = specificHeat[matID][0]
        Jmol = specificHeat[matID][2]
        Jmol_base = interp1d(T, Jmol)(baseT) 

    T_Jmol = interp1d(Jmol, T, fill_value='extrapolate')  # T as function of J/mol

    return T_Jmol(EdensityJmol + Jmol_base) - baseT


def pulseTC(matID, keV, mJ, rms_mm, density=None, baseT=25):
    """ Calculate the instant temperature (C) from single pulse
        * keV: energies in keV (vectorized)
        * mJ: pulse energy in mJ (vectorized)
        * rms_mm: beam size radius in mm (vectorized)
            -- Size of E, mJ, rms_mm must match if more than one are vectorized
        * density: in g/cm3, None=default density
        * baseT: base temperature (K), 298 K in default
    """
    return pulseT(matID, keV, mJ, rms_mm, density=density, baseT=C2K(baseT))


def spectrum_cut(spectrum, eVrange=(0.0, 0.0)):
    """ Cut spectrum to a given energy range; return [0,0] if out of range """
    if eVrange[1] == 0.0:
        return spectrum
    else:
        if spectrum[-1,0] <= eVrange[0] or spectrum[0,0] >= eVrange[1]:
            return np.array([[0, 0]], dtype=np.float)
        else:
            idx1 = np.argmax(spectrum[:,0] >= eVrange[0])
            idx2 = np.argmax(spectrum[:,0] >  eVrange[1])
            if spectrum[0,0] >= eVrange[0]:
                idx1 = 0
            if spectrum[-1,0] <= eVrange[1]:
                idx2 = -1
            return spectrum[idx1:idx2]


def spectrum_eV_power_mW(spectrum, eVrange=(0.0,0.0)):
    spec = spectrum_cut(spectrum, eVrange)
    eV = spec[:,0]
    flux = spec[:,1]
    W_bin = [(flux[i]+flux[i+1]) / 2 * (eV[i+1]-eV[i]) * (eV[i+1]+eV[i]) / 2 * c['e']
              for i in range(len(eV)-1)]
    return np.sum(W_bin) * 1000


def spectrum_eV_flux(spectrum, eVrange=(0.0,0.0)):
    spec = spectrum_cut(spectrum, eVrange)
    eV = spec[:,0]
    flux = spec[:,1]
    flux_bin = [(flux[i] + flux[i+1]) / 2 * (eV[i+1] - eV[i]) for i in range(len(eV)-1)]
    return np.sum(flux_bin)


def spectrum_flux(spectrum, eVrange=(0.0,0.0), specType='bw'):
    if specType == 'eV':
        return spectrum_eV_flux(spectrum, eVrange)
    else:
        spec = np.copy(spectrum)
        spec[:,1] = spectrum[:,1] / (spectrum[:,0] * 0.001)
        return spectrum_eV_flux(spec, eVrange)


def spectrum_power_mW(spectrum, eVrange=(0.0,0.0), specType='bw'):
    if specType == 'eV':
        return spectrum_eV_power_mW(spectrum, eVrange)
    else:
        spec = np.copy(spectrum)
        spec[:,1] = spectrum[:,1] / (spectrum[:,0] * 0.001)
        return spectrum_eV_power_mW(spec, eVrange)


def spectrum_dose(spectrum, area_cm2, eVrange=(0.0,0.0), specType='bw', particle='photon'):
    """
       spectrum: [eV, flux]
       area_cm2: area of spot
       return dose in mrem/h
    """
    spec = spectrum_cut(spectrum, eVrange)
    eV = spec[:,0]
    if specType == 'bw':
        flux_eV = spec[:,1] / (eV * 0.001)
    else:
        flux_eV = spec[:,1]

    # Load flux to dose 
    path = os.path.dirname(os.path.abspath(__file__))
    f2d = np.loadtxt(path+'/f2d_' + particle)
    f2d[:,0] *= 1.0E9    # from GeV to eV
    f2d[:,1] *= 3.6E-4   # from pSv/s to mrem/h
    
    # Linear interpolation for the log-log of f2d
    logf2d = np.log(f2d)
    f_logf2d = interp1d(logf2d[:,0], logf2d[:,1],
        bounds_error=False, fill_value='extrapolate', assume_sorted=True)
    
    d = np.exp(f_logf2d(np.log(eV))) * flux_eV
    dose_bin = [(d[i] + d[i+1]) / 2 * (eV[i+1] - eV[i]) for i in range(len(eV)-1)]
    
    return np.sum(dose_bin) / area_cm2
    

def spectrum_shield(spectrum, area_cm2, matID, density=None, eVrange=(0.0,0.0), specType='bw',
                    dose_limit=0.05, dt_limit=0.01, particle='photon', verbose=0):
    """
    Retrun required shielidng in mm
       spectrum:   [eV, flux]
       area_cm2:   area of spot
       matID:      shielding material ID
       density:    density of shielding material, use the default density if not specified
       eVrange:    cut the spectrum to the given range
       specType:   type of spectrum, 'eV', flux/eV or 'bw', flux/0.1%bw
       dose_limit: shielding goal, 0.05 mrem/h by default
       dt_limit:   smallest shielding thickness step
       particle:   type of particle
       verbose:    output details if > 0
    """
    spec = spectrum_cut(spectrum, eVrange)
    al = np.array([attenuationLength(matID, eV/1000, density=density) for eV in spec[:,0]]) * 1000
    
    d = spectrum_dose(spec, area_cm2, specType=specType)
    if verbose > 0:
        print('Dose without shielding is {:.3g} mrem/h'.format(d))
    
    if d < dose_limit:
        if verbose > 0:
            print('No shielding is needed.')
        return 0.
    
    else:
        t0 = 0.
        t1 = 1.
        while True:   # find the upper bound of shielding
            spec_trans = np.copy(spec)
            spec_trans[:,1] *= np.exp(-t1/al)
            d = spectrum_dose(spec_trans, area_cm2, specType=specType)
            if d >= dose_limit:
                t0 = t1
                t1 *= 2
            else:
                break
        
        while t1 - t0 >= dt_limit:
            t = (t0 + t1) * 0.5
            spec_trans = np.copy(spec)
            spec_trans[:,1] *= np.exp(-t/al)
            d = spectrum_dose(spec_trans, area_cm2, specType=specType)
            if verbose > 0:
                print('  Dose after {:.3f} mm {:s} is {:.3f} mrem/h'.format(t, matID, d))
            if d > dose_limit:
                t0 = t
            else:
                t1 = t
    
        if verbose > 0:
            print('Required {:s} thickness is {:.3f} mm'.format(matID, t1))
    
    return t1


def plot(x, y, xlabel=None, ylabel=None, title=None, figsize=(4.5,3), logx=False, logy=False, xmin=None, xmax=None, ymin=None, ymax=None, savefig=None):
    ''' Plot function
    '''
    plt.figure(figsize=figsize, dpi=100, facecolor='white')
    
    if not logx and not logy:
        plt.plot(x, y)
    if not logx and logy:
        plt.semilogy(x, y)
    if logx and not logy:
        plt.semilogx(x, y)
    if logx and logy:
        plt.loglog(x, y)

    if xmin is not None:
        plt.xlim(left=xmin)
    else:
        plt.xlim(left=x[0])
    
    if xmax is not None:
        plt.xlim(right=xmax)
    else:
        plt.xlim(right=x[-1])

    if ymin is not None:
        plt.ylim(bottom=ymin)
    else:
        if not logy:
            plt.ylim(bottom=0)

    if ymax is not None:
        plt.ylim(top=ymax)

    if title is not None:
        plt.title(title, fontsize=11)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)

    plt.grid(True)
    plt.tight_layout()

    if savefig is not None:
        if savefig[-4:] != '.png':
            savefig = savefig + '.png'
        plt.savefig(savefig)

    plt.show()


def C2K(degC):
    ''' Convert Celsius to Kelvin '''
    return degC + 273.15


def K2C(degK):
    ''' Convert Kelvin to Celsius '''
    return degK - 273.15


def temp3d_uniform(x, y, z, power, q, lx, ly, T_surface_avg,
                   width, height, thickness, K, T_env=300.0, verbose=0):
    ''' Calculate the temperature for radiation cooling with uniform input power '''
    Lx = width * 0.5
    Ly = height * 0.5
    h = thickness

    h_eff = power / (width * height * (T_surface_avg - T_env))
    t = 0.00001
    Ktwid = h_eff * t
    S = 0.5 * (1 - Ktwid / K)
    R = 0.5 * (1 + Ktwid / K)
    
    nx = int(Lx / lx)
    if nx < 20:
        nx = 20
    ny = int(Ly / ly)
    if ny < 20:
        ny = 20
    
    T = q * lx * ly / (2 * Lx * Ly) * (t / Ktwid + (h - z) / K)

    if verbose > 1:
        print(f'power = {power}, q = {q}, lx = {lx}, ly = {ly}, T_surface_avg = {T_surface_avg}')
        print(f'Lx = {Lx}, Ly = {Ly}, h = {h}')
        print(f'h_eff = {h_eff:.3g}, Ktwid = {Ktwid:.3g}')
        print(f'S = {S}, R = {R}')
        print(f'nx = {nx}, ny = {ny}')
        
    for m in range(1, nx+1):
        omega = m * np.pi / Lx
        gamma = omega
        temp  = -2 * q * ly * np.sin(omega * lx) * np.cos(omega * x) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= omega * gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for n in range(1, ny+1):
        phi   = n * np.pi / Ly
        gamma = phi
        temp  = -2 * q * lx * np.sin(phi * ly) * np.cos(phi * y) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for m in range(1, nx+1):
        for n in range(1, ny+1):
            omega = m * np.pi / Lx
            phi   = n * np.pi / Ly
            gamma = np.sqrt(omega * omega + phi * phi)
            temp  = -4 * q * np.sin(omega * lx) * np.sin(phi * ly) * np.cos(omega * x) * np.cos(phi * y)
            temp *= S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z))
            temp /= omega * gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
            T += temp
        
    return T + T_env

def temp3d_1dGaussian(x, y, z, power, q, lx, sigma_y, T_surface_avg,
                   width, height, thickness, K, T_env=300.0, verbose=0):
    ''' Calculate the temperature for radiation cooling with 1D Gaussian input power '''
    r2pi = np.sqrt(2 * np.pi)    
    Lx = width * 0.5
    Ly = height * 0.5
    h = thickness

    h_eff = power / (width * height * (T_surface_avg - T_env))
    t = 0.00001
    Ktwid = h_eff * t
    S = 0.5 * (1 - Ktwid / K)
    R = 0.5 * (1 + Ktwid / K)
    
    nx = int(Lx / lx)
    if nx < 20:
        nx = 20
    ny = int(Ly / sigma_y)
    if ny < 20:
        ny = 20
    
    T = r2pi * sigma_y * q * lx / (2 * Lx * Ly) * (t / Ktwid + (h - z) / K)

    if verbose > 1:
        print(f'power = {power}, q = {q}, lx = {lx}, sigma_y = {sigma_y}, T_surface_avg = {T_surface_avg}')
        print(f'Lx = {Lx}, Ly = {Ly}, h = {h}')
        print(f'h_eff = {h_eff:.3g}, Ktwid = {Ktwid:.3g}')
        print(f'S = {S}, R = {R}')
        print(f'nx = {nx}, ny = {ny}')
        
    for m in range(1, nx+1):
        omega = m * np.pi / Lx
        gamma = omega
        temp  = -r2pi * sigma_y * q * np.sin(omega * lx) * np.cos(omega * x) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= omega * gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for n in range(1, ny+1):
        phi   = n * np.pi / Ly
        gamma = phi
        temp  = -r2pi * sigma_y * q * lx * np.exp(-0.5 * (phi * sigma_y)**2) * np.cos(phi * y) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for m in range(1, nx+1):
        for n in range(1, ny+1):
            omega = m * np.pi / Lx
            phi   = n * np.pi / Ly
            gamma = np.sqrt(omega * omega + phi * phi)
            temp  = -2 * r2pi * sigma_y * q * np.sin(omega * lx) * np.exp(-0.5 * (phi * sigma_y)**2) * np.cos(omega * x) * np.cos(phi * y)
            temp *= S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z))
            temp /= omega * gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
            T += temp
        
    return T + T_env


def temp3d_2dGaussian(x, y, z, power, q, sigma_x, sigma_y, T_surface_avg,
                      width, height, thickness, K, T_env=300.0, verbose=0):
    ''' Calculate the temperature for radiation cooling with 2D Gaussian input power '''
    Lx = width * 0.5
    Ly = height * 0.5
    h = thickness

    h_eff = power / (width * height * (T_surface_avg - T_env))
    t = 0.00001
    Ktwid = h_eff * t
    S = 0.5 * (1 - Ktwid / K)
    R = 0.5 * (1 + Ktwid / K)
    
    nx = int(Lx / sigma_x)
    if nx < 20:
        nx = 20
    ny = int(Ly / sigma_y)
    if ny < 20:
        ny = 20
    
    piqxy = np.pi * q * sigma_x * sigma_y
    T = piqxy / (2 * Lx * Ly) * (t / Ktwid + (h - z) / K)

    if verbose > 1:
        print(f'power = {power}, q = {q}, lx = {sigma_x}, sigma_y = {sigma_y}, T_surface_avg = {T_surface_avg}')
        print(f'Lx = {Lx}, Ly = {Ly}, h = {h}')
        print(f'h_eff = {h_eff:.3g}, Ktwid = {Ktwid:.3g}')
        print(f'S = {S}, R = {R}')
        print(f'nx = {nx}, ny = {ny}')
    
    for m in range(1, nx+1):
        omega = m * np.pi / Lx
        gamma = omega
        temp  = -piqxy * np.exp(-0.5 * (omega * sigma_x)**2) * np.cos(omega * x) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for n in range(1, ny+1):
        phi   = n * np.pi / Ly
        gamma = phi
        temp  = -piqxy * np.exp(-0.5 * (phi * sigma_y)**2) * np.cos(phi * y) * (S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z)))
        temp /= gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
        T += temp
    
    for m in range(1, nx+1):
        for n in range(1, ny+1):
            omega = m * np.pi / Lx
            phi   = n * np.pi / Ly
            gamma = np.sqrt(omega * omega + phi * phi)
            temp  = -piqxy * np.exp(-0.5 * ((omega * sigma_x)**2 + (phi * sigma_y)**2)) * np.cos(omega * x) * np.cos(phi * y)
            temp *= S * np.sinh(gamma * (z + t - h)) + R * np.sinh(gamma * (h + t - z))
            temp /= gamma * K * Lx * Ly * (S * np.cosh(gamma * (t-h)) - R * np.cosh(gamma * (t + h)))
            T += temp
        
    return T + T_env


def radiation_cooling(power:float, beam_mode:int, q:float, lx:float, ly:float, 
                      target_width:float, target_height:float, target_thickness:float,
                      target_K:float, target_emissivity:float,
                      env_emissivity:float, T_env=300.0, back_rad_only=False, verbose=0):
    ''' Calculate the temperature from radiation cooling (Not Vectorized)
    Input:
        power: beam power in Watts
        beam_mode: 0 - uniform beam; 1 - 1D Gaussian on Y; 2 - 2D Gaussian
        q: power density in W/mm2
        lx: half beam size on X, mm (rms for Gaussian beam)
        ly: half beam size on Y, mm (rms for Gaussian beam)
        target_width: target width in mm
        target_height: target height in mm
        target_thickness: target thickness in mm
        target_K: target thermal conductivity in W/mm/K
        target_emissivity: target emissivity
        env_emissivity: environment emissivity
        T_env: environment temperature in K
        back_rad_only: whether radiation is on the back side only
        verbose: verbose level for output
        
    Return:
        Tmax: the maximum temperature on target (K)
        T_surface_avg, T000, Tmin, T_env: facilitate temperatures (K)
        ratio: if it is >>1, need to set back_rad_only=True
    '''
    # Step 1: estimate an average surface temperature using radiation cooling of the entire mask surface and calculate h_eff using this ave temperature with radiation cooling of mask back surface only
    if back_rad_only:
        area = target_width * target_height
    else:
        area = 2 * (target_width * target_height + target_width * target_thickness + target_height * target_thickness)
    emissivity_eff = target_emissivity * env_emissivity / (target_emissivity + env_emissivity - target_emissivity * env_emissivity)
    sigma = 5.67e-14  # Stefan-Boltzmann constant in W/mm2/K4
    T_surface_avg = (power / (sigma * area * emissivity_eff) + T_env**4) **0.25
    
    # Step 2: calculate the temperature distribution in the mask hot wall (front layer) using h_eff and the analytic two layer model of M302
    if beam_mode == 0:
        T000 = temp3d_uniform(0, 0, 0, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
        Tmin = temp3d_uniform(target_width/2, target_height/2, target_thickness, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
    elif beam_mode == 1:
        T000 = temp3d_1dGaussian(0, 0, 0, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
        Tmin = temp3d_1dGaussian(target_width/2, target_height/2, target_thickness, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
    elif beam_mode == 2:
        T000 = temp3d_2dGaussian(0, 0, 0, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
        Tmin = temp3d_2dGaussian(target_width/2, target_height/2, target_thickness, power, q, lx, ly, T_surface_avg,
                    target_width, target_height, target_thickness, target_K, T_env=T_env, verbose=verbose)
    else:
        raise Exception(f'Wrong beam mode [{beam_mode}]: must be 0, 1 or 2')

    Tmax = T_surface_avg + T000 - Tmin
    
    # Check if it is conservative
    # dT1 = T000 - Tmin
    # ratio1 = 2 * ((Tmax - dT1)**4 - T_env**4) / ((target_melt - dT1)**4 - T_env**4)
    ratio = (T_surface_avg - T_env) / (Tmax - Tmin)
    
    if verbose > 0:
        print(f'Tmax = {K2C(Tmax):.1f}, T000 = {K2C(T000):.1f}, Tmin = {Tmin:.1f} deg-C')

    return (Tmax, T_surface_avg, T000, Tmin, T_env, ratio)


def dSpace(ID,hkl):
    """ Computes the d spacing (m) of the specified material and reflection 
        ID is chemical fomula : 'Si'
        hkl is the reflection : (1,1,1)
    """
    ID = goodID(ID, item='crystal')
    h=hkl[0]
    k=hkl[1]
    l=hkl[2]  # noqa: E741

    lp=latticeParameters[ID]
    a=lp[0]/u['ang']
    b=lp[1]/u['ang']
    c=lp[2]/u['ang']
    alpha=lp[3]
    beta=lp[4]
    gamma=lp[5]

    ca=cosd(alpha)
    cb=cosd(beta)
    cg=cosd(gamma)
    sa=sind(alpha)
    sb=sind(beta)
    sg=sind(gamma)

    invdsqr=1/(1+2*ca*cb*cg-ca**2-cb**2-cg**2)*(h**2*sa**2/a**2 + k**2*sb**2/b**2 + l**2*sg**2/c**2 +
      2*h*k*(ca*cb-cg)/a/b+2*k*l*(cb*cg-ca)/b/c+2*h*l*(ca*cg-cb)/a/c)
      
    d=invdsqr**-0.5
    return d


def BraggAngle(ID,hkl,E=None):
    """ Computes the Bragg angle (deg) of the specified material, reflection and photon energy 
        ID is chemical fomula : 'Si'
        hkl is the reflection : (1,1,1)
        E is photon energy in eV or keV (default is LCLS value)
    """
    E = eV(E)
    d = dSpace(ID,hkl)
    theta = asind(lam(E)/2/d)
    return theta


def BraggEnergy(ID,hkl,twotheta):
    """ Computes the photon energy that satisfies the Bragg condition of the specified material, reflection and twotheta angle 
        ID is chemical fomula : 'Si'
        hkl is the reflection : (1,1,1)
        twotheta is the scattering angle in degrees
    """
    ID=goodID(ID)
    d=dSpace(ID,hkl)
    l=2*d*sind(twotheta/2.0)  # noqa: E741
    E=lam2E(l)
    return E


def StructureFactor(ID,f,hkl,z=None):
    """ Computes the structure factor
        ID is chemical fomula : 'Si'
        f is the atomic form factor
        hkl is the reflection : (1,1,1)
        z is the rhombohedral lattice parameter
    """
    ID=goodID(ID, item='crystal')
    i=complex(0,1)
    h=hkl[0]
    k=hkl[1]
    l=hkl[2]  # noqa: E741
    L=latticeType[ID]
    if L=='fcc':
      F=f*(1+np.exp(-i*np.pi*(k+l))+np.exp(-i*np.pi*(h+l))+np.exp(-i*np.pi*(h+k)))
    elif L=='bcc':
      F=f*(1+np.exp(-i*np.pi*(h+k+l)))  
    elif L=='cubic':
        F=f
    elif L=='diamond':
        F=f*(1+np.exp(-i*np.pi*(k+l))+np.exp(-i*np.pi*(h+l))+np.exp(-i*np.pi*(h+k)))*(1+np.exp(-i*2*np.pi*(h/4.0+k/4.0+l/4.0)))
    # elif L=='rhomb':
    #     z=latticeParamRhomb[ID]
    #     F=f*(1+np.exp(2*i*np.pi*(h+k+l)*z)) 
    elif L=='tetr':
        F=f
    elif L=='hcp':
        F=f*(1+np.exp(2*i*np.pi*(h/3.0+2*k/3.0+l/2.0)))
    else:
        raise Exception(f'Unrecognized L: {L}')
    return F


def StructureFactorE(ID,hkl,E=None,z=None):
    ID=goodID(ID, item='crystal')
    # E = getE(E)
    theta=BraggAngle(ID,hkl,E)
    f=FF(ID,2*theta,E)
    return StructureFactor(ID,f,hkl,z)

    
def UnitCellVolume(ID):
    """ Returns the unit cell volume in m^3
        ID is chemical fomula : 'Si'
    """   
    ID=goodID(ID)
    lp=latticeParameters[ID]
    a=lp[0]/u['ang']
    b=lp[1]/u['ang']
    c=lp[2]/u['ang']
    alpha=lp[3]
    beta=lp[4]
    gamma=lp[5]
    # L=latticeType[ID]
    ca=cosd(alpha)
    cb=cosd(beta)
    cg=cosd(gamma)
    V=a*b*c*np.sqrt(1-ca**2-cb**2-cg**2+2*ca*cb*cg)
    return V

# def DebyeWallerFactor(ID,hkl,T=293,E=None):
#     """ Computes the Debye Waller factor for a specified reflection
#         ID is chemical fomula : 'Si'
#         T is the crystal temperature in Kelvin (default is 293)
#         E is photon energy in eV or keV (default is LCLS value)
#     """
#     ID=goodID(ID)
#     # E = getE(E)
#     theta=BraggAngle(ID,hkl,E)
#     l=lam(E)*u['ang']
#     T_Debye=debyeTemp[ID]
#     mass=MolecularMass(ID)
#     y=lambda x: x/(np.exp(x)-1)
#     ratio=T_Debye/float(T)
#     intgrl,err=s.integrate.quad(y,1e-9,ratio)
#     phi=intgrl*T/T_Debye    
#     B=11492*T*phi/(mass*T_Debye**2)+2873/(mass*T_Debye)
#     M=B*(sind(theta)/l)**2
#     DWfactor=np.exp(-M)
#     return DWfactor


def DarwinWidth(ID, hkl, E, T=293):
    """ Computes the Darwin width for a specified crystal reflection (degrees)
        ID is chemical fomula : 'Si'
        hkl is the reflection : (1,1,1)
        E is photon energy in eV or keV
        T is the crystal temperature in Kelvin (default is 293)
    """
    ID = goodID(ID, item='crystal')
    E = eV(E)
    theta = BraggAngle(ID,hkl,E)
    l = lam(E)  # noqa: E741
    f = FF(ID,2*theta,E)
    F = StructureFactor(ID,f,hkl)
    V = UnitCellVolume(ID)
    dw=(2*c['eRad']*l**2*np.abs(F))/(np.pi*V*sind(2*theta))/u['rad']
    return dw


def DarwinWidthE(ID, hkl, E, T=293):
    """ Computes the delta E corresponding to the Darwin width for a specified crystal reflection (degrees)
        ID is chemical fomula : 'Si'
        hkl is the reflection : (1,1,1)
        E is photon energy in eV or keV
        T is the crystal temperature in Kelvin (default is 293)
    """
    E = eV(E)
    return DeltaEoE(ID, hkl, E, T) * E
    # ID = goodID(ID)
    # E = eV(E)
    # dw       = DarwinWidth(ID,hkl,E,T)
    # theta    = BraggAngle(ID,hkl,E)
    # DeltaEoE = 1/tand(theta) * dw * u['rad']
    # DeltaE   = DeltaEoE * E
    # return DeltaE


def DeltaEoE(ID, hkl, E, T=293):
    E  = eV(E)
    theta = BraggAngle(ID,hkl,E)
    dw = DarwinWidth(ID,hkl,E,T)
    return 1/tand(theta) * dw * u['rad']


def index(ID,E=None):
    ID=goodID(ID)
    E = keV(E)
    d=Density[ID]
    n_real=xl.Refractive_Index_Re(ID,E,d)
    n_imag=xl.Refractive_Index_Im(ID,E,d)
    n=complex(n_real,n_imag)
    return n


def FF(ID,twotheta,E=None):
    """
    Returns the atomic form factor for Rayleigh scattering
    ID is the element name
    twotheta is the scattering angle in degrees
    E is the photon energy (default us current LCLS value)
    """
    E = keV(E)
    ID=goodID(ID)
    z=elementZ[ID]
    q=MomentumTransfer(E, twotheta)
    f=xl.FF_Rayl(z,q)
    return f


def MomentumTransfer(E,twotheta):
  """Returns the momentum transfer in Ang^-1 from xraylib [sin(2theta)/lam]
     E is the photon energy (eV or KeV)
     twotheta is the scattering angle in degrees
  """
  E = keV(E)
  th = np.deg2rad(twotheta)
  p = xl.MomentTransf(E, th)
  return p


def hkl2str(hkl):
    ''' Transfer hkl index to string
    '''
    return f'({hkl[0]},{hkl[1]},{hkl[2]})'


def getX0h(matID, keV, hkl):
    """ Get the Bragg reflection for a given material, photon energy and Bragg plane.
        
        This function uses the X0h database from the X0h CGI interface of the X-ray server at the Advanced Photon Source (APS) at Argonne National Laboratory (ANL).
        https://x-server.gmca.aps.anl.gov/cgi/www_form.exe?template=x0h_form.htm

    Args:
        matID (String): Material ID. Available materials are listed in the function.
        keV (float): Photon energy in keV
        hkl (tuple, int): Bragg plane (h, k, l)

    Returns:
        Response object: returns the Response object from the request.
    """
    url = 'https://x-server.gmca.aps.anl.gov/cgi/x0h_form.exe'

    matID_list = ['ADP', 'AlAs', 'AlFe3', 'AlN', 'AlP', 'AlSb', 'Aluminium', 'Ba2ScRuO6', 'BaSnO3', 'BaTiO3', 'Beril', 'Beryllium', 'Bi12GeO20', 'BiFeO3', 'Bismuth-fcc', 'Bismuth-primitive', 'Black_Phosphorus', 'C8H5KO4', 'C9H10N2', 'Ca2RuO4_HiTemp', 'Ca2RuO4_LoTemp', 'CaCO3_R3c', 'CaMnO3', 'CaRuO3', 'CdS', 'CdSe', 'CdTe', 'CeO2', 'Co2FeSi', 'Co2TiSi', 'CoAs', 'Copper', 'CoPS3_250K', 'Cr2AlC', 'CsCl', 'CsDSO4', 'CsF', 'CsH2PO4', 'Cu2ZnSnSe4', 'Cu3Au', 'Diamond', 'DyScO3', 'Fe-alpha', 'Fe2As', 'Fe2Mo3O8', 'Fe3O4', 'Fe3Si', 'FeBO3', 'FeGe2-alpha', 'FePS3', 'FeRh', 'FeSi', 'FeTiO3', 'Forsterite', 'Ga2O3-alpha', 'Ga2O3-beta', 'Ga2O3-gamma', 'Ga2O3-kappa', 'GaAs', 'GaN', 'GaN_cubic', 'GaP', 'GaSb', 'Gd2O3', 'Gd3Ga5O12', 'Gd3Sc2Ga3O12', 'GdSb', 'GdScO3', 'GeFe_primitive', 'Germanium', 'Ge_primitive', 'Gold', 'Graphite', 'HfO2', 'HgS', 'HgSe', 'HgTe', 'InAs', 'InGaAs', 'InN', 'InP', 'InSb', 'Iron-bcc', 'KAP', 'KCl', 'KDP', 'KTiOPO4', 'La(.5)Sr(1.5)MnO4', 'La(.7)Sr(.3)MnO3', 'La2CuO4_tetragonal', 'La2O3cub', 'La2O3hex', 'LaAlO3', 'LaCuO3', 'LaFeO3', 'LaInO3', 'LaMnO3', 'Langasite', 'LaNiO2', 'LaNiO3', 'LaNiO3_cubic', 'Lead', 'Li0.5Fe2.5O4_solutn', 'Li0.5Fe2.5O4_stchmtr', 'Li2B4O7', 'LiF', 'LiH', 'LiNbO3', 'LiTaO3', 'LSAT_cubic', 'LSAT_tetragonal', 'Lu2O3cub', 'LuPtBi', 'Magnetite', 'MgAl2O4', 'MgO', 'Mica', 'MnAs', 'Molybdenite', 'NaCl', 'NaOsO3', 'NbO2', 'NdGaO3', 'NdScO3', 'Nickel', 'Paratellurite', 'PbMg.24Nb.47Ti.29O3', 'PbMg.24Nb.48Ti.28O3r', 'PbMg.25Nb.49Ti.26O3h', 'PbMoO4', 'PbSe', 'PbTe', 'PbTiO3', 'PbWO4', 'PbZrO3', 'Pentaerythritol', 'Platinum', 'Pr2O3', 'PZT_PbZr.52Ti.48O3', 'Quartz', 'RuO2', 'Rutile', 'Sapphire_hex', 'Sapphire_rhomb', 'Sc2O3', 'Si-primitive', 'SiC-3c', 'SiC-4H', 'SiC-6H', 'SiFe-primitive', 'Silicon', 'Silver', 'SmNiO3', 'Sn2P2Se6', 'SnO2_Cassiterite', 'Sr2IrO4', 'Sr2TiO4', 'Sr3Al2O6', 'Sr3Ti2O7', 'Sr4Ti3O10', 'SrRuO3', 'SrTiO3', 'SrTiO3_tetragonal', 'SrVO3', 'Tantalum', 'Tellurium-I', 'Triglycine sulfate', 'Tungsten', 'Tungsten Disulfide', 'UO2', 'VO2_Rutile', 'Y3Al5O12', 'YAlO3', 'Zincite', 'Zn3P2', 'ZnS', 'ZnSe', 'ZnTe', 'ZrO2']
    matID_alias = {
        'Si': 'Silicon',
        'YAG': 'Y3Al5O12',
        'Ge': 'Germanium'
    }
    if matID in matID_alias:
        matID = matID_alias[matID]
    if matID not in matID_list:
        raise ValueError(f"Material ID {matID} not found in the list of available materials.")
    
    # Parameters of X0h CGI interface
    # Example from X0h search: https://x-server.gmca.aps.anl.gov/cgi/x0h_form.exe?#xway=2&wave=10&coway=0&code=Silicon&i1=2&i2=2&i3=0&df1df2=-1&modeout=1
    xway = 2    # 2 - keV
    wave = keV
    coway = 0
    code = matID
    i1, i2, i3 = hkl[0], hkl[1], hkl[2]
    df1df2 = -1  # automatic database selection
    modeout = 1  # text output (no table)

    # Building address
    params = {
        'xway': xway,
        'wave': wave,
        'coway': coway,
        'code': code,
        'i1': i1,
        'i2': i2,
        'i3': i3,
        'df1df2': df1df2,
        'modeout': modeout
    }
    
    # Connect & Download
    response = requests.get(url, params=params)

    # Example to read the data
    # content = response.text
    # FWHMurad_sg: Darwin width (FWHM) in micro-radians for sigma polarization
    # DarwinWidth = float(content.split('FWHMurad_sg= ')[1].strip().split(' ')[0]) * 1e-6
    # prcnsg: relative diffraction intensity for sigma polarization
    # xh_x0 = float(content.split('prcnsg= ')[1].strip().split(' ')[0]) * 0.01

    return response


def getX0h_plane(matID, keV, hkl_min, hkl_max, qb1=0., qb2=90., prcmin=0.1, df1df2=-1, hkl_base=(1,1,1), q1=0., q2=180.):
    """ Get the relative diffraction intensity of a Bragg reflection for a given material and photon energy.
        This function scans hkl planes in the range of hkl_min to hkl_max.
        "prcmin" must be > 0 to generate results for the relative diffraction intensity.
        
        This function uses the X0h database from the X0h CGI interface of the X-ray server at the Advanced Photon Source (APS) at Argonne National Laboratory (ANL).
        https://x-server.gmca.aps.anl.gov/cgi/www_form.exe?template=x0p_form.htm

    Args:
        matID (String): Material ID. Available materials are listed in the function.
        keV (float): Photon energy in keV
        hkl_min (tuple, int): Bragg plane range: minimum.
        hkl_max (tuple, int): Bragg plane range: maximum.
        qb1 (float, optional): Bragg angle range: minimum. Defaults to 0.
        qb2 (float, optional): Bragg angle range: maximum. Defaults to 90.
        prcmin (float, optional): Report threashold in percent. Must be > 0. Report when xh/x0 is greater than the given value. Defaults to 0.1
        df1df2 (int, optional): Database for dispersion corrections. Defaults to -1 (auto select).
        hkl_base (tuple, int, optional): Bragg plane of surface. Defaults to (1,0,0).
        q1 (float, optional): Plans make angle Theta1. Defaults to 0.
        q2 (float, optional): Plans make angle Theta2. Defaults to 180.

    Returns:
        Response object: returns the Response object from the request.
        Data have 4 columns: 'hkl' in the format of (h k l), 'Angle to surface', 'Bragg angle' and 'Relative Intensity xh/x0(%)'
    """
    url = 'https://x-server.gmca.aps.anl.gov/cgi/x0p_form.exe'
    
    matID_list = ['ADP', 'AlAs', 'AlFe3', 'AlN', 'AlP', 'AlSb', 'Aluminium', 'Ba2ScRuO6', 'BaSnO3', 'BaTiO3', 'Beril', 'Beryllium', 'Bi12GeO20', 'BiFeO3', 'Bismuth-fcc', 'Bismuth-primitive', 'Black_Phosphorus', 'C8H5KO4', 'C9H10N2', 'Ca2RuO4_HiTemp', 'Ca2RuO4_LoTemp', 'CaCO3_R3c', 'CaMnO3', 'CaRuO3', 'CdS', 'CdSe', 'CdTe', 'CeO2', 'Co2FeSi', 'Co2TiSi', 'CoAs', 'Copper', 'CoPS3_250K', 'Cr2AlC', 'CsCl', 'CsDSO4', 'CsF', 'CsH2PO4', 'Cu2ZnSnSe4', 'Cu3Au', 'Diamond', 'DyScO3', 'Fe-alpha', 'Fe2As', 'Fe2Mo3O8', 'Fe3O4', 'Fe3Si', 'FeBO3', 'FeGe2-alpha', 'FePS3', 'FeRh', 'FeSi', 'FeTiO3', 'Forsterite', 'Ga2O3-alpha', 'Ga2O3-beta', 'Ga2O3-gamma', 'Ga2O3-kappa', 'GaAs', 'GaN', 'GaN_cubic', 'GaP', 'GaSb', 'Gd2O3', 'Gd3Ga5O12', 'Gd3Sc2Ga3O12', 'GdSb', 'GdScO3', 'GeFe_primitive', 'Germanium', 'Ge_primitive', 'Gold', 'Graphite', 'HfO2', 'HgS', 'HgSe', 'HgTe', 'InAs', 'InGaAs', 'InN', 'InP', 'InSb', 'Iron-bcc', 'KAP', 'KCl', 'KDP', 'KTiOPO4', 'La(.5)Sr(1.5)MnO4', 'La(.7)Sr(.3)MnO3', 'La2CuO4_tetragonal', 'La2O3cub', 'La2O3hex', 'LaAlO3', 'LaCuO3', 'LaFeO3', 'LaInO3', 'LaMnO3', 'Langasite', 'LaNiO2', 'LaNiO3', 'LaNiO3_cubic', 'Lead', 'Li0.5Fe2.5O4_solutn', 'Li0.5Fe2.5O4_stchmtr', 'Li2B4O7', 'LiF', 'LiH', 'LiNbO3', 'LiTaO3', 'LSAT_cubic', 'LSAT_tetragonal', 'Lu2O3cub', 'LuPtBi', 'Magnetite', 'MgAl2O4', 'MgO', 'Mica', 'MnAs', 'Molybdenite', 'NaCl', 'NaOsO3', 'NbO2', 'NdGaO3', 'NdScO3', 'Nickel', 'Paratellurite', 'PbMg.24Nb.47Ti.29O3', 'PbMg.24Nb.48Ti.28O3r', 'PbMg.25Nb.49Ti.26O3h', 'PbMoO4', 'PbSe', 'PbTe', 'PbTiO3', 'PbWO4', 'PbZrO3', 'Pentaerythritol', 'Platinum', 'Pr2O3', 'PZT_PbZr.52Ti.48O3', 'Quartz', 'RuO2', 'Rutile', 'Sapphire_hex', 'Sapphire_rhomb', 'Sc2O3', 'Si-primitive', 'SiC-3c', 'SiC-4H', 'SiC-6H', 'SiFe-primitive', 'Silicon', 'Silver', 'SmNiO3', 'Sn2P2Se6', 'SnO2_Cassiterite', 'Sr2IrO4', 'Sr2TiO4', 'Sr3Al2O6', 'Sr3Ti2O7', 'Sr4Ti3O10', 'SrRuO3', 'SrTiO3', 'SrTiO3_tetragonal', 'SrVO3', 'Tantalum', 'Tellurium-I', 'Triglycine sulfate', 'Tungsten', 'Tungsten Disulfide', 'UO2', 'VO2_Rutile', 'Y3Al5O12', 'YAlO3', 'Zincite', 'Zn3P2', 'ZnS', 'ZnSe', 'ZnTe', 'ZrO2']
    matID_alias = {
        'Si': 'Silicon',
        'YAG': 'Y3Al5O12',
        'Ge': 'Germanium'
    }
    if matID in matID_alias:
        matID = matID_alias[matID]
    if matID not in matID_list:
        raise ValueError(f"Material ID {matID} not found in the list of available materials.")

    # Parameters of X0h CGI interface
    # Example from X0h search: https://x-server.gmca.aps.anl.gov/cgi/x0p_form.exe?xway=2&wave=11&code=Y3Al5O12&hkl11=0&hkl12=0&hkl13=0&hkl21=5&hkl22=5&hkl23=5&qb1=0.&qb2=90.&prcmin=1&df1df2=-1&base1=1&base2=0&base3=0&modesearch=3&q1=0.&q2=180.
    xway = 2    # 2 - keV
    wave = keV
    code = matID
    hkl11, hkl12, hkl13 = hkl_min[0], hkl_min[1], hkl_min[2]
    hkl21, hkl22, hkl23 = hkl_max[0], hkl_max[1], hkl_max[2]
    base1, base2, base3 = hkl_base[0], hkl_base[1], hkl_base[2]
    modesearch = 3

    # Building address
    params = {
        'xway': xway,
        'wave': wave,
        'code': code,
        'hkl11': hkl11,
        'hkl12': hkl12,
        'hkl13': hkl13,
        'hkl21': hkl21,
        'hkl22': hkl22,
        'hkl23': hkl23,
        'qb1': qb1,
        'qb2': qb2,
        'prcmin': prcmin,
        'df1df2': df1df2,
        'base1': base1,
        'base2': base2,
        'base3': base3,
        'modesearch': modesearch,
        'q1': q1,
        'q2': q2
    }

    # Connect & Download
    response = requests.get(url, params=params)

    # Example to read the data
    # import pandas as pd
    # tables = pd.read_html(response.text)
    # _df = tables[-1]  # Return 4 columns: 'hkl' in the format of string (h k l), 'Angle to surface', 'Bragg angle' and 'Relative Intensity xh/x0(%)'
    
    # Pick data
    # data = df[(df['hkl'] == '(0 0 1)') & (df['keV'] == 10.0)]
    # bragg = df['Bragg angle'].to_numpy()

    # Convert string to tuple
    # hkl = df['hkl'].str.strip('()')
    # hkl = [tuple(map(int, x.split())) for x in hkl]
    
    return response


def loglog_negy_interp1d(x, y, xx):
    ''' 1D log-log interpolation supprting negative y values:
            linear interpolation for data before the last negative data;
            log-log for data afterwards
    '''
    if np.where(y<=0)[0].shape[0] > 0:
        idx = np.where(y<=0)[0][-1] + 1
        x1 = x[:idx+1]
        y1 = y[:idx+1]
        x2 = x[idx:]
        y2 = y[idx:]
        if xx[-1] <= x1[-1]:    # all data in linear interpolation region
            return interp1d(x1, y1)(xx)
        elif xx[0] >= x2[0]:    # all data in log-log interpolation region
            return np.exp(interp1d(np.log(x2), np.log(y2))(np.log(xx)))
        else:
            idxx = np.where(xx<x1[-1])[0][-1]
            xx1 = xx[:idxx+1]
            xx2 = xx[idxx+1:]
            yy1 = interp1d(x1, y1)(xx1)
            yy2 = np.exp(interp1d(np.log(x2), np.log(y2))(np.log(xx2)))
            return np.concatenate((yy1, yy2))
    else:   # all data are positive
        return np.exp(interp1d(np.log(x), np.log(y), bounds_error=False, fill_value='extrapolate')(np.log(xx)))

    
def get_Ef1f2(Z, datafile='f1f2_EPDL97.dat'):
    Ef1f2 = np.array([], dtype=np.float64)
    path = os.path.dirname(os.path.abspath(__file__))
    with open(path+'/'+datafile, 'r') as inp:
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


def Reflectivity_uncoated(E, theta, sub_mat, f1f2data='default', f1f2interp='linear'):
    '''
        E: photon energy in eV
        theta: incideng angle in rad
        sub_mat: material ID of substrate
        f1f2data:   'default', use 'f1f2_EPDL97.dat' file;
                    name data file (currently have 'f1f2_Windt.dat')
        f1f2interp: 'linear', default, linear interpolation for photon energy
                    else, log-log interpolation
    '''
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
    if f1f2interp == 'linear':
        f1 = interp1d(f1f2[:,0], f1f2[:,1])(E)
        f2 = interp1d(f1f2[:,0], f1f2[:,2])(E)
    else:
        f1 = loglog_negy_interp1d(f1f2[:,0], f1f2[:,1], E)
        f2 = loglog_negy_interp1d(f1f2[:,0], f1f2[:,2], E)

    lamda = c['h'] * c['c'] / (E * c['e']) + 0j
    K = 2 * np.pi / lamda
    ah = c['eRad'] * lamda**2 / np.pi
    rho = Density[sub_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho = rho / AtomicMass[sub_mat] * c['NA']   # Atoms per m3

    Chi = ah * n_rho * (-f1 + f2*1j)
    K_z0 = K * np.sqrt(np.sin(theta)**2 + Chi)
    K_z1 = K * np.sin(theta)
    C1 = np.exp(1j * K_z1 * 1 / 2)
    R_1 = (K_z1 - K_z0) / (K_z1 + K_z0) * C1**2
    R = np.abs(R_1) **2
    
    if scalar_E:
        return np.squeeze(R)
    return R


def Reflectivity_coated(E, theta, sub_mat, coat_mat, coat_thickness, f1f2data='default', f1f2interp='linear'):
    '''
        E: photon energy in eV
        theta: incideng angle in rad
        sub_mat: material ID of substrate
        coat_mat: material ID of coating
        coat_thickness: coating thickness in m
        f1f2data:   'default', use 'f1f2_EPDL97.dat' file;
                    name data file (currently have 'f1f2_Windt.dat')
        f1f2interp: 'linear', default, linear interpolation for photon energy
                    else, log-log interpolation
    '''
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
    if f1f2interp == 'linear':
        f1_sub = interp1d(f1f2_sub[:,0], f1f2_sub[:,1])(E)
        f2_sub = interp1d(f1f2_sub[:,0], f1f2_sub[:,2])(E)
    else:
        f1_sub = loglog_negy_interp1d(f1f2_sub[:,0], f1f2_sub[:,1], E)
        f2_sub = loglog_negy_interp1d(f1f2_sub[:,0], f1f2_sub[:,2], E)
    rho_sub = Density[sub_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho_sub = rho_sub / AtomicMass[sub_mat] * c['NA']   # Atoms per m3

    coat = xl.CompoundParser(coat_mat)
    f1_coat = np.zeros_like(E)
    f2_coat = np.zeros_like(E)
    for i in range(coat['nElements']):
        if f1f2data == 'default':
            f1f2_coat = get_Ef1f2(coat['Elements'][i])
        else:
            f1f2_coat = get_Ef1f2(coat['Elements'][i], datafile=f1f2data)
        if f1f2interp == 'linear':
            f1_coat += interp1d(f1f2_coat[:,0], f1f2_coat[:,1])(E) * coat['nAtoms'][i]
            f2_coat += interp1d(f1f2_coat[:,0], f1f2_coat[:,2])(E) * coat['nAtoms'][i]
        else:
            f1_coat += loglog_negy_interp1d(f1f2_coat[:,0], f1f2_coat[:,1], E) * coat['nAtoms'][i]
            f2_coat += loglog_negy_interp1d(f1f2_coat[:,0], f1f2_coat[:,2], E) * coat['nAtoms'][i]

    rho_coat = Density[coat_mat] * 1.e6  # g/cm3 --> g/m3
    n_rho_coat = rho_coat / coat['molarMass'] * c['NA']   # Atoms per m3
    
    lamda = c['h'] * c['c'] / (E * c['e']) + 0j
    K = 2 * np.pi / lamda
    ah = c['eRad'] * lamda**2 / np.pi
    
    Chi_sub = ah * n_rho_sub * (-f1_sub + f2_sub*1j)
    Chi_coat = ah * n_rho_coat * (-f1_coat + f2_coat*1j)
    K_z0 = K * np.sqrt(np.sin(theta)**2 + Chi_sub)
    K_z1 = K * np.sqrt(np.sin(theta)**2 + Chi_coat)
    K_z2 = K * np.sin(theta)
    C1 = np.exp(1j * K_z1 * coat_thickness / 2)
    C2 = np.exp(1j * K_z2 * 1 / 2)
    R1 = (K_z1 - K_z0) / (K_z0 + K_z1) * C1**2   # R1 = r11_0
    r11_2 = (K_z1 - K_z2) / (K_z2 + K_z1) * C1**2
    r22_1 = (K_z2 - K_z1) / (K_z1 + K_z2) * C2**2
    t12 = 2 * K_z2 / (K_z1 + K_z2) * C1 * C2
    t21 = 2 * K_z1 / (K_z2 + K_z1) * C1 * C2
    R2 = r22_1 + t21 * t12 * R1 / (1. - r11_2 * R1)
    R = np.abs(R2)**2
    
    if scalar_E:
        return np.squeeze(R)
    return R


def mirror_defocus_ratio(src, mirror, focus, target):
    div_ratio = (mirror - src) / (focus - mirror)
    return np.abs(target - focus) / (target - src) * div_ratio


def center_temp(power_W, rms_mm, attL_m, K):
    ''' Surface temperature at beam center when beam hits a large material
       power_W: beam power in Watt
       rms_mm:  beam rms size in mm
       attL_m:  attenuation length in meter
       K:       thermal conductivity in W/m/K
    '''
    W = 1/attL_m * rms_mm * 0.001
    Tmax = power_W / (2 * K * rms_mm*0.001 * np.sqrt(np.pi))
    if np.isscalar(W):
        if W < 50:
            N = W * dawsn(W/2) - W/(2 * np.sqrt(np.pi)) * np.exp(-W**2/4) * expi(W**2/4)
        else:
            N = 1 - 2 / np.sqrt(np.pi) / W
    else:
        W1 = W[W<50]
        W2 = W[W>=50]
        N = np.concatenate((W1 * dawsn(W1/2) - W1/(2 * np.sqrt(np.pi)) * np.exp(-W1**2/4) * expi(W1**2/4),
                            1 - 2 / np.sqrt(np.pi) / W2))
    return Tmax * N

def keV2lamda(keV):
    ''' Covert photon energy to wavelength in m (SI unit)
    '''
    return 1239.8 / (keV*1000) * 1e-9

def lamda2keV(lamda_nm):
    ''' Convert photon wavelength to energy in keV
        lamda_nm: wavelength in nm
    '''
    return 1239.8 / lamda_nm * 0.001


## Gaussian beam propogation
def RayleighLength(keV, sigma0_um, n=1):
    ''' Return the Rayleigh range in m (SI unit)
       keV: photon energy in keV
       sigma0_um: sigma in um
       n:   refractive index of the medium (default 1 for vacuum). Not tested for non-vacuum yet.
    '''
    return 4 * np.pi * (sigma0_um*1e-6)**2 * n / keV2lamda(keV)

def LensFocalLength(keV, radius_um, matID='Be'):
    ''' Return the focal length of a lens in m
        Note: Use LensFocalDistance to calculate for a specific source
        keV: photon energy in keV
        radius_um: effectiv radius of lens
        matID: material of lens, default as Be (beryllium)
    '''
    delta = 1 - xl.Refractive_Index_Re(matID, keV, Density[matID])
    return radius_um * 1e-6 / (2 * delta)

def LensFocalDistance(s, f, LR):
    ''' Return the distance from the focal point to lens
        Note: no focusing if s > f
        s: Distance from source to lens
        f: Lens focal length (refer LensFocalLength function)
        LR: Raileigh Lengh (refer RayleighRange function)
    '''
    if f < s:
        return f + (s - f) / ((s/f - 1)**2 + (LR/f)**2)
    else:
        raise Warning(f'Source distance ({s}) > focal length ({f}): Beam will NOT be focused.')
        return np.Inf

def LensFocalSigma(s, f, LR, sigma0):
    ''' Return the focal size
        Note: no focusing if s > f
        s: Distance from source to lens
        f: Lens focal length (refer LensFocalLength function)
        LR: Raileigh Lengh (refer RayleighRange function)
        sigma0: source sigma; half of waist (sigma0 = w0/2)
    '''
    if f < s:
        return sigma0 / np.sqrt((1 - s/f)**2 + (LR/f)**2)
    else:
        raise Exception(f'Source distance ({s}) > focal length ({f}): Beam will NOT be focused.')


def initGaussianBeam(keV, sigma0_um, M_square=1.0):
    """_summary_

    Args:
        keV (float, array_like): keV
        sigma0_um (float, array_like): sigma in um
        M_square(float, array_like): M^2. Default to 1.
    Returns:
        complex, array_like: q of Gaussian beam
    """    
    zR = RayleighLength(keV, sigma0_um)
    return zR / M_square * 1j

def propagateFreeSpace(q_in, d_m):
    """ free space propogation

    Args:
        q_in (complex, array_like): input Gaussian Beam, SI unit
        d_m (float, array_like): free space lengh in meters

    Returns:
        q_out (complex, array_like): q out
    """
    return q_in + d_m

def propagateThinLens(q_in, f_m):
    """ thin lens propogation

    Args:
        q_in (complex, array_like): input Gaussian Beam, SI unit
        f_m (float, array_like): lens focal length in meters

    Returns:
        q_out (complex, array_like): q out
    """
    return 1 / (1/q_in - 1/f_m)


def getGaussianSigma_um(q, keV, M_square=1.0):
    """ get sigma in um of q

    Args:
        q (complex, array_like): complex Gaussian beam
        keV (float, array_like): keV
        M_square (float, array_like): M^2, beam quality factor. Default to 1.

    Returns:
        sigma_um (float, array_like): sigma size in um
    """
    im = (1/q).imag
    lamda = keV2lamda(keV)
    return np.sqrt(-M_square * lamda / (4 * np.pi * im)) * 1e6


class GaussianBeam:
    def __init__(self, keV, sigma0_um, M_square=1.0):
        ''' Initialize a Gaussian beam
            keV: photon energy
            z0:  waist location in m
            sigma0:  sigma size in um
            M_square: beam quality factor, default to 1
        '''
        self.keV = keV
        self.sigma = sigma0_um
        self.M2 = M_square
        self.w = self.sigma * 2
        zR = RayleighLength(keV, sigma0_um)
        self.lamda = keV2lamda(keV)
        self.q0 = zR / M_square * 1j
        self.q = self.q0
        self.propogation = []

    ## Transport functions
    def addFreeSpace(self, d):
        ''' Transport in free space (n=1)
            d: transport distance in m
        '''
        self.propogation.append(("Free Space", d))
        self.q = self.q + d
        self.R = self.getR_m()
        self.w = self.getWaist_um()
        self.sigma = self.w / 2

    def addThinLens(self, f):
        ''' Transport through a think lens
            f: focal length in m
        '''
        self.propogation.append(("Thin Lens", f))
        self.q = 1 / (1/self.q - 1/f)
        self.R = self.getR_m()
        self.w = self.getWaist_um()
        self.sigma = self.w / 2

    ## Get beam information
    def getR_m(self):
        ''' Get the radius of curvature in m (positive for diverging, negative for convergin)
        '''
        # if self.q.real == 0:
        #     return np.inf
        # else:
        return 1 / (1/self.q).real

    def getWaist_um(self):
        ''' Get the waist in um
        '''
        im = (1/self.q).imag
        return np.sqrt(-self.M2 * self.lamda / (np.pi * im)) * 1e6

    def getSigma_um(self):
        ''' Get the sigma of beam: sigma = waist / 2
        '''
        return self.getWaist_um() / 2

    def print(self):
        ''' Print beam information
        '''
        print(f'{self.keV} keV, {self.sigma0} um waist:')
        print('Propogation List:')
        for i in range(len(self.propogation)):
            if self.propogation[i][0] == "Free Space":
                print(f'  {i+1}. Free Space {self.propogation[i][1]} m')
            if self.propogation[i][0] == "Thin Lens":
                print(f'  {i+1}. Thin Lens with {self.propogation[i][1]} m focal length')
        print(f'Complex beam parameter q = {self.q}')
        print(f'Radius of curvature R = {self.R} m')
        print(f'Beam waist = {self.w} um, rms = {self.sigma} um')

    def reset(nstep=0):
        ''' Reset '''
        pass

""" Define units and constants
    - Muliply to convert from SI unit to the target unit:
        e.g. a*u['cm'] to get from m to cm
    - Devide to covert from the indicated unit to SI unit
        e.g. a/u['cm'] to get from cm to m
"""
u = {'fm': 1e15,
     'pm': 1e12,
     'ang': 1e10,
     'nm': 1e9,
     'um': 1e6,
     'mm': 1e3,
     'cm': 1e2,
     'km': 1e-3,
     'kHz':1e-3,
     'MHz':1e-6,
     'GHz':1e-9,
     'THz':1e-12,
     'PHz':1e-15,
     'inch':39.370079,
     'mile':0.000621,
     'ft':3.28084,
     'yard':1.093613,
     'mil':39.370079*1000,
     'barn':1e28,  
     'fs':1e15,
     'ps':1e12,
     'ns':1e9,
     'us':1e6,
     'ms':1e3,
     'min':1/60.,
     'hour':1/3600.,
     'day':1/(3600*24.),
     'mdeg':1e3,
     'udeg':1e6,
     'ndeg':1e9,
     'rad':np.pi/180,
     'mrad':np.pi/180*1e3,
     'urad':np.pi/180*1e6,
     'nrad':np.pi/180*1e9,
     'asec':3600,
     'amin':60,
     'g':1e3,
     'eV':6.241509e+18,
     'erg':1e7,
     'cal':0.239,
     'mJ':1e3,
     'uJ':1e6,
     'nJ':1e9,
     'pJ':1e9,
     'Torr':7.5006e-3
}


# Constants in MKS units
#    NA : Avagadro Constant
#  eRad : Classical electron radius
#     e : electron charge in Coulombs
#     c : speed of light in m
#     h : Plank constant in Js
#  hbar : reduced Plank constant in Js
# emass : electron mass in kg
# pmass : proton mass in kg
#     k : Boltzman constant in J/K
#    mu : Permeability of free space in N/A^2
#    eo : Permittivity of free space in F/m
#   S-B : Stefan-Boltzmann constant
 
c = {'NA': 6.022140857e23,
     'eRad':2.81794e-15,  
     'e':1.6021766e-19,
     'c':2.99792458e8,  
     'h':6.626070e-34,
     'hbar':1.054571628e-34,
     'emass':9.10938215e-31, 
     'pmass':1.672621637e-27, 
     'k':1.3806504e-23,
     'mu':12.566370614e-7,  
     'eo':8.854187817e-12,
     'S-B':5.670373e-8,
}

# Element Name
elementName = {'H' :'Hydrogen',
     'He':'Helium',
     'Li':'Lithium',
     'Be':'Beryllium',
     'B':'Boron',
     'C':'Carbon',
     'N':'Nitrogen',
     'O':'Oxygen',
     'F':'Fluorine',
     'Ne':'Neon',
     'Na':'Sodium',
     'Mg':'Magnesium',
     'Al':'Aluminium',
     'Si':'Silicon',
     'P':'Phosphorus',
     'S':'Sulfur',
     'Cl':'Chlorine',
     'Ar':'Argon',
     'K':'Potassium',
     'Ca':'Calcium',
     'Sc':'Scandium',
     'Ti':'Titanium',
     'V':'Vanadium',
     'Cr':'Chromium',
     'Mn':'Manganese',
     'Fe':'Iron',
     'Co':'Cobalt',
     'Ni':'Nickel',
     'Cu':'Copper',
     'Zn':'Zinc',
     'Ga':'Gallium',
     'Ge':'Germanium',
     'As':'Arsenic',
     'Se':'Selenium',
     'Br':'Bromine',
     'Kr':'Krypton',
     'Rb':'Rubidium',
     'Sr':'Strontium',
     'Y':'Yttrium',
     'Zr':'Zirconium',
     'Nb':'Niobium',
     'Mo':'Molybdenum',
     'Tc':'Technetium',
     'Ru':'Ruthenium',
     'Rh':'Rhodium',
     'Pd':'Palladium',
     'Ag':'Silver',
     'Cd':'Cadmium',
     'In':'Indium',
     'Sn':'Tin',
     'Sb':'Antimony',
     'Te':'Tellurium',
     'I':'Iodine',
     'Xe':'Xenon',
     'Cs':'Caesium',
     'Ba':'Barium',
     'La':'Lanthanum',
     'Ce':'Cerium',
     'Pr':'Praseodymium',
     'Nd':'Neodymium',
     'Pm':'Promethium',
     'Sm':'Samarium',
     'Eu':'Europium',
     'Gd':'Gadolinium',
     'Tb':'Terbium',
     'Dy':'Dysprosium',
     'Ho':'Holmium',
     'Er':'Erbium',
     'Tm':'Thulium',
     'Yb':'Ytterbium',
     'Lu':'Lutetium',
     'Hf':'Hafnium',
     'Ta':'Tantalum',
     'W':'Tungsten',
     'Re':'Rhenium',
     'Os':'Osmium',
     'Ir':'Iridium',
     'Pt':'Platinum',
     'Au':'Gold',
     'Hg':'Mercury',
     'Tl':'Thallium',
     'Pb':'Lead',
     'Bi':'Bismuth',
     'Po':'Polonium',
     'At':'Astatine',
     'Rn':'Radon',
     'Fr':'Francium',
     'Ra':'Radium',
     'Ac':'Actinium',
     'Th':'Thorium',
     'Pa':'Protactinium',
     'U':'Uranium',
     'Np':'Neptunium',
     'Pu':'Plutonium'
}


# Chemical Formula Aliases
alias={'Air':'N1.562O.42C.0003Ar.0094',
       'air':'N1.562O.42C.0003Ar.0094',
       'C*':'C',
       'CVD':'C',
       'mylar':'C10H8O4',
       'Mylar':'C10H8O4',
       'polyimide':'C22H10N2O5',
       'Polyimide':'C22H10N2O5',
       'kapton':'C22H10N2O5',
       'Kapton':'C22H10N2O5',
       '304SS':'Fe.68Cr.2Ni.1Mn.02',
       'Acetone':'C3H6O',
       'acetone':'C3H6O',
       'PMMA':'C5H8O2',
       'Teflon':'C2F4',
       'teflon':'C2F4',
       'Toluene':'C7H8',
       'toluene':'C7H8',
       'FS':'SiO2',
       'GGG':'Gd3Ga5O12',
       'quartz':'SiO2',
       'Quartz':'SiO2',
       'Silica':'SiO2',
       'silica':'SiO2',
       'water':'H2O',
       'Water':'H2O',
       'Calcite':'CaCO3',
       'calcite':'CaCO3',
       'YAG':'Y3Al5O12',
       'yag':'Y3Al5O12',
       'Sapphire':'Al2O3',
       'sapphire':'Al2O3',
       'Blood':'CHN.3O7.6',
       'LMSO':'La0.7Sr0.3MnO3',
       'blood':'CHN.3O7.6',
       'Bone':'C1.5H0.3O4.3N0.4PCa2.2',
       'bone':'C1.5H0.3O4.3N0.4PCa2.2',
       'IF1':'Be0.9983O0.0003Al0.0001Ca0.0002C0.0003Cr0.000035Co0.000005Cu0.00005Fe0.0003Pb0.000005Mg0.00006Mn0.00003Mo0.00001Ni0.0002Si0.0001Ag0.000005Ti0.00001Zn0.0001',
       'PF60':'Be.994O.004Al.0005B.000003Cd.0000002Ca.0001C.0006Cr.0001Co.00001Cu.0001Fe.0008Pb.00002Li.000003Mg.00049Mn.0001Mo.00002Ni.0002N.0003Si.0004Ag.00001'
}



# Atomic Number
elementZ = {'H' :1,
     'He':2,
     'Li':3,
     'Be':4,
     'B':5,
     'C':6,
     'N':7,
     'O':8,
     'F':9,
     'Ne':10,
     'Na':11,
     'Mg':12,
     'Al':13,
     'Si':14,
     'P':15,
     'S':16,
     'Cl':17,
     'Ar':18,
     'K':19,
     'Ca':20,
     'Sc':21,
     'Ti':22,
     'V':23,
     'Cr':24,
     'Mn':25,
     'Fe':26,
     'Co':27,
     'Ni':28,
     'Cu':29,
     'Zn':30,
     'Ga':31,
     'Ge':32,
     'As':33,
     'Se':34,
     'Br':35,
     'Kr':36,
     'Rb':37,
     'Sr':38,
     'Y':39,
     'Zr':40,
     'Nb':41,
     'Mo':42,
     'Tc':43,
     'Ru':44,
     'Rh':45,
     'Pd':46,
     'Ag':47,
     'Cd':48,
     'In':49,
     'Sn':50,
     'Sb':51,
     'Te':52,
     'I':53,
     'Xe':54,
     'Cs':55,
     'Ba':56,
     'La':57,
     'Ce':58,
     'Pr':59,
     'Nd':60,
     'Pm':61,
     'Sm':62,
     'Eu':63,
     'Gd':64,
     'Tb':65,
     'Dy':66,
     'Ho':67,
     'Er':68,
     'Tm':69,
     'Yb':70,
     'Lu':71,
     'Hf':72,
     'Ta':73,
     'W':74,
     'Re':75,
     'Os':76,
     'Ir':77,
     'Pt':78,
     'Au':79,
     'Hg':80,
     'Tl':81,
     'Pb':82,
     'Bi':83,
     'Po':84,
     'At':85,
     'Rn':86,
     'Fr':87,
     'Ra':88,
     'Ac':89,
     'Th':90,
     'Pa':91,
     'U':92,
     'Np':93,
     'Pu':94
}


# Atomic Weight, unit are amu
AtomicMass = {'H' :1.00794,
     'He':4.002602,
     'Li':6.941,
     'Be':9.012182,
     'B':10.811,
     'C':12.0107,
     'N':14.0067,
     'O':15.9994,
     'F':18.9984032,
     'Ne':20.1797,
     'Na':22.98976928,
     'Mg':24.3050,
     'Al':26.9815386,
     'Si':28.0855,
     'P':30.973762,
     'S':32.065,
     'Cl':35.453,
     'Ar':39.948,
     'K':39.0983,
     'Ca':40.078,
     'Sc':44.955912,
     'Ti':47.867,
     'V':50.9415,
     'Cr':51.9961,
     'Mn':54.938045,
     'Fe':55.845,
     'Co':58.933195,
     'Ni':58.6934,
     'Cu':63.546,
     'Zn':65.38,
     'Ga':69.723,
     'Ge':72.64,
     'As':74.92160,
     'Se':78.96,
     'Br':79.904,
     'Kr':83.798,
     'Rb':85.4678,
     'Sr':87.62,
     'Y':88.90585,
     'Zr':91.224,
     'Nb':92.90638,
     'Mo':95.96,
     'Tc':98,
     'Ru':101.07,
     'Rh':102.90550,
     'Pd':106.42,
     'Ag':107.8682,
     'Cd':112.411,
     'In':114.818,
     'Sn':118.710,
     'Sb':121.760,
     'Te':127.60,
     'I':126.9044,
     'Xe':131.293,
     'Cs':132.9054519,
     'Ba':137.327,
     'La':138.90547,
     'Ce':140.116,
     'Pr':140.90765,
     'Nd':144.242,
     'Pm':145,
     'Sm':150.36,
     'Eu':151.964,
     'Gd':157.25,
     'Tb':158.92535,
     'Dy':162.500,
     'Ho':164.93032,
     'Er':167.259,
     'Tm':168.93421,
     'Yb':173.054,
     'Lu':174.9668,
     'Hf':178.49,
     'Ta':180.94788,
     'W':183.84,
     #'WC':195.85,
     'Re':186.207,
     'Os':190.23,
     'Ir':192.217,
     'Pt':195.084,
     'Au':196.966569,
     'Hg':200.59,
     'Tl':204.3833,
     'Pb':207.2,
     'Bi':208.9804,
     'Po':210,
     'At':210,
     'Rn':222,
     'Fr':223,
     'Ra':226,
     'Ac':227,
     'Th':232.03806,
     'Pa':231.03588,
     'U':238.02891,
     'Np':237,
     'Pu':244
}


# Material density in g/cm^3
Density = {'H' :0.00008988,
     'He':0.0001785,
     'Li':0.543,
     'Be':1.85,
     'B':2.34,
     'C':2.267,
     'N':0.0012506,
     'O':0.001429,
     'F':0.001696,
     'Ne':0.0008999,
     'Na':0.971,
     'Mg':1.738,
     'Al':2.698,
     'Si':2.329,
     'P':1.82,
     'S':2.067,
     'Cl':0.003214,
     'Ar':0.0017837,
     'K':0.862,
     'Ca':1.54,
     'Sc':2.989,
     'Ti':4.54,
     'V':6.11,
     'Cr':7.15,
     'Mn':7.44,
     'Fe':7.874,
     'Co':8.86,
     'Ni':8.912,
     'Cu':8.96,
     'Zn':7.134,
     'Ga':5.907,
     'Ge':5.323,
     'As':5.776,
     'Se':4.809,
     'Br':3.122,
     'Kr':0.003733,
     'Rb':1.532,
     'Sr':2.64,
     'Y':4.469,
     'Zr':6.506,
     'Nb':8.57,
     'Mo':10.22,
     'Tc':11.5,
     'Ru':12.37,
     'Rh':12.41,
     'Pd':12.02,
     'Ag':10.501,
     'Cd':8.69,
     'In':7.31,
     'Sn':7.287,
     'Sb':6.685,
     'Te':6.232,
     'I':4.93,
     'Xe':0.005887,
     'Cs':1.873,
     'Ba':3.594,
     'La':6.145,
     'Ce':6.77,
     'Pr':6.773,
     'Nd':7.007,
     'Pm':7.26,
     'Sm':7.52,
     'Eu':5.243,
     'Gd':7.895,
     'Tb':8.229,
     'Dy':8.55,
     'Ho':8.795,
     'Er':9.066,
     'Tm':9.321,
     'Yb':6.965,
     'Lu':9.84,
     'Hf':13.31,
     'Ta':16.654,
     'W':19.25,
     'WC':15.8,
     'Re':21.02,
     'Os':22.61,
     'Ir':22.56,
     'Pt':21.46,
     'Au':19.282,
     'Hg':13.5336,
     'Tl':11.85,
     'Pb':11.342,
     'Bi':9.807,
     'Po':9.32,
     'At':7,
     'Rn':0.00973,
     'Fr':1.87,
     'Ra':5.5,
     'Ac':10.07,
     'Th':11.72,
     'Pa':15.37,
     'U':18.95,
     'Np':20.45,
     'Pu':19.84,
     'CVD':3.515,
     'H2O':1.0,
     'B4C':2.52,
     'SiO2':2.2,
     'Al2O3':3.97,
     'ZnSe':5.42,
     'BeO':3.01,
     'SiC':3.21,
     'ZnTe':6.34,
     'CdS':6.749,
     'CdSe':7.01,
     'CdTe':7.47,
     'BN':3.49,
     'GaSb':5.619,
     'GaAs':5.316,
     'GaMnAs':5.316,
     'GaP':4.13,
     'InP':4.787,
     'InAs':5.66,
     'InSb':5.775,
     'TaC':13.9,
     'TiB2':4.52,
     'YAG':4.55,
     'CuBe':8.96,
     'ZnO':5.606,
     'SiC2':3.217,
     'AlN':3.3,
     'Si3N4':3.44,
     'CaF2':3.18,
     'LiF':2.635,
     'KF':2.48,
     'PbF2':8.24,
     'SrF2':4.24,
     'KBr':2.75,
     'ZrO2':5.6,
     'Gd3Ga5O12':7.08,
     'CaSiO5':2.4,
     'LaMnO3':5.7,
     'LaAlO3':6.52,
     'La0.7Sr0.3MnO3':6.17,
     'La0.5Ca0.5MnO3':6.3,
     'Fe.68Cr.2Ni.1Mn.02':8.03,
     'CaSO4H4O2':2.32,
     'C10H8O4':1.4,
     'C22H10N2O5':1.43,
     'C3H6O':0.79,
     'C5H8O2':1.19,
     'C2F4':2.2,
     'C7H8':0.867,
     'Y3Al5O12':4.56,
     'CHN.3O7.6':1.06,
     'C1.5H0.3O4.3N0.4PCa2.2':1.92,
     'Be0.9983O0.0003Al0.0001Ca0.0002C0.0003Cr0.000035Co0.000005Cu0.00005Fe0.0003Pb0.000005Mg0.00006Mn0.00003Mo0.00001Ni0.0002Si0.0001Ag0.000005Ti0.00001Zn0.0001':1.85,
     'Be.994O.004Al.0005B.000003Cd.0000002Ca.0001C.0006Cr.0001Co.00001Cu.0001Fe.0008Pb.00002Li.000003Mg.00049Mn.0001Mo.00002Ni.0002N.0003Si.0004Ag.00001':1.85,
     'Air':0.001225,  #  International Standard Atmosphere (ISA) values: 15 degC at sea level
     'air':0.001225
}

 
# Melting point in K
meltPoint = {'H' :14.175,
     'He':None,
     'Li':453.85,
     'Be':1560,
     'B':2573,
     'C':3948,
     'CVD':3948,
     'N':63.29,
     'O':50.5,
     'F':53.65,
     'Ne':24.703,
     'Na':371,
     'Mg':923,
     'Al':933.4,
     'Si':1683,
     'P':317.25,
     'S':388.51,
     'Cl':172.25,
     'Ar':83.96,
     'K':336.5,
     'Ca':1112,
     'Sc':1812,
     'Ti':1933,
     'V':2175,
     'Cr':2130,
     'Mn':1519,
     'Fe':1808,
     'Co':1768,
     'Ni':1726,
     'Cu':1357.75,
     'Zn':692.88,
     'Ga':302.91,
     'Ge':1211.45,
     'As':1090,
     'Se':494,
     'Br':266.05,
     'Kr':115.93,
     'Rb':312.79,
     'Sr':1042,
     'Y':1799,
     'Zr':2125,
     'Nb':2741,
     'Mo':2890,
     'Tc':2473,
     'Ru':2523,
     'Rh':2239,
     'Pd':1825,
     'Ag':1234,
     'Cd':594.33,
     'In':429.91,
     'Sn':505.21,
     'Sb':904.05,
     'Te':722.8,
     'I':386.65,
     'Xe':161.45,
     'Cs':301.7,
     'Ba':1002,
     'La':1193,
     'Ce':1071,
     'Pr':1204,
     'Nd':1289,
     'Pm':1204,
     'Sm':1345,
     'Eu':1095,
     'Gd':1585,
     'Tb':1630,
     'Dy':1680,
     'Ho':1743,
     'Er':1795,
     'Tm':1818,
     'Yb':1097,
     'Lu':1936,
     'Hf':2500,
     'Ta':3269,
     'W':3680,
     'WC':3168,
     'Re':3453,
     'Os':3300,
     'Ir':2716,
     'Pt':2045,
     'Au':1337.73,
     'Hg':234.43,
     'Tl':577,
     'Pb':600.75,
     'Bi':544.67,
     'Po':527,
     'At':575,
     'Rn':202,
     'Fr':300,
     'Ra':973,
     'Ac':1323,
     'Th':2028,
     'Pa':1873,
     'U':1405,
     'Np':913,
     'Pu':913,
     'SiO2':1995,
     'Gd3Ga5O12':2023,
     'LaMnO3':523,  # 523 K phase trans, melt point is 2170 K
     'La0.5Ca0.5MnO3':1300,  
     'La0.7Sr0.3MnO3':370, # ferro to para magnetic phase transition
     'LaAlO3':708,      # 708 K phase trans, melt point is 2350 K
     'B4C':2743,
     'SiC':3100,
     'YAG':2213   # https://www.crystran.co.uk/optical-materials/yttrium-aluminium-garnet-yag
}

# Boiling point in K 
boilingPoint = {'H':20.28,
     'He':4.22,
     'Li':1615,
     'Be':2742,
     'B':4200,
     'C':4300,
     'N':77.36,
     'O':90.20,
     'F':85.03,
     'Ne':27.07,
     'Na':1156,
     'Mg':1363,
     'Al':2792,
     'Si':3538,
     'P':553,
     'S':717.8,
     'Cl':239.11,
     'Ar':87.30,
     'K':1032,
     'Ca':1757,
     'Sc':3109,
     'Ti':3560,
     'V':3680,
     'Cr':2944,
     'Mn':2334,
     'Fe':3134,
     'Co':3200,
     'Ni':3186,
     'Cu':2835,
     'Zn':1180,
     'Ga':2477,
     'Ge':3106,
     'As':887,
     'Se':958,
     'Br':332,
     'Kr':119.93,
     'Rb':961,
     'Sr':1655,
     'Y':3609,
     'Zr':4682,
     'Nb':5017,
     'Mo':4912,
     'Tc':5150,
     'Ru':4423,
     'Rh':3928,
     'Pd':3236,
     'Ag':2435,
     'Cd':1040,
     'In':2345,
     'Sn':2875,
     'Sb':1860,
     'Te':1261,
     'I':457.4,
     'Xe':165.03,
     'Cs':944,
     'Ba':2170,
     'La':3737,
     'Ce':3716,
     'Pr':3793,
     'Nd':3347,
     'Pm':3273,
     'Sm':2067,
     'Eu':1802,
     'Gd':3546,
     'Tb':3503,
     'Dy':2840,
     'Ho':2993,
     'Er':3503,
     'Tm':2223,
     'Yb':1469,
     'Lu':3675,
     'Hf':4876,
     'Ta':5731,
     'W':5828,
     'Re':5869,
     'Os':5285,
     'Ir':4701,
     'Pt':4098,
     'Au':3129,
     'Hg':630,
     'Tl':1746,
     'Pb':2022,
     'Bi':1837,
     'Po':1235,
     'At':610,
     'Rn':211.3,
     'Fr':950,
     'Ra':2010,
     'Ac':3471,
     'Th':5061,
     'Pa':4300,
     'U':4404,
     'Np':4273,
     'Pu':3501
}

# Dose to reach melting point in eV/atom
meltDose = {'Li':0.043,
     'Be':0.342,
     'B':0.542,
     'C':0.23,
     'Na':0.022,
     'Mg':0.184,
     'Al':0.186,
     'Si':0.3736,
     'Ti':0.504,
     'Cr':0.734,
     'Mn':0.445,
     'Fe':0.514,
     'Co':0.554,
     'Ni':0.495,
     'Cu':0.308,
     'Zn':0.112,
     'Ga':0.0014,
     'Ge':0.2211,
     'Se':0.0407,
     'Rb':0.0046,
     'Sr':0.2388,
     'Zr':0.5611,
     'Nb':0.7986,
     'Mo':0.9321,
     'Ru':0.583,
     'Rh':0.5,
     'Ag':0.2462,
     'Cd':0.0798,
     'Te':0.1133,
     'Cs':0.0013,
     'Ba':0.2848,
     'Hf':0.7962,
     'Ta':0.9845,
     'W':1.2376,
     'WC':0.9481,
     'Pt':0.4675,
     'Au':0.2736,
     'Pb':0.0883,
     'Bi':0.0653,
     'B4C':0.6344,
     'SiC':1.006,
     'Al2O3':0.5282,
     'ZnSe':0.3781,
     'ZnTe':0.335,
     'CdS':0.359,
     'CdSe':0.307,
     'CdTe':0.272,
     'BN':0.306,
     'GaSb':0.217,
     'GaAs':0.318,
     'GaP':0.326,
     'InP':0.217,
     'InAs':0.242,
     'InSb':0.088,
     'TaC':1.252,
     'TiB2':0.871,
     'YAG':0.449,
     'Y3Al5O12':0.449,
     'CuBe':0.308,
     'ZnO':0.405,
     'SiO2':1.22,
     'AlN':0.213,
     'Si3N4':0.187,
     'ZrO2':0.257,
     'CaSiO5':0.3,
     '304SS':0.28,
     'CaSO4H4O2':0.581,
} 


# Debye temperature in K
debyeTemp = {
     'Li':344,
     'Be':1440,
     'C':2230,
     'Ne':75,
     'Na':158,
     'Mg':400,
     'Al':428,
     'Si':645,
     'Ar':92,
     'K':91,
     'Ca':230,
     'Sc':360,
     'Ti':420,
     'V':380,
     'Cr':630,
     'Mn':410,
     'Fe':470,
     'Co':445,
     'Ni':450,
     'Cu':343,
     'Zn':327,
     'Ga':320,
     'Ge':374,
     'As':282,
     'Se':90,
     'Kr':72,
     'Rb':56,
     'Sr':147,
     'Y':280,
     'Zr':291,
     'Nb':275,
     'Mo':450,
     'Ru':600,
     'Rh':480,
     'Pd':274,
     'Ag':225,
     'Cd':209,
     'In':108,
     'Sn':200,
     'Sb':211,
     'Te':153,
     'Xe':64,
     'Cs':38,
     'Ba':110,
     'La':142,
     'Gd':200,
     'Dy':210,
     'Yb':120,
     'Lu':210,
     'Hf':252,
     'Ta':240,
     'W':400,
     'Re':430,
     'Os':500,
     'Ir':420,
     'Pt':240,
     'Au':165,
     'Hg':71.9,
     'Tl':78.5,
     'Pb':105,
     'Bi':119,
     'Th':163,
     'U':207,
     'ZnSe':400,
     'ZnTe':223,
     'CdS':219,
     'CdSe':181,
     'CdTe':200,
     'BN':1900,
     'GaSb':265,
     'GaAs':344,
     'GaP':446,
     'InP':321,
     'InAs':249,
     'InSb':202,
     'LaAlO3':720,
     'LaMnO3':500,
}


# Crytal Lattice type
# cubic = cubic
# bcc = body centered cubic
# fcc = face centered cubic
# diamond = diamond
# hcp = hexaganol closed packed
# hex = hexagonal
# rhomb = rhombohedral
# zinc = zinc blende
latticeType = {'H' :'hcp',
     'He':'hcp',
     'Li':'bcc',
     'Be':'hcp',
     'B':'rhomb',
     'C':'diamond',
     'N':'cubic',
     'Ne':'fcc',
     'Na':'bcc',
     'Mg':'hcp',
     'Al':'fcc',
     'Si':'diamond',
     'Ar':'fcc',
     'K':'bcc',
     'Ca':'fcc',
     'Sc':'hcp',
     'Ti':'hcp',
     'V':'bcc',
     'Cr':'bcc',
     'Mn':'cubic',
     'Fe':'bcc',
     'Co':'hcp',
     'Ni':'fcc',
     'Cu':'fcc',
     'Zn':'hcp',
     'Ge':'diamond',
     'As':'rhomb',
     'Se':'hex',
     'Kr':'fcc',
     'Rb':'bcc',
     'Sr':'fcc',
     'Y':'hcp',
     'Zr':'hcp',
     'Nb':'bcc',
     'Mo':'bcc',
     'Tc':'hcp',
     'Ru':'hcp',
     'Rh':'fcc',
     'Pd':'fcc',
     'Ag':'fcc',
     'Cd':'hcp',
     'In':'tetr',
     'Sn':'diamond',
     'Sb':'rhomb',
     'Te':'hex',
     'Xe':'fcc',
     'Cs':'bcc',
     'Ba':'bcc',
     'La':'hex',
     'Ce':'fcc',
     'Pr':'hex',
     'Nd':'hex',
     'Eu':'bcc',
     'Gd':'hcp',
     'Tb':'hcp',
     'Dy':'hcp',
     'Ho':'hcp',
     'Er':'hcp',
     'Tm':'hcp',
     'Yb':'fcc',
     'Lu':'hcp',
     'Hf':'hcp',
     'Ta':'bcc',
     'W':'bcc',
     'Re':'hcp',
     'Os':'fcc',
     'Ir':'fcc',
     'Pt':'fcc',
     'Au':'fcc',
     'Hg':'rhomb',
     'Tl':'hcp',
     'Pb':'fcc',
     'Bi':'rhomb',
     'Po':'cubic',
     'Ac':'fcc',
     'Th':'fcc',
     'Pa':'tetr',
     'ZnSe':'zinc',
     'ZnTe':'zinc',
     'CdS':'zinc',
     'CdSe':'zinc',
     'CdTe':'zinc',
     'BN':'zinc',
     'GaSb':'zinc',
     'GaAs':'zinc',
     'GaMnAs':'zinc',
     'GaP':'zinc',
     'Gd3Ga5O12':'cubic',
     'InP':'zinc',
     'InAs':'zinc',
     'InSb':'zinc',
     'LaMnO3':'ortho',
     'LaAlO3':'rhomb',
     'La0.7Sr0.3MnO3':'rhomb',
     'GGG':'cubic',
     'YAG':'cubic',
     'Y3Al5O12':'cubic'
}


# Crystal Lattice parameters (a, b, c, alpha, beta, gamma)
# a,b,c in angstroms
# alpha, beta, gamma in degrees
latticeParameters = {
     'H' :(3.75,3.75,6.12,90,90,120),
     'He':(3.57,3.57,5.83,90,90,120),
     'Li':(3.491,3.491,3.491,90,90,90),
     'Be':(2.2866,2.2866,3.5833,90,90,120),
     'B':(5.06,5.06,5.06,58.06,58.06,58.06),
     'C':(3.567,3.567,3.567,90,90,90),
     'N':(5.66,5.66,5.66,90,90,90),
     'Ne':(4.66,4.66,4.66,90,90,90),
     'Na':(4.225,4.225,4.225,90,90,90),
     'Mg':(3.21,3.21,5.21,90,90,120),
     'Al':(4.05,4.05,4.05,90,90,90),
     'Si':(5.4310205,5.4310205,5.4310205,90,90,90),
     'Ar':(5.31,5.31,5.31,90,90,90),
     'K':(5.225,5.225,5.225,90,90,90),
     'Ca':(5.58,5.58,5.58,90,90,90),
     'Sc':(3.31,3.31,5.27,90,90,120),
     'Ti':(2.95,2.95,4.68,90,90,120),
     'V':(3.03,3.03,3.03,90,90,90),
     'Cr':(2.88,2.88,2.88,90,90,90),
     'Fe':(2.87,2.87,2.87,90,90,90),
     'Co':(2.51,2.51,4.07,90,90,120),
     'Ni':(3.52,3.52,3.52,90,90,90),
     'Cu':(3.61,3.61,3.61,90,90,90),
     'Zn':(2.66,2.66,4.95,90,90,120),
     'Ge':(5.658,5.658,5.658,90,90,90),
     'As':(4.1018,4.1018,4.1018,54.554,54.554,54.554),
     'Kr':(5.64,5.64,5.64,90,90,90),
     'Rb':(5.585,5.585,5.585,90,90,90),
     'Sr':(6.08,6.08,6.08,90,90,90),
     'Y':(3.65,3.65,5.73,90,90,120),
     'Zr':(3.23,3.23,5.15,90,90,120),
     'Nb':(3.3,3.3,3.3,90,90,90),
     'Mo':(3.15,3.15,3.15,90,90,90),
     'Tc':(2.74,2.74,4.4,90,90,120),
     'Ru':(2.71,2.71,4.28,90,90,120),
     'Rh':(3.8,3.8,3.8,90,90,90),
     'Pd':(3.89,3.89,3.89,90,90,90),
     'Ag':(4.09,4.09,4.09,90,90,90),
     'Cd':(2.98,2.98,5.62,90,90,120),
     'In':(3.25,3.25,4.95,90,90,90),
     'Sn':(6.49,6.49,6.49,90,90,90),
     'Sb':(4.4898,4.4898,4.4898,57.233,57.233,57.233),
     'Xe':(6.13,6.13,6.13,90,90,90),
     'Cs':(6.045,6.045,6.045,90,90,90),
     'Ba':(5.02,5.02,5.02,90,90,90),
     'Ce':(5.16,5.16,5.16,90,90,90),
     'Eu':(4.58,4.58,4.58,90,90,90),
     'Gd':(3.63,3.63,5.78,90,90,120),
     'Tb':(3.6,3.6,5.7,90,90,120),
     'Dy':(3.59,3.59,5.65,90,90,120),
     'Ho':(3.58,3.58,5.62,90,90,120),
     'Er':(3.56,3.56,5.59,90,90,120),
     'Tm':(3.54,3.54,5.56,90,90,120),
     'Yb':(5.45,5.45,5.45,90,90,90),
     'Lu':(3.5,3.5,5.55,90,90,120),
     'Hf':(3.19,3.19,5.05,90,90,120),
     'Ta':(3.3,3.3,3.3,90,90,90),
     'W':(3.16,3.16,3.16,90,90,90),
     'Re':(2.76,2.76,4.46,90,90,120),
     'Os':(2.74,2.74,4.32,90,90,120),
     'Ir':(3.84,3.84,3.84,90,90,90),
     'Pt':(3.92,3.92,3.92,90,90,90),
     'Au':(4.08,4.08,4.08,90,90,90),
     'Tl':(3.46,3.46,5.52,90,90,120),
     'Pb':(4.95,4.95,4.95,90,90,90),
     'Bi':(4.7236,4.7236,4.7236,57.35,57.35,57.35),
     'Po':(3.34,3.34,3.34,90,90,90),
     'Ac':(5.31,5.31,5.31,90,90,90),
     'Th':(5.08,5.08,5.08,90,90,90),
     'Pa':(3.92,3.92,3.24,90,90,90),
     'ZnSe':(5.6676,5.6676,5.6676,90,90,90),
     'ZnTe':(6.101,6.101,6.101,90,90,90),
     'CdS':(5.832,5.832,5.832,90,90,90),
     'CdSe':(6.05,6.05,6.05,90,90,90),
     'CdTe':(6.477,6.477,6.477,90,90,90),
     'BN':(3.615,3.615,3.615,90,90,90),
     'GaSb':(6.0954,6.0954,6.0954,90,90,90),
     'GaAs':(5.65315,5.65315,5.65315,90,90,90),
     'GaMnAs':(5.65,5.65,5.65,90,90,90),
     'GaP':(5.4505,5.4505,5.4505,90,90,90),
     'InP':(5.86875,5.86875,5.86875,90,90,90),
     'InAs':(6.05838,6.05838,6.05838,90,90,90),
     'InSb':(6.47877,6.47877,6.47877,90,90,90),
     'LaMnO3':(5.531,5.602,7.742,90,90,90),  
     'LaAlO3':(5.377,5.377,5.377,60.13,60.13,60.13),
     'La0.7Sr0.3MnO3':(5.4738,5.4738,5.4738,60.45,60.45,60.45),
     'Gd3Ga5O12':(12.383,12.383,12.383,90,90,90),
     'YAG':(12.006,12.006,12.006,90,90,90),  # from https://x-server.gmca.aps.anl.gov/cgi/www_form.exe?template=x0h_form.htm
     'Y3Al5O12':(12.006,12.006,12.006,90,90,90)
}

 
# specific heat capacity coefficients
# Shomate equation
# Cp = A + B*T + C*T^2 + D*T^3 + E/T^2
#   T = temperature(K)/1000
specificHeatParams = {
    'Li':(169.552,-882.711,1977.438,-1487.312,-1.609635),
    'Be':(21.20694,5.68819,0.968019,-0.001749,-0.587526),
     'B':(10.18574,29.24415,-18.02137,4.212326,-0.551),
    'Na':(72.63675,-9.491572,-730.9322,1414.518,-1.259377),
    'Mg':(26.54083,-1.533048,8.062443,0.57217,-0.174221),
    'Al':(28.0892,-5.414849,8.560423,3.42737,-0.277375),
    'Si':(22.81719,3.89951,-0.082885,0.042111,-0.354063),
    'Ti':(23.0566,5.541331,-2.055881,1.611745,-0.056075),
    'Cr':(7.489737,71.50498,-91.67562,46.0445,0.138157),
    'Mn':(27.2419,5.23764,7.78316,-2.118501,-0.282113),
    'Fe':(23.97449,8.36775,0.000277,-0.000086,-0.000005),
    'Co':(10.9943,54.375,-55.5132,25.817,0.164533),
    'Ni':(13.6916,82.49509,-174.9548,161.6011,-0.092417),
    'Cu':(17.72891,28.0987,-31.25289,13.97243,0.068611),
    'Zn':(25.60123,-4.405292,20.42206,-7.399697,-0.045801),
    'Ga':(102.3394,-347.5134,603.3621,-360.7047,-1.490304),
    'Ge':(23.3667,0,0,0,0),
    'Se':(25.363,0,0,0,0),
    'Rb':(9.44626,65.31182,45.5123,-26.78961,-0.10773),
    'Sr':(23.88223,9.297351,0.919924,0.035243,0.004934),
    'Zr':(25.3924,0.434236,4.384539,1.017123,-0.065449),
    'Nb':(22.0143,9.88816,-5.64853,1.759691,0.021839),
    'Mo':(24.72736,3.960425,-1.270706,1.153065,-0.170246),
    'Rh':(24.9,0,0,0,0),
    'Ag':(25.35,0,0,0,0),
    'Cd':(26.02,0,0,0,0),
    'Te':(25.73,0,0,0,0),
    'Cs':(57.04424,-50.0034,48.595,-16.72822,-1.223804),
    'Ba':(83.8018,-406.186,915.066,-519.805,-0.191854),
    'Hf':(22.71033,10.86443,-1.901809,0.306368,-0.007587),
    'Ta':(20.69482,17.29992,-15.68987,5.608694,0.061581),
     'W':(23.9593,2.63968,1.25775,-0.254642,-0.048407),
    'Pt':(25.86,0,0,0,0),
    'Au':(25.409,0,0,0,0),
    'Pb':(25.0145,5.441836,4.061367,-1.236214,-0.010657),
    'Bi':(25.52,0,0,0,0),
    'B4C':(95.99853,23.16513,-0.409604,0.081414,-4.395208),
    'SiC':[[298, 1000, (20.55859,64.57962,-52.98827,16.95813,-0.781847)],
           [1000,4000, (46.90222,5.845968,-1.085410,0.093021,-3.448876)]],
    'BeO':[[298,  800, (3.358974,131.5922,-140.4937,56.20953,-0.536669)],
           [800, 2780, (47.06205,5.598359,-0.495570,0.054527,-2.947373)]],
    'Al2O3':(102.429,38.7498,-15.9109,2.628181,-3.007551),
    'ZnSe':(48.9,0,0,0,0),
    'ZnTe':(50.95,0,0,0,0),
    'CdS':(47.67,0,0,0,0),
    'CdSe':(48.8,0,0,0,0),
    'CdTe':(49.2,0,0,0,0),
    'BN':(19.68,0,0,0,0),
    'GaSb':(61.3,0,0,0,0),
    'GaAs':(50.6,0,0,0,0),
    'GaP':(43.3,0,0,0,0),
    'InSb':(61.3,0,0,0,0),
    'InAs':(50.9,0,0,0,0),
    'InP':(34,0,0,0,0),
    'TaC':(44.29224,7.673707,-0.091309,0.010861,-0.875548),
    'TiB2':(52.33264,33.69484,-7.909266,0.803989,-1.540905),
    'YAG':(376,0,0,0,0),
    'Y3Al5O12':(376,0,0,0,0),
    'ZnO':(40.2,0,0,0,0),
    'CaSiO5':(111,0,0,0,0),
    'H2O':(75.327,0,0,0,0),
    'CO2':(36.94,0,0,0,0),
    'Fe.68Cr.2Ni.1Mn.02':(27.553,0,0,0,0),
    'SiO2':(-6.07659,251.6755,-324.796,168.5604,0.002548),
    'LaAlO3':(86.6,0,0,0,0),
    'LaMnO3':(89,0,0,0,0),
    'La0.5Ca0.5MnO3':(89,0,0,0,0)
}

# List of specific heat
#  [K, J/mol/K, Integrated J/mol]
specificHeat = {
    'CVD': np.array([[50, 100, 150, 200, 250, 298, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000],
                     [0.052335, 0.20850264, 1.00022652, 2.3801958, 4.15456164, 6.07797756, 6.16045752, 8.21575764, 10.17643608, 11.96252496, 13.544298, 14.92217388, 16.11206244, 17.1365724, 18.01831248, 18.77821668, 19.43470692, 20.00411172, 20.50024752, 20.93358132, 21.31416144, 21.6499428, 22.44627216, 22.8452742, 23.1697512, 23.43603168, 23.6575134, 23.99999364, 24.24952692, 24.43667688, 24.58028412, 24.69290904, 24.78292524, 24.82521192],
                     [0.00, 6.52, 36.74, 121.25, 284.62, 530.20, 542.44, 901.84, 1361.65, 1915.12, 2552.79, 3264.45, 4040.31, 4871.53, 5750.40, 6670.31, 7625.63, 8611.61, 9624.21, 10660.06, 11716.25, 13864.46, 16069.27, 18333.85, 20634.60, 22964.89, 25319.56, 30085.32, 34910.27, 39778.89, 44680.58, 49607.90, 54555.49, 59516.30]
                    ]),
    # Graphite: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782425&Units=SI&Mask=2#Thermo-Condensed
    # T = 200 to 3500 K. Least squares fit of 'best' data gives: Cp = 0.538657 + 9.11129x10-6*T - 90.2725*T^-1 - 43449.3*T^-2 + 1.59309x10^7*T^-3 - 1.43688x10^9*T^-4 cal/g*K (250 to 3000 K)
    'Graphite': np.array([[300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 3900],
                          [8.55, 10.30, 11.94, 13.42, 14.73, 15.86, 16.85, 17.71, 18.47, 19.13, 19.72, 20.25, 20.71, 21.14, 21.52, 22.18, 22.73, 23.20, 23.61, 23.96, 24.28, 24.81, 25.24, 25.61, 25.93, 26.21, 26.46, 26.69, 26.90, 27.10, 27.28, 27.45, 27.54],
                          [0, 471, 1028, 1663, 2367, 3133, 3951, 4816, 5721, 6661, 7633, 8633, 9657, 10703, 11770, 13955, 16201, 18499, 20839, 23218, 25631, 30541, 35547, 40634, 45789, 51004, 56272, 61588, 66948, 72348, 77786, 83260, 86009]
                         ])
}

# Latent heat in kJ/mol
# (Heat of Fusion, Heat of Vaporization)
# Heat of Vaporization: https://en.wikipedia.org/wiki/Enthalpy_of_vaporization
latentHeat = {
     'Fe':(13.81, 340),
    'H2O':(6.009, 40.66)
}
