# X-Ray Tools
Set of tools for x-ray analysis.

## xraytools.py
Various tools for x-ray calculations
* Dependency: [xraylib](https://github.com/tschoonj/xraylib), numpy, scipy, matplotlib, requests.
    * Maybe easier to install xraylib via Anaconda, see below or more in [the page](https://anaconda.org/conda-forge/xraylib):
      ```
      conda install conda-forge::xraylib
      ```
* f1f2_EPDL97.dat & f1f2_Windt.dat: Atomic scattering factor data. "EPDL97" is the default.
* f2d_photon: Fluence to dose coversion factor for photons.

## web
Webpage tools for x-ray analysis (under construction).
* shield11: Web version SHIELD11 code for high energy electron shielding.
* sr_sheet: Synchrotron radiation calculator.
* transmission_solid: Calculate transmission through solid materials.
