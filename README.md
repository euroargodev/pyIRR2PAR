# pyIRR2PAR
Calculation of PAR from multispectral irradiance profiles according to Argo recommendations (Python version)

**⚠️ Under Validation (Do Not Use Yet)**
This code is currently under review and validation. It is not yet ready for use. Please check back later for updates

the Matlab version is [here](https://github.com/euroargodev/matIRR2PAR)

# Objectives and Background
The Photosynthetically Available Radiation (PAR) is a radiometric parameter that has been measured on floats since the beginning of the BGC mission. It is particularly important to estimate Primary Production or to help for Chla processing. This parameter was traditionally measured using the OCR504 radiometer and was typically associated with wavelengths:  380, 412, and 490 nm. In accordance with the recommendations approved during AST-24, this configuration was replaced by: 380, 443, 490, and 555 nm (without direct PAR measurement) due to new scientific applications and the potential ability to derive PAR measurements from the four measured wavelengths [Organelli_WG_Radiometry_AST24_20March2023.pdf](https://drive.google.com/file/d/1HhDM9NZMGbhgXSlx_RA4Nooc2g_KlGUH/view?usp=sharing). 

The purpose of this Python code is to provide a wrapper for the two models evaluated within the Argo framework. It uses as input a matrix of four irradiance values at **380, 443, 490, and 555 nm**, along with the associated depth vector, and calculates the outputs of both models as well as the average of the calculated PAR values, which is the PAR estimate to be considered as recommended during ADMT26 [ADMT26-Leymarie-PAR Model.pptx](https://docs.google.com/presentation/d/1edo8Na_IFxEJGKgBwDNWRhW9zqafZP9W/edit?usp=drive_link&ouid=112140535291181030156&rtpof=true&sd=true).

# How to use this package

## Installation
`pyIRR2PAR` will soon be available on the [Python Package Index (PyPI)](https://pypi.org/project/pyIRR2PAR):
```bash
pip install pyIRR2PAR
```

## Usage
```
JTan_PAR_results, JPitarch_PAR_results, mean_PAR, mean_uncertainty = IRR2PAR(Ed, depth)

Parameters:
* Ed -- 2D array with columns [Ed380, Ed443, Ed490, Ed555] in W/m2/nm
* z  -- depth in dbar

Returns:
* Two lists: first for JTan_2025 results, second for JPitarch_2025 results,
* mean_PAR : Vector of the average PAR values from both models (excluding NA/NaN values). This is the PAR value recommended under the Argo framework.
* mean_uncertainty : vector of uncertainties
```
All PAR values are in microMoleQuanta/m^2/sec.

Here is a example:
```python
import numpy as np
from pyIRR2PAR import IRR2PAR

Ed = np.array([[3.1477e-07, 4.401e-05, 0.00031093, 2.6219e-05],
               [0.0097528, 0.061085, 0.10562, 0.087316],
               [0.20572, 0.52315, 0.54544, 0.16185]])
z = np.array([84.8, 14.2, 28.1])

JTan_PAR_results, JPitarch_PAR_results, mean_PAR, mean_uncertainty = IRR2PAR(Ed, z)
```

A more realistic example:
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyIRR2PAR import IRR2PAR

# Read example data:
# 1st column - depth, 2nd - PAR (micro E/m2/s)
# 3rd - 6th columns: Ed380 (mW/cm2/micron), Ed443 (mW/cm2/micron), Ed490 (mW/cm2/micron), Ed555 (mW/cm2/micron)

df = pd.read_csv('./example/1902685-lovuse023b-cycle12.csv', sep=";", skiprows=None)
df = df.sort_values(by='depth')

#
d = df['depth'] # Depth
par = df['RamsesPAR'] # Measured PAR (micro E/m2/s)

#
profile = df[['Ed_380', 'Ed_443', 'Ed_490', 'Ed_555']].to_numpy()
JTan_PAR_results, JPitarch_PAR_results, mean_PAR, mean_uncertainty = IRR2PAR(profile, df['depth'])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12), dpi=90, layout="constrained", sharey=True)
fig.suptitle('1902685-lovuse023b-cycle12',fontsize=16)

ax1.plot(JTan_PAR_results[0], d, 'b', label='modeled (Tan et al, 2025)')
ax1.plot(JPitarch_PAR_results[0], d, 'c', label='modeled (Pitarch et al, 2025)')
ax1.plot(mean_PAR, d, 'k', label='modeled (IRR2PAR)')
ax1.plot(par, d, 'r',label='measured')
ax1.set_ylabel('Depth, m',fontsize=14)
ax1.set_xlabel('PAR, '+r'$\mu$E/m$^{2}$/s',fontsize=14)
ax1.set_xscale('log')
ax1.legend(fontsize=12)
ax1.yaxis.set_inverted(True)
ax1.grid()

ax2.plot(JTan_PAR_results[1], d, 'b', label='modeled (Tan et al, 2025)')
ax2.plot(JPitarch_PAR_results[1], d, 'c', label='modeled (Pitarch et al, 2025)')
ax2.plot(mean_uncertainty, d, 'k', label='modeled (IRR2PAR)')
ax2.set_xlabel(r'$\Delta$'+'PAR, '+r'$\mu$E/m$^{2}$/s',fontsize=14)
ax2.set_xscale('log')
ax2.legend(fontsize=12)
ax2.grid()

fig.tight_layout()
plt.show()
```
![example/example.png](./example/example.png)


# Bibliography
We would like to warmly thank Jaime Pitarch and Jing Tan for their work and for making their code available to the community.
* Pitarch, J., Leymarie, E., Vellucci, V., Massi, L., Claustre, H., Poteau, A., Antoine, D., Organelli, E., 2025. Accurate estimation of photosynthetic available radiation from multispectral downwelling irradiance profiles. Limnology and Oceanography: Methods. [https://doi.org/10.1002/lom3.10673](https://doi.org/10.1002/lom3.10673)
GitHub original repository : [PAR_BGC_Argo](https://github.com/euroargodev/PAR_BGC_Argo)
  
* Tan, J., Frouin, R., Leymarie, E., Mitchell, B.G., 2025. Modeling underwater photosynthetically available radiation profiles from biogeochemical Argo floats using multi-spectral irradiance measurements. Opt. Express, OE 33, 44355–44377. [https://doi.org/10.1364/OE.566083](https://doi.org/10.1364/OE.566083)
GitHub original repository : [BioArgo_PAR](https://github.com/jit079/BioArgo_PAR)

# How to cite this calculation
If you use this package, please cite the original work by J. Pitarch and J. Tan and note that the PAR estimate is the average of these two models.
