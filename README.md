# Cosmic_Ray_Diffusion

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)  
[![CRPropa3](https://img.shields.io/badge/CRPropa3-3.2-green.svg)](https://crpropa.github.io/CRPropa3/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository provides tools to calculate **cosmic-ray particle diffusion** using the forward-tracking scheme implemented in [**CRPropa3**](https://crpropa.github.io/CRPropa3/).  
It is intended for researchers working on cosmic-ray propagation, diffusion studies, and comparisons with analytical models.  

📄 Please see the related publication for methodological details: [arXiv:2501.16881](https://arxiv.org/abs/2501.16881).

---

## 🚀 Features
- Forward tracking of charged particles using **CRPropa3**
- Calculation of **diffusion coefficients** in turbulent magnetic fields
- Integration with **Healpy** for sky-map projections
- FITS file support via **Astropy**
- Modular and reproducible workflow

---

## 📦 Requirements

This code requires **Python 3.10+** and the following packages:

| Package   | Recommended Version |
|-----------|----------------------|
| [CRPropa3](https://crpropa.github.io/CRPropa3/pages/Installation.html) | 3.2+ |
| numpy     | 1.23+ |
| healpy    | 1.16+ |
| astropy   | 5.3+ |
| scipy     | 1.10+ |
| matplotlib| 3.7+ |

> ⚠️ Ensure that your `numpy` and other libraries are compatible with the installed CRPropa3 version.  
> Using a virtual environment (`conda` or `venv`) is **highly recommended**.

---

## ⚙️ Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/Cosmic_Ray_Diffusion.git
cd Cosmic_Ray_Diffusion

conda create -n crdiff python=3.10
conda activate crdiff

pip install numpy==1.23.5 healpy==1.16.5 astropy==5.3 scipy==1.10 matplotlib==3.7

import crpropa
import numpy as np
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt

# Define magnetic field and simulation parameters here
# Run CRPropa simulation
# Save and analyze diffusion results

