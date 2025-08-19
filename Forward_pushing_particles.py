"""
Simulate charged-particle propagation in a combined Galactic magnetic field
using CRPropa and record arrivals at a small spherical observer at the Galactic 
Centre (GC). We use turbulent fields from CRPropa3 as an example but this can 
be adapted for any model. JF12 is also provided as an example.

Usage
-----
python this_script.py <counter>

Parameters (File counter)
----------------------
- counter (sys.argv[1]): string/number appended to the output filename so you
  can differentiate multiple runs.

Key choices (JF12 as an example):
--------------------------------------
- Regular field: JF12 with amplitude bpeak [µG].
- Turbulent field: Kolmogorov spectrum (index -11/3), RMS bTur [µG] on a 3D grid.
- Source: Lambert distribution on a sphere of radius 'radius' (kpc), centered at GC,
  emitting inward.
- Observer #1 (active): small sphere at GC (radius ≈ 1 kpc) to collect arrivals.
- Observer #2 (trash): spherical surface at (radius + 1) kpc that deactivates CRs
  leaving the Galaxy (not written to file).
- Propagation: PropagationBP with min/max step sizes of 0.1 pc and 100 pc.

Outputs
-------
A text file of detections from the GC observer:
  Ver_3_Obs_{s1}kpc_Source_Radius_{radius}_Bstr_{bpeak}_Btur_{bTur}_{counter}.txt

Dependencies
------------
- CRPropa3 Python bindings
- NumPy
"""

import numpy as np
import sys
import time  # kept for parity with your original imports
from crpropa import *

# ----------------------------- Field parameters ------------------------------

# Regular (ordered) field amplitude in microgauss (µG). 0 reproduces your default.
bpeak = 0

# Turbulent RMS field amplitude (µG). Your original scaling is preserved.
# 19.95 / sqrt(1000/150) ≈ 7.732 µG
bTur  = 19.95/np.sqrt(1000/150)

# JF12 positional offsets (kpc). Both set to 0 to center on the GC.
Rr    =  0
Z     = 0

# ----------------------------- Run identifier --------------------------------

# Counter used to disambiguate outputs on disk (taken from CLI).
# Example usage: python script.py 42
counter = float(sys.argv[1])

# ----------------------------- Magnetic fields --------------------------------

# Using a full JF12 field. Arguments:

#BToy  = JF12Field() ## Using normal JF12 field
#seed = 691342
#BToy.randomStriated(seed)
#BToy.randomTurbulent(seed) 
# Turbulence parameters
randomSeed = 10                         # RNG seed for reproducibility
BTur = bTur                             # alias to keep naming from original code
vgrid = Grid3f(Vector3d(10*kpc,10*kpc,10*kpc), 201, 100*pc)  # box, resolution, cell size

# Turbulent spectrum scales (pc)
lmin = 450*pc
lmax = 4000*pc

# Populate the turbulence grid with a Kolmogorov spectrum (index -11/3).
initTurbulence(vgrid, bTur*muG, lmin, lmax, -11./3., randomSeed)
bField_tur = MagneticFieldGrid(vgrid)

# Combine ordered + turbulent components into a single field list.
field_combined = MagneticFieldList()
field_combined.addField(BToy)
field_combined.addField(bField_tur)

# ----------------------------- I/O and run size ------------------------------

# Output directory (preserved exactly as provided).
dd = '/your_desired_location/'

# Number of injected particles.
n = 5*10**6

# Particle identity: negative nucleusId(Z=1, A=1) per your original line.
# (Using a negative Id is allowed in CRPropa to encode particle/antiparticle.)
pid = -nucleusId(1, 1)

# Injection energy (mono-energetic), in EeV.
meanEnergy = 8.5*EeV

# ----------------------------- Propagation setup -----------------------------

sim = ModuleList()

# PropagationBP(field, epsilon, minStep, maxStep)
# epsilon = 1e-4 (tolerance), min = 0.1 pc, max = 100 pc
sim.add(PropagationBP((field_combined),1e-4, 0.1 * parsec, 100 * parsec))

# ----------------------------- Observers -------------------------------------

# Primary observer: collect arriving cosmic rays at the Galactic Centre (GC).
obs_1 = Observer()

# Placeholder for an additional observer (not used). Kept for parity.
obs_out = Observer()

# Galactic Centre position (in kpc units once multiplied).
pos_GC = Vector3d(0, 0, 0) * kpc

# Source/termination geometry:
# 'radius' is the source-sphere radius (kpc). 'inward=True' emits toward the center.
center, radius, inward = Vector3d(0, 0, 0) ,30 , True

# Small-sphere observer radius at the GC.
# Your expression 0.15*(1000/150) evaluates to exactly 1.0 kpc.
s1 = 0.15*(1000/150)

# Small spherical observer at GC, radius s1 kpc.
obs_1.add(ObserverSmallSphere(pos_GC, s1 * kpc))

# Output file name mirrors your original convention.
filename_output_1 = dd+'Ver_3_Obs_'+str(s1)+'kpc_Source_Radius_'+str(radius)+'_Bstr_'+str(bpeak)+'_Btur_'+str(bTur)+'_'+str(counter)+'.txt'
obs_1.onDetection( TextOutput(filename_output_1 ) )
obs_1.setDeactivateOnDetection(False)

# Register the primary observer in the simulation chain.
sim.add(obs_1)

# ----------------------------- Escape "trash" surface ------------------------

# This observer deactivates outward-going CRs that cross the outer boundary sphere,
# preventing them from propagating further once they leave the Galaxy model volume.
obs_trash = Observer()
obs_trash.add(ObserverSurface(Sphere(Vector3d(0), (radius + 1) * kpc)))

# If you want to log these escapes, uncomment and provide a path:
# trash = 'Trash_Termi_Radius_'+str(radius * 10)+'_Bstr_'+str(bpeak)+'_Btur_'+str(bTur)+'_'+str(counter*1000)+'_set_2.txt'
# obs_trash.onDetection( TextOutput(trash ) )

# Deactivate particles when they hit the boundary (do not write by default).
obs_trash.setDeactivateOnDetection(True)
sim.add(obs_trash)

# ----------------------------- Source definition -----------------------------

source = Source()

# Lambert distribution on a sphere centered at 'center' (Vector3d) with radius 'radius' (kpc).
# 'inward=True' to emit toward the GC.
source.add(SourceLambertDistributionOnSphere(center* kpc, radius* kpc, inward))

# Particle type and energy (monoenergetic).
source.add(SourceParticleType(pid))
source.add(SourceEnergy(meanEnergy))

# ----------------------------- Run ------------------------------------------

# Execute the simulation with 'n' events.
sim.run(source, n)
