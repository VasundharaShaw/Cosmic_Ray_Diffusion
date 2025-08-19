import crpropa
import matplotlib.pyplot as plt  # not used below, kept as in your original imports
import numpy as np
import pandas as pd  # not used below, kept as in your original imports
from crpropa import *
import sys
import time
import healpy as hp  # not used below, kept as in your original imports

deg = np.degrees
rad = np.radians
pi  = np.pi

"""
Trajectory-length capped propagation in combined Galactic magnetic field.

Overview
--------
- Builds a CRPropa magnetic field: Turbulent field here.
- Propagates mono-energetic protons launched from a spherical shell inward,
  with a fixed direction (beam-like).
- Attaches a custom module that *records and deactivates* particles once their
  trajectory length exceeds a given maximum.
- (Observer lines for detections are provided but commented out to match your original.)

Command-line arguments
----------------------
argv[1] : float
    Maximum trajectory length (in CRPropa length units; see usage note below).
argv[2] : float
    Counter/index appended to the output filename.

Output
------
A text file in `dd` named:
    TrajLen_Obs_Sph_{obs_sph}_Btr_{bTur}_len_{maxLen}_{count}.txt
containing, per line:
    SerialNumber  ri(x,y,z)  pi(x,y,z)  rf(x,y,z)  pf(x,y,z)

Units & notes
-------------
- CRPropa internal units are used (e.g., kpc, pc, muG, EeV).
- The custom module uses `Candidate.getTrajectoryLength()` and compares
  against the provided `traj_len`. Pass `traj_len` in the same units as returned
  by CRPropa (typically the same unit system you build with, e.g. kpc if you divide by kpc).
"""

# Base directory for output files (kept exactly as provided)
dd = '/lustre/fs24/group/that/work_vasu/ArrivalDirections/BF_Dipole_Change_Source_BTur_GC/Source_50kpc/30kpc/19.95/SS_Data_Feb_23/'


class RecordMaxLength(Module):
    """
    Export info for particles that reach/exceed a maximum trajectory length
    and then deactivate them to stop further transport.

    Parameters
    ----------
    max_length : float
        Maximum allowed trajectory length in CRPropa length units
        (compare to `c.getTrajectoryLength()`).
    Notes
    -----
    - Output filename depends on globals `obs_sph`, `bTur`, and `count` which
      are set in this script prior to instantiation.
    - On exceeding the limit, the particle is written out and deactivated via
      `c.setActive(False)`.
    """

    def __init__(self, max_length):
        Module.__init__(self)

        self.max_length = max_length
        # Open output file (includes observer radius, turbulence RMS, length, and run count)
        self.output_file = open(
            dd + f'TrajLen_Obs_Sph_{obs_sph}_Btr_{bTur}_len_{self.max_length:3.2e}_{count}.txt',
            'w'
        )
        # Header (commented to keep file minimal; uncomment if desired)
        # self.output_file.write('SerialNumber InitialPos(x,y,z) InitialMom(x,y,z) FinalPos(x,y,z) FinalMom(x,y,z)\n')

    def process(self, c):
        """
        Called by CRPropa during propagation. If the candidate's trajectory
        length exceeds the threshold, write its state and deactivate it.
        """
        # Compare to max length (normalize by kpc here to keep check in kpc units)
        if c.getTrajectoryLength() / kpc >= self.max_length:
            # Snapshot identifiers and state
            sn = c.getSerialNumber()
            pi = c.source.getMomentum()
            ri = c.source.getPosition()
            pf = c.current.getMomentum()
            rf = c.current.getPosition()

            # Write line: serial, initial pos/mom, final pos/mom
            self.output_file.write(
                f'{sn:d} '
                f'{ri.x:6.5e} {ri.y:6.5e} {ri.z:6.5e} {pi.x:6.5e} {pi.y:6.5e} {pi.z:6.5e} '
                f'{rf.x:6.5e} {rf.y:6.5e} {rf.z:6.5e} {pf.x:6.5e} {pf.y:6.5e} {pf.z:6.5e}\n'
            )

            # Deactivate particle to exclude from further propagation
            c.setActive(False)


# ----------------------------- Magnetic field setup -----------------------------


counter = 0        # extra index (not used in filename here; 'count' from argv is)

# Turbulence grid and spectrum
randomSeed = 10
BTur = bTur
# Grid3f(origin_vector, resolution, cell_size). Using a cube spanning ~20 kpc.
vgrid = Grid3f(Vector3d(-10 * kpc, -10 * kpc, -10 * kpc), 201, 100 * pc)

# Turbulent coherence range (pc)
lmin = 200 * pc
lmax = 400 * pc

# Populate the turbulence grid with Kolmogorov spectrum (index -11/3)
initTurbulence(vgrid, bTur * muG, lmin, lmax, -11.0 / 3.0, randomSeed)
bField_tur = MagneticFieldGrid(vgrid)

# Combine ordered + turbulent field components
field_combined = MagneticFieldList()
field_combined.addField(bField_tur)

# ----------------------------- CLI arguments & run size ------------------------

# Maximum trajectory length (kpc, since process() divides by kpc) and run counter
traj_len = float(sys.argv[1])
count    = float(sys.argv[2])

# Number of injected particles
n = 10 ** 4

# ----------------------------- Particle & propagation --------------------------

# Particle ID: negative nucleusId(1,1) (proton with sign choice matching your original)
pid = -crpropa.nucleusId(1, 1)
meanEnergy = 8.5 * EeV  # injection energy

sim = crpropa.ModuleList()
# PropagationBP(field, epsilon, minStep, maxStep): 1e-4 tolerance, [0.1 pc, 100 pc] step range
sim.add(PropagationBP(field_combined, 1e-4, 0.1 * parsec, 100 * parsec))

# ----------------------------- Observer geometry -------------------------------

obs = crpropa.Observer()
pos_GC = crpropa.Vector3d(0, 0, 0) * crpropa.kpc
obs_sph = 0.15  # observer sphere radius (kpc) at GC
obs.add(crpropa.ObserverSurface(crpropa.Sphere(pos_GC, obs_sph * kpc)))

# If you want to write detections and deactivate on hit, uncomment:
# filename_output = 'HitsObserver.txt'
# detectoutput = TextOutput(filename_output)
# obs.onDetection(detectoutput)
# obs.setDeactivateOnDetection(True)
# sim.add(obs)

# ----------------------------- Source definition -------------------------------

source = crpropa.Source()
center, radius, inward = crpropa.Vector3d(0, 0, 0) * crpropa.kpc, 30 * crpropa.kpc, True

# Narrow-beam like injection: fixed direction + uniform shell radius
source.add(crpropa.SourceDirection(Vector3d(0, -1, 0)))
source.add(crpropa.SourceUniformShell(center, radius))

source.add(crpropa.SourceParticleType(pid))
source.add(crpropa.SourceEnergy(meanEnergy))

# ----------------------------- (Optional) escape surface -----------------------
# This block collects particles crossing a large spherical boundary and deactivates them.
# Kept commented-out to mirror your original.

# obs_trash = crpropa.Observer()
# obs_trash.add(crpropa.ObserverSurface(crpropa.Sphere(crpropa.Vector3d(0), (radius * 2))))
# Trashfilename = f'1e2_Trash_100pc_Grid_engine_output_Bpeak_{bpeak}_Btur_{bTur}_meanEnergy_{meanEnergy/EeV}_num_{counter}_.txt'
# trashoutput = TextOutput(Trashfilename)
# obs_trash.onDetection(trashoutput)
# obs_trash.setDeactivateOnDetection(True)
# sim.add(obs_trash)

# ----------------------------- Trajectory-length capping -----------------------

mod = RecordMaxLength(traj_len)
sim.add(mod)

# Progress bar during propagation
sim.setShowProgress(True)

# Run the simulation
sim.run(source, n)

# Close any outputs opened in custom modules
# (If you enabled the observer TextOutput above, remember to close it too.)
# detectoutput.close()
mod.output_file.close()
