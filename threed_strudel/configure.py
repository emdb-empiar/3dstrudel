from distutils.spawn import find_executable
import os

# ChimeraX executable path
# CHIMERA_PATH = '/nfs/msd/work2/andrei/soft/chimerax/18-12/redhat/opt/UCSF/ChimeraX-daily/bin/ChimeraX'
CHIMERA_PATH = '/Applications/ChimeraX-1.0.app/Contents/MacOS/ChimeraX'

out = find_executable(CHIMERA_PATH)
if out is None:
    print('WARNING: could not find ChimeraX\n please run:\n'
          'strudel_setChimeraX.py path_to_ChimeraX_executable')