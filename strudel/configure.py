from distutils.spawn import find_executable
import os

# ChimeraX executable path
CHIMERA_PATH = '/nfs/msd/work2/andrei/soft/chimerax/18-12/redhat/opt/UCSF/ChimeraX-daily/bin/ChimeraX'
# CHIMERA_PATH = '/Applications/ChimeraX_Daily.app/Contents/MacOS/Chimera'

out = find_executable(CHIMERA_PATH)
if out is None:
    print(f'WARNING: could not find ChimeraX\n please edit {os.path.abspath(__file__)}')