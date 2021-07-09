from distutils.spawn import find_executable
import os

# ChimeraX executable path
# CHIMERA_PATH = '/nfs/msd/work2/andrei/soft/chimerax/18-12/redhat/opt/UCSF/ChimeraX-daily/bin/ChimeraX'
CHIMERA_PATH = '' #'/Applications/ChimeraX-1.1.1.app/Contents/MacOS/ChimeraX'


def check_chimerax_executable():
    out = find_executable(CHIMERA_PATH)
    return out


def no_chimerax_warning():
    message = 'WARNING: could not find ChimeraX\n please run:\n'\
              'strudel_setChimeraX.py path_to_ChimeraX_executable'
    return message


def no_chimerax_error():
    message = 'ERROR: could not find ChimeraX\n please run:\n'\
              'strudel_setChimeraX.py path_to_ChimeraX_executable'
    return message


if not check_chimerax_executable():
    print(no_chimerax_warning())