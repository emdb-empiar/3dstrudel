#!/Users/andrei/PycharmProjects/for_testing/venv/bin/python

import threed_strudel
import os
import sys
from distutils.spawn import find_executable

config_path = os.path.abspath(os.path.join(os.path.dirname(threed_strudel.__file__), 'configure.py'))

CHIMERA_PATH = sys.argv[1]

out = find_executable(CHIMERA_PATH)
if out is None:
    print(f'ERROR: \n'
          f'{CHIMERA_PATH}'
          f'\n is not an executable')

else:
    with open(config_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            words = line.split()
            if len(words) > 1:
                if words[0] == 'CHIMERA_PATH':
                    words[-1] = f'"{CHIMERA_PATH}"'
                    lines[i] = ' '.join(words) + '\n'
    with open(config_path, 'w') as f:
        f.writelines(lines)