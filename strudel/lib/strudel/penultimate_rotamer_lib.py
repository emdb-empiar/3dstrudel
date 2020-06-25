"""
The Penultimate Rotamer Library
Simon C. Lovell, J. Michael Word, Jane S. Richardson, and David C. Richardson
PROTEINS: Structure, Function, and Genetics 40:389–408 (2000)
"""
import json
import os


basedir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(basedir, 'penultimate_lib.json')) as j:
    PENULTIMATE_2000 = json.load(j)
