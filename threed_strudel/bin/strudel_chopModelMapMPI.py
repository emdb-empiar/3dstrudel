#!/usr/bin/env python

from threed_strudel.chop import chop_model_map_mpi

print('from outside', __name__)
if __name__ == '__main__':
    chop_model_map_mpi.main()