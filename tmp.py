from threed_strudel.chop.chop_map import ChopMap
import time
import os
import gzip
import shutil

model = '/Volumes/data/test_gz/7c8k.cif'
map_ = '/Volumes/data/test_gz/emd_30306.map.gz'
unz_map = '/Volumes/data/test_gz/emd_30306.map'

t = time.time()
with gzip.open(map_, 'rb') as f_in:
    with open(unz_map, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
os.remove(map_)



out_map = '/Volumes/data/test_gz/test.map'


ChopMap.chop_cube(model, unz_map, zero_origin=False, out_map_path=out_map)
print(f"Time: {time.time()-t}")