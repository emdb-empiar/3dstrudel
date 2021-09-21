import mrcfile
path = '/Users/andrei/Downloads/emd_21109.map'
with mrcfile.open(path, mode='r+', permissive=True) as mrc:
    print(mrc.header.mx, mrc.header.my, mrc.header.mz)
    print(mrc.data.shape)
    print(mrc.header)