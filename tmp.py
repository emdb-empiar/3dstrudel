import threed_strudel.utils.bio_utils as bu

path = '/Volumes/data/strudel_validation/results/workdir_0.5px/cov_30342/out/input/7cec.cif'


# struct = bu.load_structure(path)
# residues = struct[0]["H"].get_list()
# print(dir(struct[0]["H"]))
# res = struct[0]["H"]
# print(res.get_unpacked_list())
# print(res.get_list())
#
# print(res[31])
# print(res[31].is_disordered())
# for r in res.get_list():
#     print(r)
#     print(r.is_disordered())
#     for atom in r:
#         print(atom.is_disordered())
#         print(atom.altloc)




            # print(r)
#     except:
#         pass
#     for atom in r.get_list():
#         if atom.is_disordered():
#             print(r, atom)


# for res in residues:
#     if res.is_disordered():
#         print(res)
# for i in range(0,36):
#     res = struct[0]['H'][i]
#     print(res)
#     print(res.is_disordered())

# print(struct[0]['H'][31])
# print(struct[0]['H'][31].is_disordered())
# res31 = struct[0]['H'][31]
# print(res31.get_unpacked_list())
# for atom in struct[0]['H'][31]:
#
#     print(atom)
#     print(atom.is_disordered())
#     print(atom.get_altloc())
#     print(atom.altloc)

# struct = bu.load_structure_label_id(path)
# for i in range(0,36):
#     print(struct[0]['H'][i])
#
# for atom in struct[0]['H'][31]:
#     print(atom)

# path1 = '/Users/andrei/Downloads/1test.cif'
# struct1 = bu.load_structure_label_id(path1)[0]
# for chain in struct1:
#     for r in chain:
#         print(r)

path = '/Volumes/data/strudel_validation/results/workdir_0.5px/T0102_err/out/input/T0102EM054_1_err.cif'

struct = bu.load_structure(path)
