
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
     name='3dstrudel',
     version='0.1.1',
     author="Andrei Istrate",
     author_email="andrei@ebi.ac.uk",
     description="Strudel package",
     include_package_data=True,
     scripts=['3dstrudel/bin/strudel_mapMotifValidation.py', '3dstrudel/bin/strudel_mapAveraging.py',
              '3dstrudel/bin/strudel_penultimateClassifier.py', '3dstrudel/bin/strudel_chopModelMapMPI.py',
              '3dstrudel/bin/strudel_setChimeraX.py'],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "Operating System :: OS Independent"],
     install_requires=['scipy', 'biopython', 'numpy', 'mrcfile', 'mpi4py', 'psutil'],
     python_requires='>=3.7',
     license="Apache License"
 )