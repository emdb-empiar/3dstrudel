
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
     name='threed_strudel',
     version='0.9.16',
     author="Andrei Istrate",
     author_email="andrei@ebi.ac.uk",
     description="Strudel package",
     include_package_data=True,
     scripts=['threed_strudel/bin/strudel_mapMotifValidation.py', 'threed_strudel/bin/strudel_mapAveraging.py',
              'threed_strudel/bin/strudel_penultimateClassifier.py', 'threed_strudel/bin/strudel_chopModelMapMPI.py',
              'threed_strudel/bin/strudel_setChimeraX.py', 'threed_strudel/bin/strudel_sequenceFinder.py',
              'threed_strudel/bin/strudel_runValidationBatch.py', 'threed_strudel/bin/strudel_ultimateClassifier.py'],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="",
     packages=setuptools.find_packages(),
     install_requires=['scipy', 'biopython', 'numpy', 'mrcfile', 'psutil'],
     python_requires='>=2.6',
     license="Apache License"
 )