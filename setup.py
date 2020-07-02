
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
     name='strudel',
     version='0.1',
     author="Andrei Istrate",
     author_email="andrei@ebi.ac.uk",
     description="Strudel package",
     include_package_data=True,
     scripts=['strudel/bin/run_mapMotifValidation.py', 'strudel/bin/run_mapAveraging.py',
              'strudel/bin/run_penultimateClassifier.py', 'strudel/bin/run_chopModelMapMPI.py'],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: Apache",
         "Operating System :: OS Independent",
     ],
 )