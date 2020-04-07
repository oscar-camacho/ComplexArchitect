# Setup.py python distribution script to install Complexbuilder

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='comparch',
      version='1.0',
      description='Application designed to generate macrocomplex structures from simple pair interactions',
      data_files = [("", ['LICENSE.txt'])],
      license = "MIT",
      scripts=['comparch/carchitect'],
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='macrocomplex structural bioinformatics',
      url='https://github.com/oscar-camacho/SBI-PYT_PROJECT',
      author='Oscar Camacho',
      author_email='osc1997@gmail.com',
      packages=find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
        ],
      install_requires=['biopython', 'matplotlib', 'numpy'],
      include_package_data=True)
