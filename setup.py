import os
import sys
from setuptools import setup
from codecs import open

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding="utf-8") as f:
        return f.read()

setup(
    name="CompGIT",
    version="0.1",
    description="Library for Geometric Invariant Theory (GIT) in Sagemath. Its algorithms describe the unstable, non-stable and stictly polystable loci of n-dimensional projective space acted upon by a simple group. Based on P. Gallardo, J. Martinez-Garcia, H.-B. Moon, D. Swinarski. Computation of GIT quotients of semisimple groups.",
    long_description=readfile("README.md"), # get the long description from the README
    url="https://github.com/Robbie-H/CompGIT",
    author="Patricio Gallardo, Jesus Martinez-Garcia, Han-Bom Moon, David Swinarski",
    author_email="jesus.martinez-garcia@essex.ac.uk",
    license="License :: GPL-3.0 license",
    classifiers=[
      "Development Status :: 5 - Production/Stable",
      "Intended Audience :: Science/Research",
      "Topic :: Mathematics",
      "License :: GNU GENERAL PUBLIC LICENSE VERSION 3",
      "Programming Language :: Python :: 3"
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = "SageMath geometry",
    packages = ["CompGIT"],
)
