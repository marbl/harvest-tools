import sys
import os
from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#ext_modules = [Extension("harvest-phylogeny", ["./src/harvest-phylogeny.py"]),Extension("harvest-snps",["./src/harvest-snps.py"]), Extension("harvest-alignments",["./src/harvest-alignments.py"]), Extension("harvest-backbone",["./src/harvest-backbone.py"]),Extension("harvest-variants",["./src/harvest-variants.py"])]

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "ParSNP",
    version = "1.0b",
    author = "Brian Ondov & Todd J. Treangen",
    author_email = "treangen+ParSNP@gmail.com",
    description = ("Ultrafast core genome alignment and WGT"),
    license = "Perl Artistic License",
    keywords = "multialignment snps bionformatics",
    url = "http://github.com/MGI/ParSNP/wiki",
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
    scripts=['src/harvest-phylogeny.py', 'src/harvest-snps.py', 'src/harvest-alignments.py', 'src/harvest-backbone.py', 'src/harvest-variants.py'],
#    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext}
)
