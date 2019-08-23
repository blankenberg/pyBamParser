import sys

if sys.version_info < (2, 5):
    print("ERROR: pyBamParser requires python 2.5 or greater", file=sys.stderr)
    sys.exit()

try:
    import setuptools
except ImportError:
    # Automatically download setuptools if not available
    from distribute_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages
from glob import glob

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = True

       
def main():
    setup(  name = "pyBamParser",
            version = "0.0.3",
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            scripts = glob( "scripts/*.py" ),
            author = "Daniel Blankenberg",
            author_email = "dan.blankenberg@gmail.com",
            description = "Tools for parsing BAM data",
            license = "GPLv2",
            url = "https://github.com/blankenberg/pyBamParser",
            zip_safe = False,
            dependency_links = [],
            classifiers=[ "Development Status :: 4 - Beta",
            "License :: OSI Approved :: GNU General Public License v2 (GPLv2)" ],
            **extra )

if __name__ == "__main__":
    main()
