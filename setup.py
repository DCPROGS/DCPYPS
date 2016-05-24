#!/usr/bin/env python
"""Notes
-----
DCPYPS is a collection of modules and scripts common to all the DCPROGS repositories. 
DCPYPS deals with file input/output, specifying kinetic mechanism, creating data structures. 
SCALCS and HJCFIT depend on DCPYPS.
"""
# Modified from numpy's setup.py

DOCLINES = __doc__.split("\n")

import os
import sys

CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU General Public License (GPL)
Programming Language :: Python
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

NAME                = 'DCPYPS'
AUTHOR              = "Remigijus Lape, Christoph Schmidt-Hieber et al."
AUTHOR_EMAIL        = "r.lape@ucl.ac.uk"
# MAINTAINER          = "DC_PyPs Developers"
# MAINTAINER_EMAIL    = "c.schmidt-hieber@ucl.ac.uk"
DESCRIPTION         = DOCLINES[0]
LONG_DESCRIPTION    = "\n".join(DOCLINES[2:])
PACKAGES            = ["dcpyps", "dcpyps.samples", "dcpyps.sccalc"]
SCRIPTS             = []
URL                 = "https://github.com/DCPROGS/DCPYPS"
DOWNLOAD_URL        = "https://github.com/dcprogs/dcpyps/tarball/master"
LICENSE             = 'GPL2'
CLASSIFIERS         = filter(None, CLASSIFIERS.split('\n'))
PLATFORMS           = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"]
MAJOR               = 0
MINOR               = 2
MICRO               = 0
ISRELEASED          = True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def write_version_py(filename='dcpyps/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM DCPYPS SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION = VERSION

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def setup_package():

    # Perform 2to3 if needed
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path

    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    write_version_py()

    # Run build
    from distutils.core import setup

    try:
        setup(
            name=NAME,
            version=VERSION,
            packages=PACKAGES,
            scripts=SCRIPTS,
            author=AUTHOR,
            author_email=AUTHOR_EMAIL,
            # maintainer=MAINTAINER,
            # maintainer_email=MAINTAINER_EMAIL,
            description=DESCRIPTION,
            long_description=LONG_DESCRIPTION,
            url=URL,
            download_url=DOWNLOAD_URL,
            license=LICENSE,
            classifiers=CLASSIFIERS,
            platforms=PLATFORMS)
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return

if __name__ == '__main__':
    setup_package()
