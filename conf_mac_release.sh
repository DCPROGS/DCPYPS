#! /bin/bash

arch_flags=""

./configure CXX="/usr/bin/g++-4.2" CC="/usr/bin/gcc-4.2" LD="/usr/bin/g++-4.2" CFLAGS="$arch_flags" CXXFLAGS="$arch_flags -I/opt/local/include" LDFLAGS="$arch_flags -headerpad_max_install_names -L/opt/local/lib -L/usr/lib" --with-blitz=/opt/local PYTHON=/opt/local/bin/python
