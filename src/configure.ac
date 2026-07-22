#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TelSeq], [0.3.0], [https://github.com/michtrofimov/telseq-parallel/issues])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Telseq/telseq.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for header files.
AC_LANG([C++])

AC_ARG_WITH(bamtools,
        AS_HELP_STRING([--with-bamtools=DIR], [specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]),
        [ if test -d "$withval" ; then
            bamtools_ldflags="-Wl,-rpath,$withval/lib -L$withval/lib -Wl,-rpath,$withval/lib/bamtools -L$withval/lib/bamtools"
    		bamtools_include="-I$withval/include -I$withval/include/bamtools"
		else
          	AC_MSG_ERROR([Cannot stat directory $withval])
      	fi ]
        )

AC_ARG_WITH(htslib,
        AS_HELP_STRING([--with-htslib=DIR], [specify directory containing HTSlib (https://github.com/samtools/htslib)]),
        [ if test -d "$withval" ; then
            htslib_ldflags="-Wl,-rpath,$withval/lib -L$withval/lib"
            htslib_include="-I$withval/include"
        else
            AC_MSG_ERROR([Cannot stat directory $withval])
        fi ]
        )

# Set compiler flags.
# Keep normal floating-point semantics for output compatibility. In
# particular, do not add -ffast-math here.
AC_SUBST(AM_CXXFLAGS, "-std=c++11 -Wall -Wextra -Wno-unknown-pragmas -g -O3 -pthread")
#AC_SUBST(CXXFLAGS, "-O3")
#AC_SUBST(CFLAGS, "-O3")

AC_SUBST(CXXFLAGS, "-std=c++11")
AC_SUBST(CFLAGS, "-std=c++11")

AC_SUBST(CPPFLAGS, "$CPPFLAGS $bamtools_include $htslib_include $gzstream_include")
AC_SUBST(LDFLAGS, "$bamtools_ldflags $htslib_ldflags $gzstream_ldflags $LDFLAGS")

#AC_MSG_NOTICE(["CPPFLAGS is $CPPFLAGS"])
#AC_MSG_NOTICE(["LDFLAGS is $LDFLAGS"])

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS -lbamtools -lhts -lz -pthread")

#AC_MSG_NOTICE(["LIBS is $LIBS"])

# Make sure the bamtools headers can be found
AC_CHECK_HEADER([api/BamReader.h],,[AC_MSG_ERROR(
[The bamtools library must be installed (http://github.com/pezmaster31/bamtools). \
You can specify its path with the --with-bamtools=PATH option])])

AC_CHECK_HEADER([htslib/sam.h],,[AC_MSG_ERROR(
[HTSlib must be installed (https://github.com/samtools/htslib). \
You can specify its path with the --with-htslib=PATH option])])

AC_CONFIG_FILES([Makefile
    Util/Makefile
    Telseq/Makefile
    Test/Makefile])

AC_OUTPUT
