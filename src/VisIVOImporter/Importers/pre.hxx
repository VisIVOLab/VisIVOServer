// file      : xsd/cxx/pre.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2007 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file


#ifdef _MSC_VER
#  if (_MSC_VER >= 1400)
#    include <xsd/cxx/compilers/vc-8/pre.hxx>
#  elif (_MSC_VER >= 1300)
#    include <xsd/cxx/compilers/vc-7/pre.hxx>
#  else
#    error Microsoft Visual C++ 6 is not supported.
#  endif
#endif
