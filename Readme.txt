XLPack SDK (for XLPack 7.0) (September 5, 2024)

1. What is "XLPack SDK"

XLPack SDK includes the library software and the sample codes for using the
numerical computation module of XLPack (and XLPack Lite) from development
environments other than Excel VBA such as C/C++. The program developed by using
the SDK can run on the personal computer in which XLPack is installed.

Note - The sample codes are covered by the basic feature set of XLPack without
addons, and can be used if XLPack Basic is installed.

2. Terms of use of this software

Please follow the following conditions when using.

- All the rights to this product are the sole property of K Technologies.
- This software provides the internal interface in as-is condition. It has not
  been thoroughly tested under all conditions. Therefore, the function,
  performance, or reliability of this software are not guaranteed. The
  specifications are subject to change without notice.
- This software must be used under the responsibility of the user. There is no
  warranty as to any damage caused as a result of using the software. Support
  services are not provided for this software.
- Redistribution of this software in any format is prohibited.

3. Contents of this software

The following files are included in this SDK.

Readme.txt - Read me file (this file)
Readme_ja.txt - Read me file (Japanese)
include\
  cnumlib.h - C interface header
  cnumlib_mangling.h - C interface compatibility header (see item 6)
  cnumlib_complex.h - C interface header file (complex definition)
  cblas.h - CBLAS header
  lapacke.h - LAPACKE header
  spcnumlib.h - C interface header (sparse matrices)
  spcnumlib_mangling.h - C interface compatibility header
  pdecnumlib.h - C interface header (PDE)
  pdecnumlib_mangling.h - C interface compatibility header
lib\
  XLPack.lib - C interface library (64 bit version)
  XLPack_32.lib - C interface library (32 bit version)
  Lapacke.lib - LAPACKE/CBLAS library (64 bit version)
  Lapacke_32.lib - LAPACKE/CBLAS library (32 bit version)
samples\
  Sample codes: txt files are the command line example output.
    C_C++ - Sample code in C/C++
      To run testcpp_Matrix, Matrix.h and MatrixIO.h, which are described in
      Stroustrup "Programming -- Principles and Practice Using C++, 1st Ed."
      Addison-Wesley, 2009, are required.
    Python - Sample code in Python
      numpy is required.
      pyd file is the 64 bit binary module for Windows versions for Python 3.7
      or later. Use XLPack.py (ctypes version with equivalent functions) for 32
      bit Python and if the binary module does not work well.
    C# - Sample code in C#

4. System environment

This software has been tested by using the following environment. The
modification of sample codes may be required depending on the installation
conditions.

- XLPack 7.0.0
- Windows 11 (23H2)
  - Visual Studio 2022 (17.11.0)
  - Python 3.12.5 (numpy 2.0.1)

5. Documents

Online manuals are available at our web site(*).

(*) https://www.ktech.biz


6. Change of C interface

In this version, an underscore (_) is added to the beginning of C/C++ function
names in XLPack 6.0 and former. cnumlib_mangling.h header defines the
compatible function names.

---
(C) 2014-2024  K Technologies
