XLPack SDK (for XLPack 6.0) (July 11, 2022)

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

The following files are included in this ZIP file.

Readme.txt - Read me file (this file)
Readme_ja.txt - Read me file (Japanese)
include\
  cnumlib.h - C interface header file (C language)
  cnumlib - C interface header file (C++)
  cblas.h - CBLAS header file
  lapacke.h - LAPACKE header file
lib\
  XLPack.lib - C interface library (64 bit version)
  XLPack_32.lib - C interface library (32 bit version)
  Lapacke.lib - LAPACKE/CBLAS library (64 bit version)
  Lapacke_32.lib - LAPACKE/CBLAS library (32 bit version)
samples\
  Sample codes: txt files are the command line example output files.
    C_C++ - Sample code in C/C++
      To run testcpp_Matrix, Matrix.h and MatrixIO.h, which are described in
      Stroustrup "Programming -- Principles and Practice Using C++, 1st Ed."
      Addison-Wesley, 2009, are required.
    Python - Sample code in Python
      numpy is required.
      pyd files are the binary modules for Windows versions for Python 3.7 or
      leter. Binary modules may not work if the software environment is
      different from the environment below. In such cases, try XLPack.py (ctypes
      version with equivalent functions).
    C# - Sample code in C#
    VB - Sample code in VB.NET
    F# - Sample code in F#
    Julia - Sample code in Julia
    Pascal - Sample code in Pascal

4. System environment

This software has been tested by using the following environment. The
modification of sample codes may be required depending on the installation
conditions.

- XLPack 6.0.5
- Windows 10 (21H2)
  - Visual Studio 2022 (17.2.5)
  - Python 3.7.9 (numpy 1.21.6), 3.10.5 (numpy 1.22.4)
  - Julia 1.7.3
  - Free Pascal 3.2.2

4. Documents

Online manuals are available at our web site(*).

(*) https://www.ktech.biz

---
(C) 2014-2022  K Technologies
