Version  1.0
------------

This work was funded by subcontract NIHGMS grant: GM076692 
from James Glazier, Indiana University.

Authors: Deepak Chandran and Herbert M Sauro

==================
     ABOUT
==================

This is a simplified C API to the Copasi C++ library. The API uses 
a few data structures that wrap Copasi's C++ classes as well as 
provide some additional features, such as being able to use simple
names, e.g. "A", instead of Copasi's complete names, e.g. <model1.compartmentX.A>

==================
  DEPENDENCIES
==================

lapack -- linear algebra package, used by Copasi C++ code
libSBML -- Systems Biology Markup Language parser, used by Copasi C++ code
expat -- XML parser used by libSBML and Copasi C++
curl -- required by Raptor library
raptor -- RDF parser required by Copasi C++
muparser -- used within the C API for parsing mathematical expressions and getting the tokens
boost -- the boost regex library is used by the C API to substitute strings inside formulas

===============================
  GETTING THE SOURCE CODE
===============================
First install Subversion. In Linux and Mac, just use apt-get svn.
In Windows, download TortoiseSVN (graphical version) or SilkSVN (command-line version)

Once Subversion is installed, use the following command to get the source:
svn checkout https://copasi-simple-api.googlecode.com/svn/trunk/ copasi-simple-api

If you are using TortoiseSVN, first create an empty folder called copasi-api. Then
right-click on that newly created folder and select "SVN Checkout". Enter 
https://copasi-simple-api.googlecode.com/svn/trunk/ as the source URL

===================================================
  BUILDING THE LIBRARY FOR WINDOWS, MAC AND LINUX
===================================================

Lets assume that the source code is located in a folder named "copasi-simple-api"

1) Install Cmake from cmake.org
    for Linux, just type 'apt-get cmake cmake-gui'

2) Run cmake-gui
    In the Cmake program, find two text boxes labeled, "where is the source code" and "where to build"

3) Type the full copasi-simple-api/ folder path for the line that says "where is the source code"

4) Type copasi-simple-api/BUILD folder path for the line that says "where to build"

5) Click the "Configure" button

6) Cmake will create the BUILD folder and then ask you to select the compiler. 
    It is recommended that you use standard compilers, e.g. gcc for Unix and 
    Visual Studio or MinGW for MS Windows.

7) Once configuration is done, click the "Generate" button

8) Go to the copasi-simple-api/BUILD folder. 
    If you used Visual Studio as the compiler, you will find a project file (.sln) -- open it and build all
    If you used GCC, you will find a makefile. cd into this folder and run "make"

==========================
  USING THE LIBRARY
==========================

See test_copasi_api.c for example





