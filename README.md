IPS Project on the local density in a nuclear system <br/>
To see more details about the problem and its solution, go to see the documentation

Prerequisite:<br/>
-armadillo version>=8.200.2<br/>
-cxxtest <br/>
-gnuplot <br/>
-doxygen <br/>

Table of contents
==================

1. [Compilation](#Compilation) 
2. [Documentation](#Documenation)
3. [Run](#Run)
4. [Plot](#Plot)
5. [Test](#Test)


# Compilation <a name="Compilation"></a>
Execute the makefile at the root of the project, to compile source code and unitary test.
The executable of the application will be in the folder src, whereas the executable for test is in test.
There is also a target clean, to clean folders src and test.


# Documentation <a name="Documentation"></a>
The makefile at the root have a target doc, so you can type "make doc" to generate all the documentation doxygen in the folder html 


# Run <a name="Run"></a>
In the folder src, you can execute the executable file "main" with 2 parameters:<br/>
   1) s_z=size of the vector Z <br/>
   2) s_r=size of the vector R <br/>
It will generate 1 new file "result.dat" which contains value of the result, in a format to plot with gnuplot


# Plot <a name="Plot"></a>
In the folders src, you can type "make plot" to generate "reslut.png" which represents solution associated to your last run. <br/>

# Test <a name="Test"></a>
The executable to run test is in the folder test, and you just need to type "./test" to launch existing tests
