# CMOST_v3

This is a repository for CMOST: an open-source framework for the microsimulation of colorectal cancer screening strategies.

It is a newest implementation of the model described here:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5460500/

Full source code is available in "Source" folder. 
Compilation requires Boost library.

In order to compile for Windows machines one can use:
cl.exe /c midaco.c
cl.exe /Zi /EHsc /openmp /O2 /Ot /Oi "/IC:PATH_TO_BOOST_HEADERS" SimulationParameters.cpp Polyp.cpp Cancer.cpp Person.cpp Output.cpp Screening.cpp Optimizer.cpp Evaluate.cpp Stratification.cpp main.cpp midaco.obj /link /out:simulateCMOST.exe "/LIBPATH:PATH_TO_BOOST_LIBRARIES"

Simulation settings are defined in .ini files (in main directory and "InputFiles" folder).
Attached is document explaining main parameters in settings.ini file.

You can find compiled executables for Windows and Linux in the main folder (simulateCMOST files).

There is also attached exemplary script in R language showing how to plot some of the simulation results.


-Jan Poleszczuk
