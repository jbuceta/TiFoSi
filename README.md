# TiFoSi
TiFoSi (Tissues: Forces and Signaling), pronounced [tiˈfoːzi], is an efficient simulation code of epithelia (multiple cell populations) that allows to include feedback between tissue growth/mechanics and gene regulation.

TiFoSi simulations rely on a C++ code. Yet, the user doesn't need to have any knowledge of this programming language. Thus, the characteristics and properties of the tissue to simulate are prescribed in a XML configuration file (config.xml).

System Requirements

TiFoSi should work in any system with a standard installation of the following required packages. The code has been tested in Mac OS X and Linux (Ubuntu) systems but in principle should run also on Microsoft Windows systems or any other OS, e.g. Solaris.:

• gcc (version 4.2 or above) compiler: http://gcc.gnu.org

• python (version 2.7.2 or above): http://www.python.org

These packages are (usually!) installed by default in standard Linux distributions. We note that gcc and python must be in the PATH.

Installation

TiFoSi is distributed in a single ZIP compressed file; its name indicates the version of the code, e.g. V3.5.zip. Copy this file in any directory, e.g. TIFOSI, and decompress. The ZIP file contains the following files and subdirectories:

• Directories: bin, pylib, src, and templates

• Files: config.xml and compile.py

pylib, src, and templates directories contain the source files of the code and the directory bin is where, after a successful compilation, the executable of the code (tifosi) will be copied. The files at directories pylib, src, and templates, with the exception of the Makefile file contained in the directory src, must NOT be edited by the user unless you really know what you´re doing.

Before using the code for the first time, edit the Makefile file and indicate the flag for the variable SYS depending on your system (UNIX, MAC or WIN)In all cases edit the files with a not formatting editor, e.g. emacs.. For example, in the case of a Linux or other UNIX system (but not Mac OS X)Note that in the case of a Windows system it is also necessary to specify the complete path to the gcc compiler executable.:

# Set this variable to either UNIX, MAC or WIN

SYS = UNIX

THE COMPILER PATH IS ALSO REQUIRED!!! (NEXT LINES IN THE MAKEFILE)

Configuration/Execution

A TiFoSi simulation relies on a configuration file, config.xml, where the properties of the simulation are specified. Edit that file following the rules given below. Once the configuration file has been defined, execute the following command in a terminal window:

python compile.py_{{\color{red}\hookleftarrow}} 

This command will first check that the configuration file has been properly defined according to the rules given below. Afterwards it will generate/modify the C^{++} source files and finally it will compile the code and copy the executable, tifosi, in the bin directory. If all these processes have been successful, the last line of the standard output of the command will read,

The process has been successfully completed!

This means that you are ready to go! Just execute tifosi and wait for the simulation to finishThe result of a simulation (output files) will be saved in the directory where the executable tifosi is run.. For example, in a Linux system, and assuming that you are at the directory where the file tifosi is,

./tifosi_{{\color{red}\hookleftarrow}} 

When the simulation finishes the following message is shown at the standard output,

***********************Smile! the simulation is over!***********************

