# TiFoSi
TiFoSi (Tissues: Forces and Signaling), pronounced [tiˈfoːzi], is an efficient simulation code of epithelia (multiple cell populations) that allows to include feedback between tissue growth/mechanics and gene regulation.

TiFoSi simulations rely on a C++ code. Yet, the user doesn't need to have any knowledge of this programming language. Thus, the characteristics and properties of the tissue to simulate are prescribed in a XML configuration file (config.xml).

System Requirements

TiFoSi should work in any system with a standard installation of the following required packages. The code has been tested extensively in Linux (Ubuntu) systems but in principle should run also on any other OS assuming the following packages are installed:

• g++ (version 4.2 or above) compiler: http://gcc.gnu.org

• python (version 2.7.2 or above): http://www.python.org

• make

We note that g++ and python must be in the PATH.

Installation

After downloading, TiFoSi comprises the following files and folders:

• Folders: bin, pylib, pylib3, src, and templates

• Files: config.xml, compile.py, compile3.py, buildandrun, buildandrun3, and makePlot.py

pylib, pylib3, src, and templates folders contain the source files of the code and the directory bin is where, after a successful compilation, the executable of the code (tifosi) will be copied. The files at directories pylib, pylib3, src, and templates, with the exception of the Makefile file contained in the directory src, must NOT be edited by the user unless you really know what you´re doing.

Before using the code for the first time, edit the Makefile file and indicate the flag for the variable SYS depending on your system (UNIX, MAC or WIN). In all cases edit the files with a non-formatting editor, e.g. emacs.

Configuration/Execution

A TiFoSi simulation relies on a configuration file, config.xml, where the properties of the simulation are specified. Edit that file following the rules indicated in the TiFoSi manual. Once the configuration file has been defined, execute the following command in a terminal window:

• python compile.py (in case you are using python2.7)

• python compile3.py (in case you are using python3.x)

This command will first check that the configuration file has been properly defined. Afterwards it will generate/modify the C^{++} source files and finally it will compile the code and copy the executable, tifosi, in the bin directory. If all these processes have been successful, the last line of the standard output of the command will read,

The process has been successfully completed!

This means that you are ready to go! Just execute tifosi (inside the bin folder) and wait for the simulation to finish. The result of a simulation (output files) will be saved in the directory where the executable tifosi is run. For example, in a Linux system, and assuming that you are at the directory where the executable file, tifosi, is,

./tifosi

When the simulation finishes the following message is shown at the standard output (console window),

***********************Smile! the simulation is over!***********************

Automatization of the compilation and execution process

The script files "buildandrun" and "buildandrun3" (UNIX) allows to automatization the compilation and execution in a single step. Basically this file combines an initial check of the requirements to use TiFoSi, a call to the python file processing script, and running the executable file together with other features: the possibility of creating a directory for the executable, an error detection step, the possibility to visualize the configuration file before compilation...
In order to use this file the user must first change its attribute an make it executable:

• chmod +x ./buildandrun

• chmod +x ./buildandrun3

After that the script can be used by single executing:

• ./buildandrun (in case you are using python2.7)

• ./buildandrun3 (in case you are using python3.x)


Output Visualization

The python (3.x) script makePlot.py generates images to visualize the data files generated during a simulation (protein concentration and cell types). Alternative visualization tools can be found in the webpage of the project http://tifosi.thesimbiosys.com (visualization tools)

makePlot.py requires matplotlib:

python -m pip install -U pip
python -m pip install -U matplotlib

In case pip install fails:

sudo apt-get install python3-pip (python3.x)

Once matplotlib is installed, to execute the makePlot.py script:

python makePlot.py --ifolder INPUT_FOLDER --ofolder OUTPUT_FOLDER

The above command will generate PNG images for all frames. As a preview of a simulation the following command will generate PNG images for the final frame of each simulation stage (all proteins and cell types):

python makePlot.py --ifolder INPUT_FOLDER --ofolder OUTPUT_FOLDER --preview

Further options can be passed in the command line to makePlot.py. For a full list of options type:

python makePlot.py --help


