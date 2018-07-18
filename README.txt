Code for simulating the LOLX prototype detector.

Compiling Code:

Before attempting to compile the code make sure you have the following codes:

GEANT4 (https://github.com/Geant4/geant4.git)
ROOT (http://github.com/root-project/root.git)
CADMesh (https://github.com/christopherpoole/CADMesh.git)

Once dependancies are installed, build data structure library

$cd ds
$make

Build the simulation in a new build directory using cmake:

$mkdir build
$cmake ../
$cmake --build ./

To run the simulation:

$cd build
$./LXe ../runSr90.mac output.root

To process output file into useful histograms:

$cd analysiscode
$./ProcessRun In.root Out.root 0,1 (0 = Cylindrical, 1 = Sphere)
