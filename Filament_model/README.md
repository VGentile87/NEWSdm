CREATED BY V.GENTILE (2020/05/27)

The project makes a simulation of AgBr filaments produced by ionizing particles in NIT emulsions.

The project contains several files:
1) coord_crystal.txt  --> NIT 70nm crystal framework 1um3
2) C30V_segmentedfilaments.txt  -->  Data from SEM of C30V ions in NIT
3) SRIM_data.h  --> Example of header file to read SRIM data
4) nkinks_weights.txt  -->  Table of weigths to assign number of kinks according to the filament length
5) filament_model.C  --> Main script for the filament simulation. It gives in output the file "filament.root" with all MC info.
6) drawEvent.C  --> It reads data from "filament.root" file and draw an event given as input.

Usage:

root -l
.L filament_model.C
arun()

-----------------

root -l drawEvent.C (and type the event ID)
