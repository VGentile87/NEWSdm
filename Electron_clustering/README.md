CREATED BY V.GENTILE (2020/02/19)

The script performs two times a density based clustering of all energy deposit due to electrons stopped in nuclear emulsions.
The first iteration is meant to form the small clusters which represent the sensitize grains; the second one to produce the big
clusters which represent the visible grains as seen by optical microscopes.
The script  provides a TTree with position and total energy deposit of big clusters. 
A threshold for the visibility of small clusters is required.

Usage:

root -l 
.L Make_Tree_test.C+
myrun()
