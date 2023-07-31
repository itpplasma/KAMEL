# KiLCA_interface 
This is the KiLCA interface in python written to easily communicated with KiLCA. KiLCA takes confiuration and input files in the form of plain text files. For convenience, the interface automates the reading and writing of these files, while changing an input variable is simply done by chaging a class variable in the interface. Tutorials in the /tutorial directory show hot this is done.

# KiLCA_postprocessor
This is a postprocessing class that handles the output of KiLCA. Most importantly, it handles the EB.dat file that contains the electromagnetic fields. For what this file contains, look in ../../Documentation/kilca-2.pdf. Further, this class includes plotting mehtods.
