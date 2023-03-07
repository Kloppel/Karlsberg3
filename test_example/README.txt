Karlsberg2+ (version 2.0)
Created by: Tim Meyer and Jovan Dragelj 
Property of Freie Universität Berlin


Instructions and Example

*Instructions:

Full documentation of Karlsebrg2+ is not yet available.
For the most relevant modules, functions are described within the Python script.
In the kbp2 module, general overview of Python modules that belong to Karlsberg2+ can be found.
for more information about Karlsberg2+ background, please consult the following literature:

1. T. Meyer, E.W. Knapp, pKa Values in Proteins Determined by Electrostatics Applied to Molecular Dynamics Trajectories, J. Chem. Theory Comput. 11 (2015) 2827–2840.
2. G. Kieseritzky, E.W. Knapp, Optimizing pKA computation in proteins with pH adapted conformations, Proteins Struct. Funct. Genet. 71 (2008) 1335–1348.
3. B. Rabenstein, E.W. Knapp, Calculated pH-Dependent Population and Protonation of Carbon-Monoxy-Myoglobin Conformers, Biophys. J. 80 (2001) 1141–1150.

If Karlsberg2+ is used to generate data that will be publicly presented or published, please give credits to authors by citation. 
Thank you!

*Example:

Karlsberg2+ example script is "get_pkas_example.py", made for titration of the protein with the PDB code 4pti. Crystal structure coordinates are provided in the file 4pti.pdb1.

More details on the required parameters and settings are given within the script as comments after the "#" usually above the line.

IMPORTANT NOTE: the script "get_pkas_example.py" requires editing by the user, mostly regarsing filepaths and work folders. Please do that before usage in any case (be it OPTION 1 or OPTION 2).

Script can be executed as described in OPTION 1 or OPTION 2 in INSTALL_README.txt.
In case of an OPTION 1, from an interpreter, like Pycharm, the execution is strightforward.
In case of an OPTION 2, from the terminal, "run_pka.sh" script can be utilizied which takes care of previously described PYTHONPATH environment settings. 





