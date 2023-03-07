#!/usr/bin/tcsh

echo 'Processing "combine_md_tools"'
echo 'Creating html file..'
cython -a analyse_md_pkas_tools.pyx
echo 'Compiling analyse_md_pkas_tools.pyx..'
python setup_analyse_md_pkas_tools.py build_ext --inplace


#echo 'Processing "combine_md_tools"'
#echo 'Creating html file..'
#cython -a combine_md_tools.pyx
#echo 'Compiling combine_md_tools.pyx..'
python setup_combine_md_tools.py build_ext --inplace

#echo 'Processing "setup_all_vs_all_atoms"'
#echo 'Creating html file..'
#cython -a all_vs_all_atoms.pyx
#echo 'Compiling combine_md_tools.pyx..'
#python setup_all_vs_all_atoms.py build_ext --inplace


