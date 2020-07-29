#!/bin/bash
# this file takes one parameter: the file name of the list of vdbs that will be surveyed 
toBeSurveyed='listOfVDBs.list' # [1] 
# for each file, the filename (sans extension) is the pdb code
while IFS= read -r lineA 
do 
    PoA=( $lineA ) # points of attachment
	PDB=`echo "${PoA[-1]}"| cut -d'.' -f1` 
    unset PoA[-1]
	
	# check if the PDB code + "spacefilling.pdb" file exists.
	# if it does, add a number, or add onto the number for new file name
    for i in "${PoA[@]}"
    do
        name=${PDB}_$i'_fsp'
        if [[ -e spacefills/${PDB}_${PoA}_fsp.pdb || -L spacefills/${PDB}_${PoA}_fsp.pdb ]]
        then
            j=2
            while [[ -e spacefills/${PDB}_${PoA}_fsp-${j}.pdb || -L spacefills/${PDB}_${PoA}_fsp-${j}.pdb ]]
            do
                let j++
            done
            name=`echo "${name%%-*}"`'-'$j 

        fi
        # open pymol
        /mnt/c/Users/janis/PyMOL/PyMOLWin.exe spacesurvey.py -- $PDB $i $name # use option -c to suppress PyMOL GUI
    done

done < $toBeSurveyed
