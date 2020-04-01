#!/bin/bash

# Rui's copy

# parameters to be changed by user:

structure='<folders>/<vlp file>.vdb'
HETATM=0 # set to 0 for false and 1 for true to include HETATMs like glycans, cofactors, solvent

# there should be a distinction made for finding N term in VDB files and in 
# PDB files
# in PDB files, the next N-term immediately follows the "TER" term, whereas in 
# VDB files, there are matrix headers before the next N-term. 
# the above was just a generalization, but to the best of my knowledge, that 
# is how the files are set up

# because of this distinction, this program should be able to figure out the 
# file extension and act dependently on that. In addition, we must be more 
# careful about assigning certain file extensions to certain files and not 
# others.

bash findCoordinates.sh $structure $HETATM;

if [[ -e <folders>/blacklist.list ]]; then rm <folders>/blacklist.list; fi
if [[ -e <folders>/orgCoord.list ]]; then rm <folders>/orgCoord.list; fi

# determine if the file is a pdb or vdb
PDB=`echo "${structure##*/}"| cut -d'.' -f1`
extension=`echo "${structure##*.}"`

if [[ $extension == "vdb" ]]
then 
	# find the first instance of the string "TER"
	TERline=`grep -n -m 1 "TER" "$structure" | cut -d':' -f1`
	# line number that the first "TER" shows up
	TERline=$((TERline-1))
	# write all the lines from the ogMatrix to file
	`head -n +$TERline "$structure" | sed -n '/ATOM/p' > ogMatrix.txt` 
	# find first instance lines of N-termini
	`sort -u -k 5,5 ogMatrix.txt > ogMatrix2.txt` 
	
else
	# find the first line with the tag "ATOM" and store the line
	`grep -m 1 "ATOM      1" "$structure" > ogMatrix.txt`
	# find all the line numbers that the string "TER" shows up
	TERline=($(grep -n "TER   " "$structure" | cut -d':' -f1))
	for ((i = 0; i != ${#TERline[@]}; i++)); 
	do
		TERline[$i]=$(( ${TERline[i]} + 1 ))
		# write the line after the "TER" to file
		`tail -n +${TERline[i]} "$structure" | head -n 1 >> ogMatrix.txt`
		
	done
	# find the unique chains and keep
	cat ogMatrix.txt | sed -n '/ATOM/p' | sort -u -k 5,5 > ogMatrix2.txt

fi

resi=()
while IFS= read -r lineA
do
	# keep the formatting of the N-termini and copy only the chain and 
	# residue numbers
	printf '%s\n' "$lineA" | cut -c21-26 >> NsearchTerms.list
	chain=`printf '%s' "$lineA" | cut -c22`
	number=`printf '%s' "$lineA" | cut -c23-26 | tr -d '[:space:]'`
	residues+=("$number$chain")
done < ogMatrix2.txt
echo ${residues[@]} > residues.list

rm ogMatrix2.txt



echo "$structure" >> <folders>/orgCoord.list
while IFS= read -r lineC
do
 	if [[ "$lineC" != "" ]]
	then 
		echo "$lineC" >> OsearchTerms.list
		searchTerm=("$(grep -m 20 "$lineC" "$structure")")
		newOrg=`printf '%s\n' "${searchTerm[@]}" | grep -m 1 "  CA  "`
		if [[ "$newOrg" != "" ]]
		then 
	 		orgCoord=(`echo "$newOrg" | cut -c31-38 | tr -d '[:space:]'`)
	 		orgCoord+=(`echo "$newOrg" | cut -c39-46 | tr -d '[:space:]'`)
	 		orgCoord+=(`echo "$newOrg" | cut -c47-54 | tr -d '[:space:]'`)
			echo "${orgCoord[@]}" >> <folders>/orgCoord.list
		fi
		
	fi	
done < <folders>/NsearchTerms.list
echo -ne "\n" >> <folders>/orgCoord.list

rm ogMatrix.txt
rm NsearchTerms.list

while IFS= read -r lineE
do
	echo "$(grep -m 100 "$lineE" "$structure")" > bl1.list
		sort -u -k 2,2 bl1.list > bl2.list 
		rm bl1.list
		while IFS= read -r lineD
		do
			adjacent=(`echo "$lineD" | cut -c31-38 | tr -d '[:space:]'`)
			adjacent+=(`echo "$lineD" | cut -c39-46 | tr -d '[:space:]'`)
            adjacent+=(`echo "$lineD" | cut -c47-54 | tr -d '[:space:]'`)
            echo ${adjacent[@]} >> blacklist.list
		done < bl2.list
		rm bl2.list

done < <folders>/OsearchTerms.list

rm OsearchTerms.list


if [[ ! -e figures/ ]]; then mkdir figures/; fi
python occlusion_map.py "$PDB"

