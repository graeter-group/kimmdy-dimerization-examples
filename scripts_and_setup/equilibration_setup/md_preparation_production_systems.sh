#!/bin/bash

cd ../PDBs
ls
run_num=3

for file in *.pdb; do
	current_folder="${file%.pdb}"
	mkdir -p "../$current_folder"
	
	echo "Working on simulation based on ${file%.pdb}"

	sequence_name="${file%.pdb}"
	for (( run=1; run<=run_num; run++ )); do
		mkdir -p "../$current_folder/R$run"
		cd "../$current_folder/R$run" || continue 
		cp "../../PDBs/$file" .

	
		# Copy input files (mdps, resubmit scripts_and_setup, forcefield) and name jobs
		cp -r ../../InputFiles/* .
		sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${sequence_name}_R${run}_NPT|" resubmit_npt.sh
        	sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${sequence_name}_R${run}_NVT|" resubmit_nvt.sh
	
		# Create gromacs files
		echo -e "1\n1\n" | GMXLIB="$(pwd)" gmx pdb2gmx -f "${sequence_name}.pdb" -o "${sequence_name}.gro" -p "${sequence_name}.top" -i "${sequence_name}.itp" -merge "all" 

		# Adjust box size (HARDCODED, depends on system see SI)
		gmx editconf -f "${sequence_name}.gro" -o "${sequence_name}_newbox.gro" -c -box "18" "18" "18"

		# Solvate
		gmx solvate -cp "${sequence_name}_newbox.gro" -cs "spc216.gro" -o "${sequence_name}_solvated.gro" -p "${sequence_name}.top"

		# Add ions 
		gmx grompp -f "ions.mdp" -c "${sequence_name}_solvated.gro" -p "${sequence_name}.top" -o "ions.tpr"
		echo -e "3\n" | gmx genion -s "ions.tpr" -o "${sequence_name}_ions.gro" -p "${sequence_name}.top" -pname "NA" -nname "CL" -neutral

		# Energy minimization
		gmx grompp -f "minim.mdp" -c "${sequence_name}_ions.gro" -p "${sequence_name}.top" -o "${sequence_name}_EM.tpr" -maxwarn 1
		gmx mdrun -v -deffnm "${sequence_name}_EM"
	
		# Renaming
		mv "${sequence_name}.top" "${sequence_name}_EM.top"
		rm ./"#"*
		rm ./step*

        	cd "../../PDBs/"
	done
done
