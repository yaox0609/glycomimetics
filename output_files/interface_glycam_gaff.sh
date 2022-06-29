#!/bin/bash
glycomimetics_output_dir="$1"
amber_gaff_dat="/home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/glycomimetics/output_files/gaff.dat"

cd ${glycomimetics_output_dir}
for i in $(ls -d */); do
    echo ${i}
    cd ${i}

    moiety_name=$(echo ${i} | sed 's/\///')
    echo "Moiety name: ${moiety_name}"

    if [[ ${i} == "natural_ligand"* ]]; then
        ligand_pdb="natural_ligand.pdb"
    else
        ligand_pdb="${moiety_name}_ligand.pdb"
    fi


    mol2=corona.mol2
    pdb2glycam_log=${moiety_name}_pdb2glycam.log
    antechamber_frcmod=corona.frcmod
    output_glycam_gaff_frcmod=${moiety_name}_glycam_gaff.frcmod
    output_glycam_gaff_off=${moiety_name}_glycam_gaff.off

    #echo "/home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/glycomimetics/glycam_gaff_interfacing/main.exe ${mol2} ${ligand_pdb} ${pdb2glycam_log} ${amber_gaff_dat} ${antechamber_frcmod} ${output_glycam_gaff_frcmod} ${output_glycam_gaff_off}"
    /home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/glycomimetics/glycam_gaff_interfacing/main.exe ${mol2} ${ligand_pdb} ${pdb2glycam_log} ${amber_gaff_dat} ${antechamber_frcmod} ${output_glycam_gaff_frcmod} ${output_glycam_gaff_off}
    cd ..
done

