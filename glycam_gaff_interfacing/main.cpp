#include "../../gmml/includes/gmml.hpp"
#include "../../gmml/includes/MolecularModeling/assembly.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include "../../gmml/includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../gmml/includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../gmml/includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../gmml/includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../gmml/includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../gmml/includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../gmml/includes/utils.hpp"

#include "../src/vina_bond_by_distance_for_pdb.hpp"
#include "mainparm.cpp"
#include "frcmod.cpp"
#include "mapping.cpp"

#include <fstream>
#include <string>
#include <vector>
#include <map>
typedef std::vector<MolecularModeling::Atom*> AtomVector;

//argv:1, mol2 file; 2, ligand pdb file; 3, ligand pdb2glycam log file; 4,AMBER gaff dat file. 5, Input frcmod from antechamber. 6, output frcmod file. 7. Output off file. 
int main(int argc, char* argv[])
{
    std::vector<mol2_atom> mol2_atoms = ParseMol2File(argv[1]);
    std::vector<pdb2glycam_entry> pdb2glycam_info = ParsePdb2GlycamLogFile(argv[3]);

    MolecularModeling::Assembly ligand_pdb(std::string(argv[2]), gmml::InputFileType::PDB);
    ligand_pdb.SetName("corona");
    VinaBondByDistanceForPDB(ligand_pdb, 0);
    AtomVector pdb_atoms = ligand_pdb.GetAllAtomsOfAssembly();

    if (mol2_atoms.size() != pdb2glycam_info.size() || pdb2glycam_info.size() != pdb_atoms.size() || mol2_atoms.size() != pdb2glycam_info.size()){
        std::cout << "Error. Number of atoms in mol2, pdb2glycam log, and pdb file are not equal." << std::endl;
	    std::exit(1);
    }

    OverWriteGAFFParametersWithGLYCAM(pdb_atoms, mol2_atoms, pdb2glycam_info);
    std::vector<std::pair<MolecularModeling::Atom*, MolecularModeling::Atom*>> glycam_gaff_bonds = GetGlycamGaffBonds(pdb_atoms, pdb2glycam_info);

    ChargeAdjustment(mol2_atoms, pdb_atoms, glycam_gaff_bonds);

    std::vector<angle> glycam_gaff_angles;
    std::vector<torsion> glycam_gaff_torsions; 
    ProfileGlycamGaffAnglesAndTorsions(glycam_gaff_bonds, glycam_gaff_angles, glycam_gaff_torsions);
    std::vector<torsion> glycam_gaff_improper_torsions = DetectImproperTorsions(pdb_atoms, pdb2glycam_info);

    for (unsigned int i = 0; i < glycam_gaff_improper_torsions.size(); i++){
        std::cout << "Improper: " << glycam_gaff_improper_torsions[i].atom1_->GetName() << "-" << glycam_gaff_improper_torsions[i].atom2_->GetName() << "-" << glycam_gaff_improper_torsions[i].atom3_->GetName() << "-" << glycam_gaff_improper_torsions[i].atom4_->GetName() << std::endl;
    }

    MainParm amber_gaff_dat(argv[4]);
    Frcmod   antechamber_frcmod(argv[5]);

    std::vector<MassLine*> interface_masses;//Leave empty, should be irrelevant
    std::vector<BondLine*> interface_bond_lines;
    std::vector<AngleLine*> interface_angle_lines;
    std::vector<TorsionBlock*> interface_torsions;
    std::vector<TorsionBlock*> interface_improper_torsions;
    SixTwelvePotential* interface_nonbon = NULL;

    MapInterface2GaffParm(glycam_gaff_bonds, glycam_gaff_angles, glycam_gaff_torsions, glycam_gaff_improper_torsions, interface_bond_lines, interface_angle_lines, interface_torsions, 
		          interface_improper_torsions, interface_nonbon, amber_gaff_dat, antechamber_frcmod, pdb_atoms, mol2_atoms);

    Frcmod glycam_gaff_frcmod("Interfacing glycam and gaff", interface_masses, interface_bond_lines, interface_angle_lines, interface_torsions, interface_improper_torsions, interface_nonbon);
    std::string output_frcmod_file = std::string(argv[6]);
    glycam_gaff_frcmod.Write(output_frcmod_file);

    std::string output_off_file = std::string(argv[7]);
    ligand_pdb.CreateOffFileFromAssembly(output_off_file, 0);
    return 0;
}
