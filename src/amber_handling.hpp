#ifndef AMBER_HANDLING_HPP
#define AMBER_HANDLING_HPP

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
//#include "../../gmml/src/MolecularMetadata/guesses.cc"


void RenameDisulfideCYS2CYX(std::vector<MolecularModeling::Residue*>& assembly_residues){
   for (unsigned int i = 0; i < assembly_residues.size(); i++){
       MolecularModeling::Residue* residue = assembly_residues[i];
       std::string resname = residue->GetName();

       if (resname == "CYS"){
           AtomVector res_atoms = residue->GetAtoms();    
	   for (unsigned int j = 0; j < res_atoms.size(); j++){
	       MolecularModeling::Atom* this_res_atom = res_atoms[j];

	       if (this_res_atom->GetName() == "SG"){
	           AtomVector sg_atom_neighbors = this_res_atom->GetNode()->GetNodeNeighbors();
		   for (unsigned int k = 0; k < sg_atom_neighbors.size(); k++){
		       MolecularModeling::Atom* this_sg_neighbor = sg_atom_neighbors[k];
		       MolecularModeling::Residue* neighbor_residue = this_sg_neighbor->GetResidue();

		       //If this SG neighbor is bonded to another SG atom in another residue, then both CYS gets renames to CYX
		       if (neighbor_residue != residue && neighbor_residue->GetName() == "CYS" && this_sg_neighbor->GetName() == "SG"){
		           residue->SetName("CYX");
			   neighbor_residue->SetName("CYX");
		       }
		   }
	       }
	   }
       }
   }

}

#endif
