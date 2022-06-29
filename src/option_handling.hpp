#ifndef OPTION_HANDLING_HPP
#define OPTION_HANDLING_HPP

#include <unistd.h> //getopt
#include <string>
#include "open_valence_derivative_moiety.hpp"
#include <sstream>

bool ValidateOptions(int argc, char* argv[], bool& use_input_file){
    std::string input_file_argument = "-f";
    std::vector<std::string> required_options = {"-c", "-a", "-i", "-o", "-t", "-l"};
    std::vector<std::string> parsed_arguments; 
    bool option_valid = true;

    for (unsigned int i = 0; i < argc; i++){
        parsed_arguments.emplace_back(std::string(argv[i]));
    }

    if (std::find(parsed_arguments.begin(), parsed_arguments.end(), input_file_argument) != parsed_arguments.end()){
        use_input_file = true;
        return option_valid;
    } 

    for (unsigned int i = 0; i < required_options.size(); i++){
        std::string& option = required_options[i];
        if (std::find(parsed_arguments.begin(), parsed_arguments.end(), option) == parsed_arguments.end()){
            std::cout << "Required option " << option << " missing." << std::endl;
            option_valid =  false;
        }
    }

    return option_valid;

}

bool IsInteger(std::string& str){
    std::stringstream stream(str);
    int i;
    stream >> i;
    return (stream.eof() && !stream.fail());
}

void ProcessInputFile(int argc, char* argv[], int& interval, int& max_num_threads, int& num_threads, std::string& cocrystal_pdb_path, std::string& output_pdb_path, std::string& logfile_path, std::vector<open_valence_option>& open_valence_options){
    std::vector<std::string> required_keywords = {"ComplexPdb", "OpenValence", "Interval", "NumThreads", "OutputPath", "LogFile"};
    std::vector<bool> keywords_present(required_keywords.size(), false);
    std::string input_file_argument = "-f";

    int stop_index = argc -1;
    for (unsigned int i = 0; i < stop_index; i++){
        std::string this_argument(argv[i]);
        std::cout << "This argument: " << this_argument << std::endl;
        if (this_argument == input_file_argument){
            std::string input_file_path(argv[i+1]);
            std::cout << "Input file path: " << input_file_path << std::endl;

            std::ifstream input_file(input_file_path);
            if (input_file.fail()){
                std::cout << "Failed to open input file: " << input_file_path << " for reading. Aborting." << std::endl;
                std::exit(1);
            }

            std::string line;
            while (std::getline(input_file, line)){
                std::cout << "This line: " << line << std::endl;
                if (line.find(":") == std::string::npos){
                    std::cout << "Bad line: " << line << std::endl;
                    std::cout << "Each line must be key and value, separated by colon(:)" << std::endl;
                    std::exit(1);
                }

                std::vector<std::string> colon_split_token = gmml::Split(line, ":");
                std::string& keyword = colon_split_token[0];
                std::string& value = colon_split_token[1];

                if (std::find(required_keywords.begin(), required_keywords.end(), keyword) == required_keywords.end()){
                    std::cout << "Unrecognized keyword " << keyword << " in input file. Ignored." << std::endl;
                }
                else{
                    std::vector<std::string>::iterator it = std::find(required_keywords.begin(), required_keywords.end(), keyword);
                    int index = std::distance(required_keywords.begin(), it);
                    keywords_present[index] = true;

                    if (keyword == "ComplexPdb"){
                        cocrystal_pdb_path = value;
                    }
                    else if (keyword == "OpenValence"){
                        std::vector<std::string> dash_split_token = gmml::Split(value, "-");
                        open_valence_options.emplace_back(open_valence_option(dash_split_token));
                    }
                    else if (keyword == "Interval"){
                        if (IsInteger(value)){
                            interval = std::stoi(value);
                        }
                        else{
                            std::cout << "Not a valid integer: " << value << std::endl;
                            std::exit(1);
                        }

                    }
                    else if (keyword == "NumThreads"){
                        if (IsInteger(value)){
                            num_threads = std::stoi(value);
                        }
                        else{
                            std::cout << "Not a valid integer: " << value << std::endl;
                            std::exit(1);
                        }
                    }
                    else if (keyword == "OutputPath"){
                        output_pdb_path = value; 
                    }
                    else if (keyword == "LogFile"){
                        logfile_path = value;
                    }
                }
            }
            break;
        }     
    }

    return;
}

void ProcessOptions(int argc, char* argv[], int& interval, int& max_num_threads, int& num_threads, std::string& cocrystal_pdb_path, std::string& output_pdb_path, std::string& logfile_path, std::vector<open_valence_option>& open_valence_options, std::string requested_combinations){

    int option = 0;
    while ((option = getopt(argc, argv, "c:a:i:m:o:t:e:l:")) != -1){
        switch (option) {
            case 99: { //-c complex pdb path
                //std::cout << "Process -c" << std::endl;
                cocrystal_pdb_path = std::string(optarg);
                break;
            }
            case 97: { //-a open valence options.
                //Format: resnum_atom1_openvalence_atomToReplace:linkage_value(if present)-moiety_path-moiety_name_pattern-resnum'_atom1'_atom2'_atom3'_atom4':fixed value.
	            //Resnum: residue number this atom belongs to. Atom1-2, first 2 atoms of the linkage torsion. Atom3, the open valence atom. Atom4, the bonded neighbor of atom 3 (and its downstrem atoms) to be
	            //        replaced by an R group. True/false, indicate whether this linkage torsion is rotatable. Resnum' and atom', other torsion the user species for rotation upstream of the linkage torsion,
	            //        and if : is present, it indicates that this torsion is fixed to a value rather than freely rotatable.
	            //example: 3_O5_C1_O1_H1O-2_C5_C6_C7_C8:120
	            //example: 3_O5_C1_O1_H1O:150-2_C5_C6_C7_C8:120
	    
                //std::cout << "Process -a" << std::endl;
                std::string optarg_str = std::string(optarg);
                std::vector<std::string> dash_split_token = gmml::Split(optarg_str, "-"); 
                open_valence_options.emplace_back(open_valence_option(dash_split_token));
                break;
            }
            case 105: { //-i
                //std::cout << "Process -i" << std::endl;
                interval = std::stoi(std::string(optarg));
                break;
            }
            case 111:{ //-o
                //std::cout << "Process -o" << std::endl;
                output_pdb_path = std::string(optarg);
                break;
            }
            case 116:{ //-t
                //std::cout << "Process -t" << std::endl;
                num_threads = std::stoi(std::string(optarg));
                break;
            }
            case 108:{//-l for log file
                //std::cout << "Process -l" << std::endl;
                logfile_path = std::string(optarg);
                break;
            }
            /*case 114:{//-r for requested combintations to write, if >=2 open positions.
                requested_combinations = std::string(optarg);        
                break;
            }*/
        }
    }
}

std::vector<OpenValence*> GetOpenValenceSideGroupSp3Oxygens (CoComplex* cocomplex, std::vector<open_valence_option>& open_valence_options, int num_threads){
    std::vector<OpenValence*> open_valences;
    std::vector<MolecularModeling::Residue*> ligand_residues = cocomplex->GetLigandResidues();

    for(std::vector<open_valence_option>::iterator vec_it = open_valence_options.begin(); vec_it != open_valence_options.end(); vec_it++){
        bool this_open_option_found = false;

        for (unsigned int k = 0; k < ligand_residues.size(); k++){
            MolecularModeling::Residue* this_ligand_residue = ligand_residues[k];
            std::string this_ligand_residue_id = this_ligand_residue->GetId();

            std::vector<std::string> underscore_split_token = gmml::Split(this_ligand_residue_id, "_");
            int this_residue_number = 999999;
            if (underscore_split_token.size() == 6){ //Residue Id
                this_residue_number = std::stoi(underscore_split_token[2]);
            }
            else if (underscore_split_token.size() == 8){ //Atom Id
                this_residue_number = std::stoi(underscore_split_token[4]);
            }

            if (this_residue_number != 999999 && vec_it->residue_number == this_residue_number){
                AtomVector this_ligand_residue_atoms = this_ligand_residue->GetAtoms();
                for (unsigned int l = 0; l < this_ligand_residue_atoms.size(); l++){
                    MolecularModeling::Atom* atom = this_ligand_residue_atoms[l];
                    if (vec_it->atom_name == atom->GetName()){ //If same residue and same atom
                        OpenValence* new_open_valence = new OpenValence(cocomplex, atom, num_threads, vec_it->linkage_torsion_atom1, vec_it->derivatization_type, vec_it->linkage_preset, vec_it->linkage_preset_torsion_value_str, vec_it->explicit_torsion_str, vec_it->explicit_torsion_str_preset, vec_it->moiety_path, vec_it->moiety_name_pattern);
                        open_valences.push_back(new_open_valence);
                        this_open_option_found = true;
                        break;
                    }
                }
            }

        }

        if (!this_open_option_found){
            std::cout << "Cannot build open valence object for " << vec_it->residue_number << " and " << vec_it->atom_name << std::endl;
            std::exit(1);
        }
    }

    return open_valences;
}

#endif
