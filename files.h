#include <fstream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include "global.h"
#include "functions.h"
using namespace Eigen;

void get_thread_contacts(int step, int threadnum, bool norm) {
    if (norm) {
        std::ostringstream fn;
        fn << forward_output_dir << "contacts/"<< step << "/contacts_" << threadnum << ".txt";
        std::ofstream final_cont;
        final_cont.open(fn.str().c_str(), std::ios_base::binary); //write contact frequencies
        std::vector<std::vector<double>> norm_map(bin_num, std::vector<double>(bin_num, 0));
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j ++){
                norm_map[i][j] += total_contacts[threadnum][i][j];
            }
        }
        normalize(norm_map);
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                double contact;
                if (std::abs(i-j) <= 1 || (i==0 && j == bin_num - 1) || (j==0 && i == bin_num - 1)) { contact = 0; }
                else {
                    contact =norm_map[i][j];
                }
                final_cont << contact << ' ';
            }
            final_cont << '\n';
        }
        final_cont.close();
    }
    else {
        std::ostringstream fn;
        fn << forward_output_dir << "contacts/"<< step << "/contacts_" << threadnum << ".txt";
        std::ofstream final_cont;
        final_cont.open(fn.str().c_str(), std::ios_base::binary); //write contact frequencies
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                total_contacts[threadnum][j][i] = total_contacts[threadnum][i][j];
            }
        }
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                double contact;
                if (i==j) { contact = 0; }
                else {
                    contact = total_contacts[threadnum][i][j]/mc_moves;
                }
                final_cont << contact << ' ';
            }
            final_cont << '\n';
        }
        final_cont.close();
    }
}


void check_input_configuration_compatibility(std::string& old_simulation_folder){
    std::ifstream input_file;
    std::string word = " ";
    input_file.open(dir + old_simulation_folder + "sim_params.txt");
    if(input_file.fail()){
        throw std::invalid_argument( "Input configuration file doesnt exist (or sim_params.txt missing)." );}
    while(word != "Thread"){
        input_file >> word;
    }
    input_file >> word;
    input_file >> word; // word = number of threads
    if(number_of_threads > std::stoi(word)){
        throw std::invalid_argument( "Incompatible configuration input (not enough configurations)." );
    }
}


void read_configuration(int thread_num, int step) {
    std::ostringstream filename; 
    filename << dir << configuration_data_folder << "Configurations/thread_"<<thread_num<<"/config_" << step <<".txt";
    std::ifstream monomers;
    monomers.open(filename.str());  //read in monomer positions
    if(monomers.fail()){
        std::cerr<<"configuration's "+filename.str()+" input file doesnt exist"<<std::endl;
        throw std::invalid_argument( "configuration's "+filename.str()+" input file doesnt exist" );
    }

    for (int i = 0; i < pol_length; i++) {
        Vector3i location = {0,0,0};
        for (int j = 0; j < 3; j++) {
            monomers >> location[j];
        }
        polymer[thread_num].push_back(location);
    }
    monomers.close();
}

void read_input_data(){
    std::ifstream input_file;
    std::string word;
    int site;
    float exp_data;
    std::ostringstream filename;
    filename << dir << "Input/" << bacteria_name << "/" << bacteria_name << ".txt";
    input_file.open(filename.str());
    if(input_file.fail()){throw std::invalid_argument( "Missing input data for "+ bacteria_name +"." );}
    while (!input_file.eof()) {
        input_file >> word;
        if(word == "length"){
            input_file >> length;
        }
    }
    input_file.close();
}


void get_present_configs(int thread_num, int sample_num) {
    std::ostringstream fn1;
    fn1 << forward_output_dir << output_folder << "/Configurations/thread_" << thread_num << "/config_" << sample_num << ".txt";
    std::ofstream config_file;
    config_file.open(fn1.str().c_str(), std::ios_base::binary);

    for (int i = 0; i < pol_length; i++) {
        for (int j = 0; j < 3; j++) {
            config_file << polymer[thread_num][i][j] << " ";
        }
        config_file << '\n';
    }
    config_file.close();
 
    //also save the positions of the crosslinks
    std::ostringstream fn3;
    fn3 << forward_output_dir << output_folder << "/Configurations/thread_" << thread_num << "/crosslinks_" << sample_num << ".txt";
    std::ofstream crosslink_file;
    crosslink_file.open(fn3.str().c_str(), std::ios_base::binary);

    for (int i = 0; i < bin_num; i++) {
        if (is_crosslinked[thread_num][i]){
            for (auto other : crosslinked_with[thread_num][i]) {
                crosslink_file<< other<< ", ";
            }
        }
        crosslink_file << '\n';
    }
    crosslink_file.close();
}

void open_safely(std::string filename, std::ofstream &streamname){
    streamname.open(filename);
    if(!streamname){
        throw std::invalid_argument("Couldn't open output file "+filename+". Missing directory?");
    }
}

//create file for saving contacts
void save_contacts(int thread_num, std::unordered_map<std::vector<int>, std::vector<int>, vec_hash> locations_dict){
    std::ostringstream filename;
    filename<<forward_output_dir<<output_folder<<"/Configurations/thread_"<<thread_num<<"/locations.txt";
    std::ofstream final_contact_list;
    open_safely(filename.str(), final_contact_list);
    for (auto & entry : locations_dict){ //key: location, item: all the monomers there
        if (entry.second.size()>1){ //only store monomers in contact
            for (auto & index : entry.second){
                final_contact_list<<index<<" "; //add every item in contact at this point
            }
            final_contact_list<<std::endl; //start new line before next entry
        }
    }
    final_contact_list.close();
}

//append to file saving contacts
void add_contacts(int thread_num, std::unordered_map<std::vector<int>, std::vector<int>, vec_hash> locations_dict){
    std::ofstream final_contact_list;
    std::ostringstream filename;
    filename<<forward_output_dir<<output_folder<<"/Configurations/thread_"<<thread_num<<"/locations.txt";
    final_contact_list.open(filename.str(),std::ios_base::app);

    for (auto & entry : locations_dict){ //key: location, item: all the monomers there

        if (entry.second.size()>1){ //only store monomers in contact
            for (auto & index : entry.second){
                final_contact_list<<index<<" "; //add every item in contact at this point
            }
            final_contact_list<<std::endl; //start new line before next entry       
        }
    }       

    final_contact_list<<std::endl; //empty line to separate the next config
    final_contact_list.close();
}

void get_affinities(){
    std::ostringstream fn;
    fn << forward_output_dir << output_folder << "/affinities.txt";
    std::ofstream affinities;
    affinities.open(fn.str().c_str(), std::ios_base::binary);
    for (int i = 0; i < bin_num; i++) {
        affinities<<binding_affinities[i] << '\n';
    }
    affinities.close();
}

void get_final_contacts(const int& step) {
    //all contacts from all threads
    std::ostringstream fn;
    fn << forward_output_dir << output_folder << "/" << "Contacts/contacts_" << step << ".txt";
    std::ofstream final_cont;
    final_cont.open(fn.str().c_str(), std::ios_base::binary); //write contact frequencies
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            double contact;
            if (std::abs(i-j) <= 1 || (i==0 && j == bin_num - 1) || (j==0 && i == bin_num - 1)) { contact = 0; }
            else {
                contact = final_contacts[i][j];
            }
            final_cont << contact << ' ';
        }
        final_cont << '\n';
    }
    final_cont.close();
}


void get_sim_params() {
    std::ostringstream fn;
    fn << forward_output_dir << output_folder << "/" << "sim_params.txt";
    std::ofstream params;
    params.open(fn.str().c_str(), std::ios_base::binary);
    params << "Thread number: " << number_of_threads << '\n';
    params << "Bin number: " << bin_num << '\n';
    params << "Polymer length: " << pol_length << '\n';
    params << "MC moves: " << mc_moves << '\n';
    params << "Bacteria: " << bacteria_name << '\n';
    params << "Periodic potential: " << periodic_en << '\n';
    params << "crosslink_mu: " << crosslink_mu << '\n';
    params << "crosslink energy amplitude: " << energy_amp<< '\n';
    
    params << '\n';

    if(initConfig) {
        params << "Input configurations folder: " << configuration_data_folder << '\n';
    }

    params << '\n';
    params << "Confinement/cell shape:" << '\n';
    params << "cell radius: " << radius << '\n';
    params << "cell length: "<< length << "\n";
    params << "offset: " << offset[0] << " "<< offset[1] << '\n';
}
