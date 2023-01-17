#include <thread>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include "RandomGenerator.h"
#include "global.h"
#include "functions.h"
#include "files.h"
#include "initialization.h"
#include "boost/lexical_cast.hpp"

////// Directories for data //////////
std::string forward_output_dir = std::filesystem::current_path().string()+"/out/";//where results stored
std::string dir = std::filesystem::current_path().string()+"/initial_configs/"; //where you want to load previous configurations from, if initConfig=true
std::string output_folder;//time-stamp as directory name in ./out, defined on run

////// Simulation properties //////
const int number_of_threads = 30; //each thread runs own simulation; faster convergence with more 
const std::string bacteria_name = "crescentus_separations_3";
const int bin_num = 405; //Crescentus
const int reduction_factor = 4; //monomers per bin

const int pol_length = reduction_factor*bin_num;
const int mc_moves = 2e9;
const int burn_in_steps = 1e7;
const int config_sample_freq = 1000000;


////// cell size /////
bool boundary_cond = true; //set to false for no confinement
double radius{ 3.2 }; // Crescentus
double length{18}; // JM: cell lengths. Ecoli: 8-21
std::vector<double> offset {0.5,0.5}; //JM: offsets in x and y directions, determines the center of the cylinder
double offset_z{(int(length) % 2)/2.};

////////// initial configurations //////////
bool initConfig;
std::string configuration_data_folder;
int init_config_number = 200; //iteration number of the initial configuration used


/////// variables to be used later ///////

std::vector<RandomGenerator> generators;
std::vector<std::vector<std::vector<double>>> total_contacts(number_of_threads, std::vector< std::vector<double>>(bin_num, std::vector<double>(bin_num, 0))); // This stores the contact frequencies during the simulation
std::vector<std::vector<Eigen::Vector3i>> polymer(number_of_threads);
std::vector<std::vector<double>> final_contacts(bin_num, std::vector<double>(bin_num, 0));

//Crosslink related variables
double crosslink_mu;//command line
bool periodic_en;
double energy_amp;
double p_cl;
double p_det;
std::vector<int> current_crosslinks(bin_num, 0);
std::vector<std::vector<bool>> is_crosslinked(number_of_threads, std::vector<bool>(bin_num, false));
std::vector<double> binding_affinities(bin_num, 0);

/////// functions to be used later ///////
void run(int thread_num, int move_num);//monte carlo step
void clean_up();//clear arrays
void create_output_dirs();


/////// script ///////
int main(int argc, char **argv) {
    if (argc==4){
        initConfig=false;
        periodic_en=std::atoi(argv[1]);
        crosslink_mu=std::atof(argv[2]);
        energy_amp=atof(argv[3]);
    }
    else if(argc==5){
        initConfig=true;
        periodic_en=std::atoi(argv[1]);
        crosslink_mu=std::atof(argv[2]);
        energy_amp=atof(argv[3]);
        configuration_data_folder=argv[4];
        //folder should end in /
        if (configuration_data_folder.back() != '/'){
            configuration_data_folder+='/';
        }
    }
    else{
        std::cerr<<"Wrong number of arguments!\nPlease call ./Forward_crosslinks [periodic 1/random 0] [crosslink mu] [energy_amplitude] ([input_config_folder])"<<std::endl;
        exit(EXIT_FAILURE);
    }
    //set the probabilities
    p_cl=std::min(0.5, 0.5*exp(-crosslink_mu));
    p_det=std::min(0.5, 0.5*exp(crosslink_mu));
    //Print details as a check
    std::cout << "bin number: " << bin_num << '\n' << "polymer length: " << pol_length << '\n' <<"crosslink mu: "<<crosslink_mu<<"\nperiodic energies: "<<periodic_en<<"\nenergy amplitude: "<<energy_amp<<"\nmc moves: " << mc_moves <<'\n' << "boundaries: " << boundary_cond << '\n' << '\n' << "saved in " << forward_output_dir <<'\n';

    std::cout<<"creating folders"<<std::endl;
    create_output_dirs(); //empty directories in .out/

    // create and seed random generators. Change seed to get different results! //
    for (int i=0; i<number_of_threads; i++){
        generators.push_back( RandomGenerator(i+int(time(0)), pol_length) );
    }

    generate_affinities(periodic_en);
    get_affinities(); //save to file
    //Save the input parameters of the simulation in "sim_params.txt"//
    get_sim_params();

    std::cout<<"initialising"<<std::endl;
    //Initialize configurations//
    std::vector<std::thread> iniThreads(number_of_threads);
    for (auto l = 0; l < number_of_threads; l++) {
        iniThreads[l] = std::thread(initialize, l,init_config_number);
    }
    for (auto&& l : iniThreads) {
        l.join();
    }


    std::cout<<"start burning in"<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    //Burn in  configurations//
    std::vector<std::thread> burnThreads(number_of_threads);
    for (auto l = 0; l < number_of_threads; l++) {
        burnThreads[l] = std::thread(burn_in, l,burn_in_steps);
    }
    for (auto&& l : burnThreads) {
        l.join();
    }
    
    //create files for saving contacts
    for (int i=0; i<number_of_threads;i++){
        save_contacts(i, locations[i]); //save contacts by saving the locations container
    }

    //Run the Monte Carlo algorithm//
    int step = 0;
    std::cout<<"start simulations"<<std::endl;

    std::vector<std::thread> threads(number_of_threads);
    for (auto l = 0; l < number_of_threads; l++) {
        threads[l] = std::thread(run, l, std::round(mc_moves * sqrt(step+1)));//saves configurations during run
    }
    for (auto&& l : threads) {
        l.join();
    }

    //create a final Hi-C map
    for (int l = 0; l < number_of_threads; l++) {
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                final_contacts[i][j] += total_contacts[l][i][j];
            }
        }
    }
    //symmetrise:
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            final_contacts[j][i] = final_contacts[i][j];
        }
    }
    std::cout<<"Done! Saving final contacts."<<std::endl;

    //save final contact matrix//
    get_final_contacts(step);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}

//Function definitions

void run(int thread_num, int move_num) {
    int m = 0;
    int old_m;
    while (m < move_num) { //perform polymer simulation
        old_m = m;
        move(thread_num, m); //GG: "m" is passed by reference now
        crosslink_moves(thread_num);
        if(m > old_m) { //GG: not to include the same data multiple times when no move was made!
            //sample configurations and crosslink positions
            if (m%config_sample_freq==0){
                get_present_configs(thread_num, m/config_sample_freq);
                add_contacts(thread_num, locations[thread_num]); //save contacts by saving the locations container
            }
        }
    }
    
    for (const auto& elem : contacts[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
}

void create_output_dirs(){
    time_t time_now = time(0);   // get time now
    struct tm * now = localtime( & time_now );
    char buffer [20];
    strftime(buffer, 20, "%Y-%m-%d", now);
    std::string buffer_str = buffer; //converts to string
    if (periodic_en){
        output_folder = buffer_str+"_mu_"+boost::lexical_cast<std::string>(crosslink_mu)+"_"+boost::lexical_cast<std::string>(energy_amp)+"_cos";
    }else{
        output_folder = buffer_str+"_mu_"+boost::lexical_cast<std::string>(crosslink_mu)+"_"+boost::lexical_cast<std::string>(energy_amp)+"_rand";
    }

    std::string command = "mkdir " + forward_output_dir + output_folder;
    std::system(command.c_str());

    command = "mkdir " + forward_output_dir + output_folder + "/Configurations";
    std::system(command.c_str());
    for (int l=0; l<number_of_threads; l++){
        command = "mkdir " + forward_output_dir + output_folder + "/Configurations/thread_" + std::to_string(l);
        std::system(command.c_str());
    }
    command = "mkdir " + forward_output_dir + output_folder + "/Contacts";
    std::system(command.c_str());
}

void clean_up() {//empty vectors and arrays
    std::vector<double> zeroVec(bin_num, 0);
    std::fill(final_contacts.begin(), final_contacts.end(), zeroVec);
    //remove all crosslinks
    std::fill(is_crosslinked.begin(),is_crosslinked.end(),std::vector<bool>(bin_num,false));
    std::fill(current_crosslinks.begin(),current_crosslinks.end(), 0);

    for (int l = 0; l < number_of_threads; l++) {
        std::fill(total_contacts[l].begin(), total_contacts[l].end(), zeroVec);
        contacts[l].clear();
        crosslinked_with[l].clear();
        for (int i = 0; i < pol_length; i+=reduction_factor) {
            int red_i=i/reduction_factor;
            if (locations[l].find({polymer[l][i][0], polymer[l][i][1], polymer[l][i][2]}) != locations[l].end()) {
                for (auto elem : locations[l][{polymer[l][i][0], polymer[l][i][1],polymer[l][i][2]}]) {
                    contacts[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                }
            }
        }
    }
}
