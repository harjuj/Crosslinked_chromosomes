#ifndef INIT_H
#define INIT_H

#include "global.h"
#include "moves.h"
using namespace Eigen;

/*
All of these maps store only physical monomers, ie polymer index/reduction_factor
Makes maps 4*smaller, meaning searches are a lot faster.
Note that the polymer array uses full indices; keeps track of unphysical monomers
*/
std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations(number_of_threads);

//Given physical monomer index, return vector of crosslinked monomers
std::vector<std::unordered_map<int, std::vector<int>>> crosslinked_with(number_of_threads);

void generate_affinities(bool periodic_en){
    if (periodic_en){
        for (int i=0; i<bin_num; i++){
            binding_affinities[i]=energy_amp*cos(i*2*M_PI/bin_num);
        }
    }
    else{
        for (int i=0; i<bin_num; i++){
            binding_affinities[i]=energy_amp*2*(generators[0].disReal()-0.5);
        }
    }
}

void initialize(int thread_num, int step) {
    // Construct initial configurations if these are not taken from saved configurations//
    if(!initConfig) {
        int xi{0}, yi{0}, zi{0};
        Vector3i monomer;
        for (int i = 0; i < pol_length; i++) {
            if (i % 2 == 0) {
                zi -= 1;
            } else {
                zi += 1;
            }
            monomer = {xi, yi, zi};
            polymer[thread_num].push_back(monomer);
        }
    }
    else{
        // Construct initial configurations from saved configurations//
        read_configuration(thread_num, step);
    }


    // set contacts and locations. ONLY TRACK PHYSICAL CONTACTS. //
    for (int i = 0; i < pol_length; i+=reduction_factor) { //Ring ring contacts
        int red_i = i/reduction_factor;
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(red_i);
    }
}


void burn_in(int thread_num, int n_steps) {
    for (int m = 0; m < n_steps ; m++) {
        move(thread_num, m);
    }

    //Empty the contacts found during burn-in
    std::vector<double> zeroVec(bin_num, 0);
    std::fill(total_contacts[thread_num].begin(), total_contacts[thread_num].end(), zeroVec);
    
    locations[thread_num].clear();
    contacts[thread_num].clear();

    //Fill in starting contacts and locations
    for (int i = 0; i < pol_length; i+=reduction_factor) {
        int red_i = i/reduction_factor;
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(red_i);
    }
}

#endif
