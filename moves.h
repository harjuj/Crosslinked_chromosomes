#ifndef MV_H
#define MV_H

#include <algorithm>
#include "global.h"
#include "crosslinks.h"

//Note locations tracks physical monomers, ie pol index/reduction_factor. Same as indices of total_contacts and contacts
void update_contacts(int thread_num, int moved_site, Eigen::Vector3i prop_move, int m){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;
    //Throw away old contacts_lin[thread_num], at the same time update contact frequency map
    for (auto elem : locations[thread_num][{polymer[thread_num][moved_site][0], polymer[thread_num][moved_site][1], polymer[thread_num][moved_site][2]}]) {
        if (elem != red_site) {
            int i = std::min(elem, red_site);
            int j = std::max(elem, red_site);
            total_contacts[thread_num][i][j] += double(m) - contacts[thread_num][{i,j}];
            // GG: here one can save other total contacts
            contacts[thread_num].erase({i, j});
        }
    }
    //put in new contacts[thread_num]
    if (locations[thread_num].find({prop_move[0], prop_move[1], prop_move[2]}) != locations[thread_num].end()) {
        for (auto elem : locations[thread_num][{prop_move[0], prop_move[1], prop_move[2]}]) {
            contacts[thread_num][{std::min(elem, red_site),std::max(elem, red_site)}] = m;
        }
    }
}



void update_locations(int thread_num, int moved_site, Eigen::Vector3i prop_move){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;

    //update hash map locations[thread_num]
    std::vector<int> site_coords = { polymer[thread_num][moved_site][0], polymer[thread_num][moved_site][1], polymer[thread_num][moved_site][2] };
    std::vector<int> prop_move_vec = {prop_move[0], prop_move[1], prop_move[2]};
    // remove old locations
    if (locations[thread_num][site_coords].size() == 1) {
        locations[thread_num].erase(site_coords);
    }
    else {
        locations[thread_num][site_coords].erase(std::find(locations[thread_num][site_coords].begin(), locations[thread_num][site_coords].end(), red_site));
    }
    // save new locations
    if (locations[thread_num].find(prop_move_vec) != locations[thread_num].end()) {
        locations[thread_num][prop_move_vec].push_back(red_site);
    }
    else {
        locations[thread_num][prop_move_vec] = { red_site };
    }
}


bool check_boundary_rest(Eigen::Vector3i prop_move1, int thread_num) {
    if (boundary_cond == 1) {
        bool accept_1;
        // GG: check for prop_move1
        if (std::abs(prop_move1[2] - offset_z) <= length / 2) { //GG: is inside the cylinder, not in a cap (z-coord)
            accept_1 = (pow(prop_move1[0]+offset[0], 2) + pow(prop_move1[1]+offset[1], 2) <= pow(radius, 2));
        }
        else if (std::abs(prop_move1[2] - offset_z) <= length / 2 + radius) { // GG: is in one of the caps (z-coord)
            accept_1 = (pow(prop_move1[0]+offset[0], 2) + pow(prop_move1[1]+offset[1], 2) + pow(std::abs(prop_move1[2] - offset_z) - length / 2, 2) <= pow(radius, 2));
        }
        else {
            accept_1 = 0;
        }

        return accept_1;
    }
    else { return 1; }
}

void kink_move(int thread_num, int site, int &m) {
    if ((polymer[thread_num][(site + 2) % pol_length] != polymer[thread_num][site]) && (polymer[thread_num][(site + 2) % pol_length] != 2 * polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site])) {
        Eigen::Vector3i prop_move1;
        prop_move1 = polymer[thread_num][site] + polymer[thread_num][(site + 2) % pol_length] - polymer[thread_num][(site + 1) % pol_length];
        if (check_boundary_rest(prop_move1, thread_num) == 1) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}

void crankshaft_move(int thread_num, int site, int &m) {
    if ((polymer[thread_num][(site + 2) % pol_length] != polymer[thread_num][site]) && (polymer[thread_num][(site + 2) % pol_length] != 2 * polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site]) && (polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][(site + 2) % pol_length] == polymer[thread_num][site] - polymer[thread_num][(site + 1) % pol_length])) {
        int direction = generators[thread_num].unidir();
        Eigen::Vector3i prop_move1; Eigen::Vector3i prop_move2;
        if (direction == 1) { //180 degree flip
            prop_move1 = 2 * polymer[thread_num][site] - polymer[thread_num][(site + 1) % pol_length];
            prop_move2 = 2 * polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][(site + 2) % pol_length];
        }
        else { //90 degree flip
            prop_move1 = polymer[thread_num][site] + (direction - 1) * (polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site]).cross(polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][site]);
            prop_move2 = polymer[thread_num][(site + 3) % pol_length] + (direction - 1) * (polymer[thread_num][(site + 2) % pol_length] - polymer[thread_num][(site + 3) % pol_length]).cross(polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][site]);;
        }
        if (check_boundary_rest(prop_move1, thread_num) == 1 && check_boundary_rest(prop_move2, thread_num) == 1) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            update_contacts(thread_num, (site + 2)%pol_length, prop_move2, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);
            update_locations(thread_num, (site + 2)%pol_length, prop_move2);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
            polymer[thread_num][(site + 2) % pol_length] = prop_move2;
        }
        m++;
    }
}

void loop_move(int thread_num, int site, int &m) {
    if (polymer[thread_num][site] == polymer[thread_num][(site + 2) % pol_length]) {
        int direction = generators[thread_num].unidir_loop();
        Eigen::Vector3i rotated_vector(3); Eigen::Vector3i prop_move1;
        for (int i = 0; i < 3; i++) { // rotate loop in one of 5 possible new directions
            rotated_vector[(i + direction / 2) % 3] = (-2 * (direction % 2) + 1) * (polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site])[i];
        }
        prop_move1 = polymer[thread_num][site] + rotated_vector;
        if (check_boundary_rest(prop_move1, thread_num) == 1) {
            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}


void crosslink_moves(int thread_num) {
    //First try detaching
    int site = generators[thread_num].unisitering();

    try_detach(site, thread_num);

    //Then try attaching
    site = generators[thread_num].unisitering();
    try_attach(site, thread_num);
}

void move(int thread_num, int &m) {
    //Pick a non-crosslinked site on chromosome: 
    int site;
    do{
        site= generators[thread_num].unisitering();
    }while(site%reduction_factor==0 && is_crosslinked[thread_num][site/reduction_factor]);

    //Now attempt a polymer move on the non-crosslinked site
 
    int action = generators[thread_num].unimove();
    if (action == 0) {
        kink_move(thread_num, site, m);
    } else if (action == 1) {
        crankshaft_move(thread_num, site, m);
    } else {
        loop_move(thread_num, site, m);
    }
}

#endif
