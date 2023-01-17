#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

template <typename T> void hash_combine(std::size_t& seed, const T& v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& pair) const noexcept {
        std::size_t seed = 0;
        hash_combine(seed, pair.first);
        hash_combine(seed, pair.second);
        return seed;
    }
};
struct vec_hash {
    std::size_t operator()(const std::vector<int> vec) const noexcept {
        std::size_t seed = 0;
        for (auto el : vec) {
            hash_combine(seed,el);
        }
        return seed;
    }
};

extern std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts;
extern std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations;

extern double crosslink_mu; 
extern bool periodic_en;
extern double energy_amp;
extern double p_cl;
extern double p_det;
extern std::vector<int> current_crosslinks;
extern std::vector<std::unordered_map<int, std::vector<int>>> crosslinked_with; //mon -> all crosslinks
extern std::vector<std::vector<bool>> is_crosslinked;
extern std::vector<double> binding_affinities;

extern const int number_of_threads;
extern std::vector<RandomGenerator> generators;
extern const int mc_moves;


extern std::string dir;
extern std::string output_folder;
extern std::string forward_output_dir;
extern std::string output_folder;

extern bool initConfig;
extern std::string configuration_data_folder;

extern const std::string bacteria_name;
extern const int bin_num;
extern const int pol_length;
extern const int reduction_factor;
extern std::vector<std::vector<Eigen::Vector3i>> polymer;

extern const int res;

extern bool boundary_cond;
extern double radius;
extern double length;
extern const int x_axis;
extern const int y_axis;
extern const int z_axis;
extern std::vector<double> offset;
extern double offset_z;

extern std::vector<std::vector< std::vector<double>>> total_contacts;
extern std::vector< std::vector<double>> final_contacts;


#endif
