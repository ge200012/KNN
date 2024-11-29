#ifndef INDEXFILTER_H
#define INDEXFILTER_H

#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "distance.h"

//extern std::map<std::pair<int, int>, std::set<int>> invert_index;

//void add_index(const path& trajectory, int idx, double grid_size);
//std::set<int> run_get_index(const std::vector<path>& paths, const std::vector<value>& result_grid, double grid_size);
std::pair<std::set<int>, std::map<std::pair<int, int>, std::set<int>>> run_get_index(const std::vector<path>& paths, const std::vector<value>& result_grid, double grid_size); 
    
// std::pair<std::set<int>, double> get_candi_osf(
//     const path& x,
//     const std::map<std::pair<int, int>, std::set<int>> invert_index,
//     double rate,
//     double grid_size,
//     const std::vector<path>& paths,
//     const std::vector<value>& result_grid,
//     const std::string& dataset,
//     const std::string& metric
// );
// std::pair<std::set<int>, double> osf_filter(
//     const path& current_trajectory,
//     double radius,
//     double grid_size,
//     const std::vector<path>& paths,
//     const std::vector<value>& result_grid,
//     const std::string& dataset,
//     const std::string& metric
// );
std::pair<std::set<int>, double> get_candi_osf(
    const path& current_trajectory,
    const std::map<std::pair<int, int>, std::set<int>> invert_index,
    double rate,
    double grid_size,
    const std::vector<path>& paths,
    const std::vector<value>& result_grid,
    const std::string& dataset,
    const std::string& metric);

std::pair<std::set<int>, double> osf_filter(
    const path& current_trajectory,
    const std::map<std::pair<int, int>, std::set<int>> invert_index,
    double radius,
    double grid_size,
    const std::vector<path>& paths,
    const std::vector<value>& result_grid,
    const std::string& dataset,
    const std::string& metric);

#endif // INDEXFILTER_H
