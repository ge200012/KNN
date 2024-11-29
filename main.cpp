#include <iostream>
#include <vector>
#include "util/matrix.h"
#include "algorithms/knn_graph.h"
#include "algorithms/gnns_index.h"
#include "util/dist.h"
#include "io.h"
#include "define.h"
#include "util/params.h"
#include "evaluation/evaluation.h"

#include "mostSimilar.h"
//自己添加的trajectory内容
#include "utils.h"
#include <csignal>

#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cstdlib>
#include <string>
#include <Python.h>
#include <tuple>
#include <chrono>

#include <nlohmann/json.hpp>
#include <fstream>
#include <utility>
namespace py=pybind11;



using namespace std;
using json = nlohmann::json;
//using namespace gnns;
//自己添加的trajectory内容


//xian数据集保存位置
// const string graph_index_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/saved_graph_index30000_test";
// const string graph_dist_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/saved_graph_dist30000_test";
// //分类
// const string graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/saved_graph_start_end_position30000_test";

// const string top_level_graph_points_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/top_level_graph_points30000_test";
// const string top_level_graph_indices_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/top_level_graph_indices30000_test";
// const string top_level_graph_distances_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/top_level_graph_distances30000_test";
// //分类
// const string top_level_graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/test30000/top_level_graph_start_end_position30000_test";

// const string indices_queryies_groundtruth_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_xian/groundtruth_results/dtw/graph_CMA/groundtruth_results30000.txt";

// const string top_k_results_to_RR_saved_path="/home/KNN/orders_data2000query_notN/top_k_results/xian/DTW/graph_CMA/top_k_results_to_RR8start_test.csv";

//chendu数据集保存位置
const string graph_index_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/saved_graph_index30000_test";
const string graph_dist_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/saved_graph_dist30000_test";
//分类
const string graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/saved_graph_start_end_position30000_test";

const string top_level_graph_points_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/top_level_graph_points30000_test";
const string top_level_graph_indices_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/top_level_graph_indices30000_test";
const string top_level_graph_distances_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/top_level_graph_distances30000_test";
//分类
const string top_level_graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/test30000/top_level_graph_start_end_position30000_test";

const string indices_queryies_groundtruth_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_chengdu/groundtruth_results/dtw/graph_CMA/groundtruth_results30000.txt";

const string top_k_results_to_RR_saved_path="/home/KNN/orders_data2000query_notN/top_k_results/chengdu/DTW/graph_CMA/top_k_results_to_RR8start_test.csv";

//proto数据集保存位置
// const string graph_index_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/saved_graph_index3000_test";
// const string graph_dist_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/saved_graph_dist3000_test";
// //分类
// const string graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/saved_graph_start_end_position3000_test";

// const string top_level_graph_points_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/top_level_graph_points3000_test";
// const string top_level_graph_indices_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/top_level_graph_indices3000_test";
// const string top_level_graph_distances_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/top_level_graph_distances3000_test";
// //分类
// const string top_level_graph_start_end_position_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/test3000/top_level_graph_start_end_position3000_test";

// const string indices_queryies_groundtruth_saved_path="/home/KNN/orders_data2000query_notN/trajectory_new_porto/groundtruth_results/CMA_dtw/groundtruth_results3000.txt";

// const string top_k_results_to_RR_saved_path="/home/KNN/orders_data2000query_notN/top_k_results/porto/DTW/graph_CMA/top_k_results_to_RR1.csv";

//创建graph时计算data-data的minSubTrajectory2中的参数algorithm,但没有一点用，空作为一个参数
const string algorithm="efficientWED";

void signalHandler(int signum)
{
    std::cerr<<"Caught signal"<<signum<<std::endl;
}

struct Config {
    std::string data_path;
    int query_min_len;
    int query_max_len;
    int data_min_len;
    int data_max_len;
    std::string metric;
    std::string task;
};
void visualize_trajectories(const std::vector<std::vector<std::pair<double,double>>>& trajectories,size_t num_to_visualize)
{
    py::scoped_interpreter guard{};

    py::module plt=py::module::import("matplotlib.pyplot");

    size_t limit=std::min(trajectories.size(),num_to_visualize);

    for(size_t i=0;i<limit;++i)
    {
        const auto& trajectory=trajectories[i];
        std::vector<double> latitudes;
        std::vector<double> longitudes;

        std::cout<<"第"<<i<<"条轨迹"<<std::endl;
        for(const auto& point:trajectory)
        {
            latitudes.push_back(point.first);
            longitudes.push_back(point.second);
            std::cout<<point.first<<" "<<point.second<<std::endl;
        }

        plt.attr("plot")(longitudes,latitudes,"b-",py::arg("linewidth")=2);

    }

    plt.attr("xlabel")("Longitude");
    plt.attr("ylabel")("Latitude");
    plt.attr("title")("Trajectory Visualization (First " + std::to_string(limit)+" Trajectories)");

    plt.attr("savefig")("porto100_trajectories.png");
    plt.attr("clf")();
}

void write_results_to_file(const string& filename, const vector<vector<int>>& indices, const vector<vector<float>>& dists)
{
    std::ofstream file(filename);

    if(!file.is_open())
    {
        std::cerr<<"无法打开文件: "<<filename<<std::endl;
        return;
    }

    file<<"Query Index,Neighbor Index,Distance\n";

    for(size_t q=0;q<indices.size();++q)
    {
        for(size_t k=0;k<indices[q].size();++k)
        {
            file<<q<<","<<indices[q][k]<<","<<dists[q][k]<<"\n";
        }
    }

    file.close();

}
void saveTables(const std::vector<std::vector<std::pair<double,double>>>& t_table,
                const std::vector<std::vector<std::pair<double,double>>>& q_table)
{
    json t_json(t_table);
    json q_json(q_table);

    std::ofstream t_file("3w_data/origin/t_table_xian.json");
    std::ofstream q_file("3w_data/origin/q_table_xian.json");

    t_file << t_json.dump(4); // pretty-print with 4 spaces of indentation
    q_file << q_json.dump(4);

    t_file.close();
    q_file.close();
}

std::tuple<std::vector<std::vector<std::pair<double, double>>>,
           std::vector<std::vector<std::pair<double, double>>>>
loadTables() {
    std::ifstream t_file("3w_data/origin/t_table_chengdu.json");
    std::ifstream q_file("3w_data/origin/q_table_chengdu.json");

    std::vector<std::vector<std::pair<double, double>>> t_table_loaded, q_table_loaded;

    if (t_file.is_open() && q_file.is_open()) {
        json t_json, q_json;
        t_file >> t_json;
        q_file >> q_json;

        t_table_loaded = t_json.get<std::vector<std::vector<std::pair<double, double>>>>();
        q_table_loaded = q_json.get<std::vector<std::vector<std::pair<double, double>>>>();
    } else {
        // Handle case where files don't exist or can't be opened
        // You might want to re-fetch or generate the tables here
        std::cerr << "Files t_table.json or q_table.json not found. Fetching data...\n";
    }

    t_file.close();
    q_file.close();

    return std::make_tuple(t_table_loaded, q_table_loaded);
}




void scale_saveTables(const std::vector<std::vector<std::pair<double,double>>>& t_table)
{
    json t_json(t_table);

    std::ofstream t_file("scale_data/origin/q_table_chengdu_reverse.json");

    t_file << t_json.dump(4); // pretty-print with 4 spaces of indentation

    t_file.close();
}

std::tuple<std::vector<std::vector<std::pair<double, double>>>>
scale_loadTables() {
    std::ifstream t_file("scale_data/normal/t_table_chengdu.json");

    std::vector<std::vector<std::pair<double, double>>> t_table_loaded;

    if (t_file.is_open()) {
        json t_json;
        t_file >> t_json;

        t_table_loaded = t_json.get<std::vector<std::vector<std::pair<double, double>>>>();
    } else {
        // Handle case where files don't exist or can't be opened
        // You might want to re-fetch or generate the tables here
        std::cerr << "Files t_table.json or q_table.json not found. Fetching data...\n";
    }

    t_file.close();

    return std::make_tuple(t_table_loaded);
}




struct DistandIndex {
    double dist;
    int index;
    bool operator<(const DistandIndex& other) const {
        return dist < other.dist;
    }
};


// 计算并存储groundtruth结果到文件
void computeAndStoreGroundtruth(const std::vector<path>& query_trajectories,
                                const std::vector<path>& data_trajectories,
                                const std::string path)
{

    auto start_build=std::chrono::high_resolution_clock::now();

    std::ofstream outfile(path); // 指定文件路径，当前目录下
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }

    #pragma omp parallel for
    for (int i = 0; i < query_trajectories.size(); ++i) 
    {
        std::vector<DistandIndex> dists_and_indices_queryies_groundtruth(data_trajectories.size());

        #pragma omp parallel for
        for (int j = 0; j < data_trajectories.size(); ++j) 
        {
            subResult result = execute("pss",query_trajectories[i], data_trajectories[j]);
            dists_and_indices_queryies_groundtruth[j] = { result.second, j };
            //std::cout<<"query与第 "<<j<<" 条data轨迹的groundtruth"<<" start from "<<result.first.first<<" to "<<result.first.second<<", distances "<<result.second<<std::endl;
        }

        std::sort(dists_and_indices_queryies_groundtruth.begin(), dists_and_indices_queryies_groundtruth.end());

        
        for (int t = 0; t < dists_and_indices_queryies_groundtruth.size(); ++t) {
            outfile << dists_and_indices_queryies_groundtruth[t].index << " ";
        }
        outfile << std::endl;
        
        //std::cout<<"第"<<i<<"条轨迹的groundtruth"<<std::endl;

        //std::cout<<"这条query的groundtruth ";
        //for (int t = 0; t < 50; ++t) {
        //    std::cout<<dists_and_indices_queryies_groundtruth[t].index<<" ";
        //}
        //std::cout<<std::endl;
    }

    outfile.close();

    auto end_build=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_build=end_build-start_build;
    std::cout<<"Elpased time: "<<elapsed_build.count()<<" seconds\n"<<endl;
}

// 从文件中读取并比较groundtruth结果
void readAndCompareGroundtruth(const std::string path, const int knn, std::vector<std::vector<int>>& indices_queryies) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::ifstream infile(path); // 指定文件路径，当前目录下
    if(!infile.is_open())
    {
        std::cerr<<"Error: Unable to open file for reading."<<std::endl;
        return;
    }

    double count_HR10_sum=0;
    double count_HR20_sum=0;
    double count_HR30_sum=0;
    double count_HR40_sum=0;
    double count_HR50_sum=0;

    for(int i=0;i<indices_queryies.size();++i)
    {
        std::vector<int> result_indices(knn);
        int index;

        for(int t=0;t<knn;++t)
        {
            if(!(infile>>result_indices[t]))
            {
                std::cerr<<"Error: Unable to read expected number of results from file."<<std::endl;
                infile.close();
                return;
            }
        }

        int HR10_size=static_cast<int>(0.2*(indices_queryies[i].size()));
        std::unordered_set<int> set1(indices_queryies[i].begin(),indices_queryies[i].begin()+HR10_size);

        int HR20_size=static_cast<int>(0.4*(indices_queryies[i].size()));
        std::unordered_set<int> set2(indices_queryies[i].begin(),indices_queryies[i].begin()+HR20_size);

        int HR30_size=static_cast<int>(0.6*(indices_queryies[i].size()));
        std::unordered_set<int> set3(indices_queryies[i].begin(),indices_queryies[i].begin()+HR30_size);

        int HR40_size=static_cast<int>(0.8*(indices_queryies[i].size()));
        std::unordered_set<int> set4(indices_queryies[i].begin(),indices_queryies[i].begin()+HR40_size);

        int HR50_size=static_cast<int>(1.0*(indices_queryies[i].size()));
        std::unordered_set<int> set5(indices_queryies[i].begin(),indices_queryies[i].begin()+HR50_size);

        std::cout<<"第 "<<i<<" 条query执行groundtruth查询完成"<<std::endl;
        int count_HR10=0;
        for(int t_HR10=0;t_HR10<int(0.2*knn);++t_HR10)
        {
            if(set1.count(result_indices[t_HR10]))
            {
                count_HR10++;
            }
        }
        std::cout<<"HR10 "<<(count_HR10/(0.2*knn))<<std::endl;
        count_HR10_sum+=(count_HR10/(0.2*knn));

        int count_HR20=0;
        for(int t_HR20=0;t_HR20<int(0.4*knn);++t_HR20)
        {
            if(set2.count(result_indices[t_HR20]))
            {
                count_HR20++;
            }
        }
        std::cout<<"HR20 "<<(count_HR20/(0.4*knn))<<std::endl;
        count_HR20_sum+=(count_HR20/(0.4*knn));

        int count_HR30=0;
        for(int t_HR30=0;t_HR30<int(0.6*knn);++t_HR30)
        {
            if(set3.count(result_indices[t_HR30]))
            {
                count_HR30++;
            }
        }
        std::cout<<"HR30 "<<(count_HR30/(0.6*knn))<<std::endl;
        count_HR30_sum+=(count_HR30/(0.6*knn));

        int count_HR40=0;
        for(int t_HR40=0;t_HR40<int(0.8*knn);++t_HR40)
        {
            if(set4.count(result_indices[t_HR40]))
            {
                count_HR40++;
            }
        }
        std::cout<<"HR40 "<<(count_HR40/(0.8*knn))<<std::endl;
        count_HR40_sum+=(count_HR40/(0.8*knn));

        int count_HR50=0;
        for(int t_HR50=0;t_HR50<int(1.0*knn);++t_HR50)
        {
            if(set5.count(result_indices[t_HR50]))
            {
                count_HR50++;
            }
        }
        std::cout<<"HR50 "<<(count_HR50/(1.0*knn))<<std::endl;
        count_HR50_sum+=(count_HR50/(1.0*knn));

        std::string dummy;
        getline(infile,dummy);

    }
    infile.close();
    
    std::cout<<"average HR10 "<<count_HR10_sum/indices_queryies.size()<<std::endl;
    std::cout<<"average HR20 "<<count_HR20_sum/indices_queryies.size()<<std::endl;
    std::cout<<"average HR30 "<<count_HR30_sum/indices_queryies.size()<<std::endl;
    std::cout<<"average HR40 "<<count_HR40_sum/indices_queryies.size()<<std::endl;
    std::cout<<"average HR50 "<<count_HR50_sum/indices_queryies.size()<<std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    std::ifstream infile_(path); // 指定文件路径，当前目录下
    if (!infile_.is_open())
    {
        std::cerr << "Error: Unable to open file for reading." << std::endl;
        return;
    }

    double count_hr20_sum = 0;
    double count_r10_sum = 0;
    double count_hr50_sum = 0;

    for (int i = 0; i < indices_queryies.size(); ++i)
    {
        std::vector<int> result_indices(knn);
        int index;

        for (int t = 0; t < knn; ++t) 
        {
            if (!(infile_ >> result_indices[t])) 
            {
                std::cerr << "Error: Unable to read expected number of results from file." << std::endl;
                infile_.close();
                return;
            }
        }

        int hr20_size=static_cast<int>(0.4*(indices_queryies[i].size()));
        std::unordered_set<int> set1(indices_queryies[i].begin(),indices_queryies[i].begin()+hr20_size);

        std::unordered_set<int> set2(indices_queryies[i].begin(),indices_queryies[i].end());
        std::cout<<"第 "<<i<<" 条query执行groundtruth查询完成"<<std::endl;
        int count_hr20 = 0;
        for(int t_hr20=0; t_hr20 < int(0.4*knn); ++t_hr20)
        {
            if(set1.count(result_indices[t_hr20]))
            {
                count_hr20++;
            }
        }
        std::cout<<"HR20 "<<(count_hr20/(0.4*knn))<<std::endl;
        count_hr20_sum += (count_hr20/(0.4*knn));

        int count_r10 = 0;
        for (int t_r10 = 0; t_r10 < 0.2 * knn; ++t_r10) 
        {
            if (set2.count(result_indices[t_r10])) 
            {
                count_r10 += 1;
            }
        }
        std::cout<<"R10_50 "<<(count_r10/(0.2*knn))<<std::endl;
        count_r10_sum += (count_r10 / (0.2 * knn));

        int count_hr50 = 0;
        for (int t_hr50 = 0; t_hr50 < knn; ++t_hr50) 
        {
            if (set2.count(result_indices[t_hr50])) 
            {
                count_hr50++;
            }
        }
        std::cout<<"HR50 "<<(count_hr50/(1.0*knn))<<std::endl;
        count_hr50_sum += (count_hr50 / (1.0 * knn));

        // std::cout<<"真正的top-k: "<<std::endl;
        // for(int j=0;j<result_indices.size();j++)
        // {
        //     std::cout<<result_indices[j]<<" ";
        // }
        // std::cout<<endl;

        // std::cout<<"预测的top-k: "<<std::endl;
        // for(int j=0;j<indices_queryies[i].size();j++)
        // {
        //     std::cout<<indices_queryies[i][j]<<" ";
        // }
        // std::cout<<endl;
        
        std::string dummy;
        getline(infile_,dummy);
    }

    infile_.close();

    std::cout << "average HR20 " << count_hr20_sum / indices_queryies.size() << std::endl;
    std::cout << "average R10_50 " << count_r10_sum / indices_queryies.size() << std::endl;
    std::cout << "average HR50 " << count_hr50_sum / indices_queryies.size() << std::endl;
}

int main()
{

    signal(SIGTERM,signalHandler);
    std::cout<<"hello,world"<<endl;
    
     Config config {
        "/home/KNN/data/cp_data/chengdu",  // data_path
        30,                                 // query_min_len
        89,                                // query_max_len
        90,                                 // data_min_len
        300,                                // data_max_len
        "dtw",                              // metric
        "subtra"                               // task
    };

    std::vector<std::tuple<int, int, int, int, double>> score_table;
    std::vector<std::vector<std::pair<double, double>>> t_table;
    std::vector<std::vector<std::pair<double, double>>> q_table;

    //py::scoped_interpreter guard{}; // Start the interpreter



//数据处理和存储的注释**************************************************************************************//
    // try {
    //     // Import the Python module and function
    //     py::module data_loader = py::module::import("get_data");
    //     py::object get_data_func_score = data_loader.attr("get_data_score");
    //     py::object get_data_func_t = data_loader.attr("get_data_t");
    //     py::object get_data_func_q = data_loader.attr("get_data_q");

    //     // Convert Config struct to Python dict
    //     py::dict config_dict;
    //     config_dict["data_path"] = config.data_path;
    //     config_dict["query_min_len"] = config.query_min_len;
    //     config_dict["query_max_len"] = config.query_max_len;
    //     config_dict["data_min_len"] = config.data_min_len;
    //     config_dict["data_max_len"] = config.data_max_len;
    //     config_dict["metric"] = config.metric;
    //     config_dict["task"] = config.task;

    //     // Call Python function and get returned tuple
    //     auto result_tuple_score = get_data_func_score(config_dict);
    //     auto result_tuple_t = get_data_func_t(config_dict);
    //     auto result_tuple_q = get_data_func_q(config_dict);

    //     std::cout<<"执行成功"<<endl;
    //     // Extract elements from Python tuple to C++ vectors
        

    //     // Extract score_table
    //     py::object score_table_py = result_tuple_score;
    //     py::object score_table_columns=score_table_py.attr("columns");

    //     std::cout<<"执行成功"<<endl;
    //     int i=0;
    //     for (py::handle score_row : score_table_py.attr("itertuples")()) {
    //         score_table.push_back(std::make_tuple(
    //             score_row.attr("que_index").cast<int>(),
    //             score_row.attr("tra_idx").cast<int>(),
    //             score_row.attr("s").cast<int>(),
    //             score_row.attr("e").cast<int>(),
    //             score_row.attr("score1").cast<double>()
    //         ));
    //         if(i==0)
    //         {
    //             std::cout<<score_row<<endl;
    //             i++;
    //         }
    //     }
    //     //std::cout<<"执行成功"<<endl;

    //     for (const auto& sublist : result_tuple_t) 
    //     {
    //         std::vector<std::pair<double, double>> inner;
    //         for (const auto& tup : sublist) 
    //         {
    //             inner.emplace_back(py::cast<std::pair<double, double>>(tup));
    //         }

    //         t_table.push_back(inner);
    //         //std::cout<<"执行成功"<<endl;
    //     }

    //     for (const auto& sublist : result_tuple_q) 
    //     {
    //         std::vector<std::pair<double, double>> inner;
    //         for (const auto& tup : sublist) 
    //         {
    //             auto lat_lon=py::cast<std::pair<double,double>>(tup);
    //             inner.emplace_back(lat_lon.second,lat_lon.first);

    //             //inner.emplace_back(py::cast<std::pair<double, double>>(tup));
    //         }
    //         q_table.push_back(inner);
    //         //std::cout<<"执行成功"<<endl;
    //     }
        
       
    //     std::cout<<"执行成功"<<endl;
    // }
    // catch(const py::error_already_set& e)
    // {
    //     std::cerr<<"Python error:\n"<<e.what()<<std::endl;
    // }
    // exit(1);

    // Print some data to verify

    //saveTables(t_table,q_table);
    //scale_saveTables(q_table);
//数据存储的注释**************************************************************************************//
    std::cout<<"开始读取"<<std::endl;

    
    std::tie(t_table,q_table)=loadTables();
    //std::tie(t_table)=scale_loadTables();
    std::cout<<"读取结束"<<std::endl;
    std::cout<<t_table.size()<<std::endl;
    std::cout<<q_table.size()<<std::endl;
    for(int i=0;i<1;i++)
    {
        auto trajectory_t=t_table[i];
        for(int j=0;j<1;j++)
        {
            auto point_t=trajectory_t[j];
            std::cout<<std::get<0>(point_t)<<" "<<std::get<1>(point_t)<<std::endl;
        }
    }
    for(int i=0;i<1;i++)
    {
        auto trajectory_q=q_table[i];
        for(int j=0;j<1;j++)
        {
            auto point_q=trajectory_q[j];
            std::cout<<std::get<0>(point_q)<<" "<<std::get<1>(point_q)<<std::endl;
        }
    }
    


    std::cout << "Score Table:" << std::endl;
    std::cout<<"Loaded trajectories: "<<t_table.size()<<std::endl;
    std::cout<<q_table.size()<<std::endl;
    // for (auto it = t_table.begin(); it != t_table.begin() + 3&& it != t_table.end(); ++it) 
    // {
    //     const auto& vec = *it;
    //     path ls;
        
    //     for (const auto& pair : vec) {
    //         std::cout<<pair.first<<" "<<pair.second<<std::endl;
    //     }
    // }

    //size_t num_trajectories_to_visualize=100;
    //visualize_trajectories(t_table,num_trajectories_to_visualize);   

    //std::cout<<paths.size()<<endl;
    //std::cout<<paths[0].size()<<endl;
    //std::cout<<paths[1].size()<<endl;
    //std::cout<<t_table.size()<<endl;

    //exit(1);

    // std::vector<path> data_trajectories(t_table.begin(),t_table.begin()+100);
    // 重构赋值data_trajectories的方式
    std::vector<path> data_trajectories;
    data_trajectories.reserve(30000); 

    // Iterate through t_table and convert each path to a linestring
    for (auto it = t_table.begin(); it != t_table.begin() + 30000 && it != t_table.end(); ++it) 
    {
        const auto& vec = *it;
        path ls;
        
        for (const auto& pair : vec) {
            rtree_point pt(pair.first, pair.second);
            bg::append(ls, pt);
        }
        data_trajectories.push_back(ls);
    }
    
    
    //std::vector<path> first_10_paths(10);
    //std::cout<<"执行了吗"<<endl;
    //std::copy_n(t_table.begin(),10,std::back_inserter(first_10_paths));
    //std::cout<<"执行了吗"<<endl;
    std::cout<<data_trajectories.size()<<endl;
    //std::cout<<data_trajectories[0].size()<<endl;
    //std::cout<<data_trajectories[1].size()<<endl;
    //exit(1);



    //构建gnns index
    //Gnns_Index<L2Distance<float>> gnns_index(3000);
    Gnns_Index<L2Distance<float>> gnns_index(data_trajectories);

    std::cout<<"执行成功"<<endl;
    //vector<int> querys{649};

    // build the gnns index
    //Gnns_Index<L2Distance<float> > gnns_index(data);
    //最后一个参数是rebuild，为true表示从头开始执行 knn_graph，计算data-to-data之间的top-k，为false表示从文件中读取之前已经计算好的knn_graph，
    auto start_build=std::chrono::high_resolution_clock::now();
    //gnns_index.build_index(algorithm, graph_index_saved_path, graph_dist_saved_path, top_level_graph_points_saved_path, top_level_graph_indices_saved_path, top_level_graph_distances_saved_path, true);
    //分类
    gnns_index.build_index(algorithm,graph_index_saved_path,graph_dist_saved_path,graph_start_end_position_saved_path,top_level_graph_points_saved_path,top_level_graph_indices_saved_path,top_level_graph_distances_saved_path,top_level_graph_start_end_position_saved_path,false);

    auto end_build=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_build=end_build-start_build;
    std::cout<<"Elpased time: "<<elapsed_build.count()<<" seconds\n"<<endl;
    std::cout<<"knn_build执行完成"<<endl;
    std::cout<<"创建成功"<<endl;
    //exit(1);



    // std::vector<path> query_trajectories(q_table.begin(),q_table.begin()+2000);
    // 重构赋值data_trajectories的方式
    
    std::vector<path> query_trajectories;
    
    //根据索引提取查询轨迹
    //xian
    //std::vector<int> indices_to_extract={1955};
    //chengdu
    //std::vector<int> indices_to_extract={1955, 1561, 821, 1129, 1810, 1320, 1473, 1668, 1446, 642, 1771, 103, 673, 774, 36, 1301, 
    //                                     638, 931, 1578, 1884, 953, 573, 1283, 656, 1455, 23, 1040, 441, 732, 1060, 749, 818, 1277, 
    //                                     950, 1294, 1547, 1042, 1300, 192, 729};
    //query_trajectories.reserve(indices_to_extract.size());
    // for(const auto& index: indices_to_extract)
    // {
    //     if(index<q_table.size())
    //     {
    //         const auto& vec=q_table[index];
    //         path ls;
            
    //         for(const auto& pair: vec)
    //         {
    //             rtree_point pt(pair.first, pair.second);
    //             bg::append(ls, pt);
    //         }
    //         query_trajectories.push_back(ls);
    //     }
    //     else 
    //     {
    //         std::cout<<"Index "<<index<<" is out of range."<<std::endl;
    //     }
    // }
    // std::cout<<query_trajectories.size()<<std::endl;


    query_trajectories.reserve(8); 

    // Iterate through t_table and convert each path to a linestring
    for (auto it = q_table.begin(); it != q_table.begin() + 8 && it != q_table.end(); ++it) 
    {
        const auto& vec = *it;
        path ls;
        
        for (const auto& pair : vec) {
            rtree_point pt(pair.first, pair.second);
            bg::append(ls, pt);
        }
        query_trajectories.push_back(ls);
    }
    std::cout<<query_trajectories.size()<<endl;


    std::cout<<"开始执行knn_search"<<endl;
    std::cout<<query_trajectories.size()<<endl;
    // //std::cout<<query_trajectories[0].size()<<endl;
    // //std::cout<<query_trajectories[1].size()<<endl;

    



    //knn搜索
    int knn=50;
    //获取query轨迹

    //存储query轨迹邻居的indexs
    int queryies_num=query_trajectories.size();
    
    std::vector<std::vector<int>>indices_queryies(queryies_num,std::vector<int>(knn));
    std::vector<std::vector<float>>dists_queryies(queryies_num,std::vector<float>(knn));

    auto start_search=std::chrono::high_resolution_clock::now();
    gnns_index.knn_search(query_trajectories,indices_queryies,dists_queryies,knn,Search_Params());
    auto end_search=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_search=end_search-start_search;
    std::cout<<"Elpased time: "<<elapsed_search.count()<<" seconds\n"<<endl;
    std::cout<<"knn_search执行完成"<<endl;

    //write_results_to_file(top_k_results_to_RR_saved_path,indices_queryies,dists_queryies);
    //std::cout<<"top-k结果已保存到: "<<top_k_results_to_RR_saved_path<<std::endl;




    //std::vector<std::vector<int>> indices_queryies_groundtruth(queryies_num,std::vector<int>(knn));
    //std::vector<std::vector<float>>dists_queryies_groundtruth(queryies_num,std::vector<float>(knn));


    //computeAndStoreGroundtruth(query_trajectories,data_trajectories,indices_queryies_groundtruth_saved_path);

    readAndCompareGroundtruth(indices_queryies_groundtruth_saved_path,knn,indices_queryies);
}