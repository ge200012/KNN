#include "IndexFilter.h"
#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <numeric>

//std::map<std::pair<int, int>, std::set<int>> invert_index;

class ThePoint
{
    public:
        vector<int> N;
        int C;
        double w;
        double v;

        ThePoint(vector<int> n, int c):N(n),C(c),w(0),v(0){}
};

std::pair<std::set<int>, std::map<std::pair<int, int>, std::set<int>>> run_get_index(const std::vector<path>& paths, const std::vector<value>& result_grid, double grid_size) 
{
    std::map<std::pair<int, int>, std::set<int>> invert_index;
    std::set<int> selected_ids;
    for (size_t i = 0; i < result_grid.size(); ++i) 
    {
        for (const auto& point : paths[result_grid[i].second]) 
        {
            int grid_x=static_cast<int>((point.get<0>() - 30.72693)/grid_size);
            int grid_y=static_cast<int>((point.get<1>() - 104.04761)/grid_size);

            std::pair<int,int> val={grid_x,grid_y};

            invert_index[val].insert(result_grid[i].second);
        }
        selected_ids.insert(result_grid[i].second);
    }
    return {selected_ids, invert_index};
}
// // 添加索引
// void add_index(const path& trajectory, int idx, double grid_size) 
// {
//     for (const auto& point : trajectory) 
//     {
//         int grid_x=static_cast<int>((point.get<0>() - 30.72693)/grid_size);
//         int grid_y=static_cast<int>((point.get<1>() - 104.04761)/grid_size);

//         std::pair<int,int> val={grid_x,grid_y};

//         invert_index[val].insert(idx);
//     }
// }

// // 运行索引生成
// std::set<int> run_get_index(const std::vector<path>& paths, const std::vector<value>& result_grid, double grid_size) 
// {
//     std::set<int> selected_ids;
//     for (size_t i = 0; i < result_grid.size(); ++i) 
//     {
//         add_index(paths[result_grid[i].second], result_grid[i].second, grid_size);
//         selected_ids.insert(result_grid[i].second);
//     }
//     return selected_ids;
// }
// std::set<int> run_get_index(const std::vector<path>& paths, double grid_size) 
// {
//     for (size_t i = 0; i < paths.size(); ++i) 
//     {
//         add_index(paths[i], i, grid_size);
//     }
// }

// 获取候选项
std::pair<std::set<int>, double> get_candi_osf(
    const path& current_trajectory,
    const std::map<std::pair<int, int>, std::set<int>> invert_index,
    double rate=0.1,
    double grid_size=0.0005,
    const std::vector<path>& paths={},
    const std::vector<value>& result_grid={},
    const std::string& dataset="",
    const std::string& metric="") 
{
    double eta=0.0;
    if(dataset=="chengdu")
    {
        eta=0.000242;
    }
    else if(dataset=="xian")
    {
        eta=0.00023;
    }
    else if(dataset=="porto")
    {
        eta=0.00079;
    }
    
    if(metric=="edr"||metric=="lcss")
    {
        eta=0.001;
        if(dataset=="porto")
        {
            eta=0.003;
        }
    }
    //std::cout<<"执行到5"<<std::endl;

    std::set<int> candis;
    std::vector<std::pair<std::vector<int>,double>> candi_list;
    std::pair<int,int> last_val={-1,-1};
    double temp_time=0.0;

    
    //使用result_grid作为过滤条件
    std::set<int> valid_indices; //用于存储result_grid中的轨迹索引
    for(const auto& rg:result_grid)
    {
        valid_indices.insert(rg.second);
    }


    // 示例逻辑，替换为实际计算
    for (const auto& point : current_trajectory) 
    {
        int grid_x=static_cast<int>((point.get<0>()-30.72693)/grid_size);
        int grid_y=static_cast<int>((point.get<1>()-104.04761)/grid_size);
        auto start_time=std::chrono::steady_clock::now();
        std::set<int> candi;
        //std::cout<<"执行到6"<<std::endl;

        auto val=std::make_pair(grid_x,grid_y);
        if(last_val==val)
        {
            continue;
        }
        else 
        {
            last_val=val;
            for(const auto& j: {
                    std::make_pair(grid_x - 1, grid_y - 1),
                    std::make_pair(grid_x, grid_y - 1),
                    std::make_pair(grid_x + 1, grid_y - 1),
                    std::make_pair(grid_x - 1, grid_y),
                    std::make_pair(grid_x, grid_y),
                    std::make_pair(grid_x + 1, grid_y),
                    std::make_pair(grid_x - 1, grid_y + 1),
                    std::make_pair(grid_x, grid_y + 1),
                    std::make_pair(grid_x + 1, grid_y + 1)
                })
                {
                    auto it=invert_index.find(j);
                    if(it!=invert_index.end())
                    {
                        candi.insert(it->second.begin(),it->second.end());
                    }
                }
        }
        //std::cout<<"执行到7"<<std::endl;

        auto end_time=std::chrono::steady_clock::now();
        temp_time+=std::chrono::duration<double>(end_time-start_time).count()/9;

        std::vector<int> f_candi;
        double m_dis=0.01;
        for(int k:candi)
        {


            if(valid_indices.find(k)==valid_indices.end())
            {
                continue;
            }


            const auto& td=paths[k];
            for(const auto& td_p:td)
            {
                double dis_p=std::sqrt(std::pow(point.get<0>()-td_p.get<0>(),2)+std::pow(point.get<1>()-td_p.get<1>(),2));
                if(dis_p<=eta)
                {
                    f_candi.push_back(k);
                    break;
                }
                else if(dis_p<m_dis)
                {
                    m_dis=dis_p;
                }
            }
        }
        candi_list.emplace_back(f_candi,m_dis);
    }
    //std::cout<<"执行到8"<<std::endl;

    double cq = 0;
    double all_cost = (metric == "edr" || metric == "lcss") ?
        (candi_list.size() * rate) :
        (std::accumulate(candi_list.begin(), candi_list.end(), 0.0,
            [](double sum, const auto& c) { return sum + c.second; }) * rate);
    //std::cout<<"执行到9"<<std::endl;

    std::vector<ThePoint> all_point;
    for (const auto& c : candi_list) 
    {
        all_point.emplace_back(c.first, c.second); // 需要根据实际情况设置初始化参数
    }
    //std::cout<<"执行到10"<<std::endl;

    auto start_time = std::chrono::steady_clock::now();
    while (cq < all_cost) 
    {
        //std::cout<<"执行到12"<<std::endl;
        for (auto& p : all_point)
        {
            p.v = (p.N.size() - p.w) / std::min(static_cast<double>(p.C), all_cost - cq);
        }
        //std::cout<<"执行到13"<<std::endl;
        std::sort(all_point.begin(), all_point.end(), [](const ThePoint& a, const ThePoint& b) { return a.v < b.v; });
        //std::cout<<"这里可以执行"<<std::endl;

        if(all_point.empty())
        {
            break;
        }
        ThePoint the_candi = all_point.front();
        //std::cout<<all_point.size()<<std::endl;
        all_point.erase(all_point.begin());
        //std::cout<<all_point.size()<<std::endl;
        //std::cout<<"执行到14"<<std::endl;

        for (auto& p : all_point) 
        {
            p.w += std::min(static_cast<double>(p.C), all_cost - cq) * the_candi.v;
        }
        //std::cout<<"执行到15"<<std::endl;

        cq += the_candi.C;
        candis.insert(the_candi.N.begin(), the_candi.N.end());
        //std::cout<<"执行到16"<<std::endl;
    }
    //std::cout<<"执行到11"<<std::endl;

    auto end_time = std::chrono::steady_clock::now();
    return {candis, temp_time + std::chrono::duration<double>(end_time - start_time).count()};
}

// OSF过滤
std::pair<std::set<int>, double> osf_filter(
    const path& current_trajectory,
    const std::map<std::pair<int, int>, std::set<int>> invert_index,
    double radius,
    double grid_size,
    const std::vector<path>& paths,
    const std::vector<value>& result_grid,
    const std::string& dataset="",
    const std::string& metric="") 
{
    auto candi = get_candi_osf(current_trajectory, invert_index, radius, grid_size, paths, result_grid, dataset, metric);
    return candi;
}
