#ifndef GNNS_KNN_GRAPH_H
#define GNNS_KNN_GRAPH_H

#include <iostream>
#include <string>
#include <exception>
#include <algorithm>
#include <vector>
#include "../util/matrix.h"
#include "../io.h"
#include "../util/dist.h"
#include "../define.h"
#include "../general.h"


//自己添加的内容
#include "../similarTrajectory.h"
#include "mostSimilar.h"
#include <unordered_set>
#include "io.h"
#include <random>
#include <chrono>
#include <omp.h>
#include "constant.h"
#include "distance.h"
#include "IndexFilter.h"

//namespace gnns
//{
    /*
     * the build graph method
     * naive: use an brute way to build O(n*n*d)
     */
    enum BUILD_GRAPH_METHOD
    {
        NAIVE
    };

    template <typename Distance>
    class Knn_Graph
    {
    public:
        typedef typename Distance::ElementType ElementType;
        typedef typename Distance::DistanceType DistanceType;

    private:
        /*
         * @brief find the first n minimum elements and their index
         * @param dist: a vector which records the distance of each two points
         * @param length: the length of the dist vector
         * @param elements: the distance of first n points with minimum distances
         * @param indices: the index of the first n points with minimum distances
         * @param n: the param n
         */
        //top-level-graph选择上层graph的节点集合
        void top_level_points_build(std::vector<std::vector<double>>& dist,const int length,int drop_neighbors)
        {
            //用两个结构是为了空间换时间，points_visitpool是为了在加最远点时判断是否已被visited更快，points_set是为了在while循环时判断是否还有点可加时更快，
            std::vector<char> top_level_points_visitpool(length,false);
            
            std::unordered_set<int> top_level_points_set;
            for(int i=0;i<length;i++)
            {
                top_level_points_set.insert(i);
            }
            std::cout<<length<<endl;

            int add_top_level_point=rand()%length;
            std::cout<<"来到这里"<<endl;
            top_level_points.push_back(add_top_level_point);

            //因为是开始，所以加了这两句
            top_level_points_visitpool[add_top_level_point]=true;
            top_level_points_set.erase(add_top_level_point);

            //std::cout<<"top_level_points_set的size "<<top_level_points_set.size()<<endl;
            while(!top_level_points_set.empty())
            {
                //std::cout<<"add_top_level_point为 "<<add_top_level_point<<endl;

                //这里注释掉是因为可以在下面for循环中删除邻居时删除自己
                top_level_points_visitpool[add_top_level_point]=true;
                top_level_points_set.erase(add_top_level_point);

                std::vector<std::pair<double,int>> top_level_v;

                for(int i=0;i<length;++i)
                {
                    top_level_v.push_back(std::make_pair(dist[add_top_level_point][i],i));
                }
                //std::cout<<"第一个for循环"<<endl;
                std::nth_element(top_level_v.begin(),top_level_v.begin()+drop_neighbors,top_level_v.end());
                for(auto it=top_level_v.begin();it!=top_level_v.begin()+drop_neighbors;++it)
                {
                    top_level_points_visitpool[it->second]=true;
                    top_level_points_set.erase(it->second);
                }

                //std::cout<<"第二个for循环"<<endl;
                std::sort(top_level_v.begin()+drop_neighbors,top_level_v.end());
                for(auto it=top_level_v.end()-1;it!=top_level_v.begin()+drop_neighbors-1;--it)
                {
                    if(top_level_points_visitpool[it->second]==false)
                    {
                        add_top_level_point=it->second;
                        top_level_points.push_back(add_top_level_point);
                        break;
                    }
                }
                //std::cout<<"第二个for循环"<<endl;
                //std::cout<<"top_level_points_set的size "<<top_level_points_set.size()<<endl;
            }
        }


        // 自动计算网格大小
        std::tuple<double,double,double,double,std::pair<double,double>> calculate_grid_size(const std::vector<path>& trajectories, int num_grids) 
        {
            // 计算轨迹数据集的范围
            double min_x = std::numeric_limits<double>::max();
            double min_y = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double max_y = std::numeric_limits<double>::min();

            for (const auto& traj : trajectories) 
            {
                for (const auto& point : traj) 
                {
                    double x = bg::get<0>(point);
                    double y = bg::get<1>(point);
                    if (x < min_x) min_x = x;
                    if (y < min_y) min_y = y;
                    if (x > max_x) max_x = x;
                    if (y > max_y) max_y = y;
                }
            }

            double graph_width=max_x-min_x;
            double graph_height=max_y-min_y;
            double grid_width=graph_width/num_grids;
            double grid_height=graph_height/num_grids;
            std::pair<double,double> grid_size(grid_width,grid_width);
            return std::make_tuple(min_x,min_y,max_x,max_y,grid_size);
        }


        //将轨迹点转换为网格索引
        rtree_point to_centpoint_grid_index(const rtree_point& p, double grid_size) 
        {
            int x_index = static_cast<int>(std::floor(bg::get<0>(p) / grid_size));
            int y_index = static_cast<int>(std::floor(bg::get<1>(p) / grid_size));
            return rtree_point(x_index, y_index);
        }


        //计算轨迹点的grid_size大小周边网格
        bg::model::box<rtree_point> to_centpoint_grid_range(const rtree_point& p,double grid_size)
        {
            double half_grid=grid_size/2.0f;
            rtree_point min_corner(bg::get<0>(p)-half_grid,bg::get<1>(p)-half_grid);
            rtree_point max_corner(bg::get<0>(p)+half_grid,bg::get<1>(p)+half_grid);
            return bg::model::box<rtree_point>(min_corner,max_corner);
        }


        bg::model::box<rtree_point> to_traj_grid_range(const path& trajectory)
        {
            // 计算轨迹数据集的范围
            double min_x = std::numeric_limits<double>::max();
            double min_y = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double max_y = std::numeric_limits<double>::min();

            for (const auto& point : trajectory) 
            {
                double x = bg::get<0>(point);
                double y = bg::get<1>(point);
                if (x < min_x) min_x = x;
                if (y < min_y) min_y = y;
                if (x > max_x) max_x = x;
                if (y > max_y) max_y = y;
            }
            
            rtree_point min_corner(min_x,min_y);
            rtree_point max_corner(max_x,max_y);
            return bg::model::box<rtree_point>(min_corner,max_corner);
        }


        bg::model::box<rtree_point> to_traj_Mul_grid_range(const path& trajectory, const double Mul_grid)
        {
            // 计算轨迹数据集的范围
            double min_x = std::numeric_limits<double>::max();
            double min_y = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double max_y = std::numeric_limits<double>::min();

            for (const auto& point : trajectory) 
            {
                double x = bg::get<0>(point);
                double y = bg::get<1>(point);
                if (x < min_x) min_x = x;
                if (y < min_y) min_y = y;
                if (x > max_x) max_x = x;
                if (y > max_y) max_y = y;
            }
            
            // 计算当前范围的中心点
            double center_x = (min_x + max_x) / 2.0;
            double center_y = (min_y + max_y) / 2.0;

            // 计算新的最小和最大角点，增大两倍
            double new_min_x = center_x - (center_x - min_x) * Mul_grid;
            double new_min_y = center_y - (center_y - min_y) * Mul_grid;
            double new_max_x = center_x + (max_x - center_x) * Mul_grid;
            double new_max_y = center_y + (max_y - center_y) * Mul_grid;

            // 构造并返回新的 box<rtree_point>
            rtree_point new_min_corner(new_min_x, new_min_y);
            rtree_point new_max_corner(new_max_x, new_max_y);
            return bg::model::box<rtree_point>(new_min_corner, new_max_corner);
        }


        bg::model::box<rtree_point> to_traj_far_grid_range(const path& trajectory,const std::tuple<double,double,double,double,double> total_grid)
        {
            // 计算轨迹数据集的范围
            double min_x = std::numeric_limits<double>::max();
            double min_y = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double max_y = std::numeric_limits<double>::min();

            for (const auto& point : trajectory) 
            {
                double x = bg::get<0>(point);
                double y = bg::get<1>(point);
                if (x < min_x) min_x = x;
                if (y < min_y) min_y = y;
                if (x > max_x) max_x = x;
                if (y > max_y) max_y = y;
            }
            
            double center_x=(min_x+max_x)/2.0;
            double center_y=(min_y+max_y)/2.0;

            //找最远grid
            //左下角
            double left_bottom_x=std::get<0>(total_grid)+(std::get<4>(total_grid))/2.0;
            double left_bottom_y=std::get<1>(total_grid)+(std::get<4>(total_grid))/2.0;
            //左上角
            double left_top_x=std::get<0>(total_grid)+(std::get<4>(total_grid))/2.0;
            double left_top_y=std::get<3>(total_grid)-(std::get<4>(total_grid))/2.0;
            //右上角
            double right_top_x=std::get<2>(total_grid)-(std::get<4>(total_grid))/2.0;
            double right_top_y=std::get<3>(total_grid)-(std::get<4>(total_grid))/2.0;
            //右下角
            double right_bottom_x=std::get<2>(total_grid)-(std::get<4>(total_grid))/2.0;
            double right_bottom_y=std::get<1>(total_grid)+(std::get<4>(total_grid))/2.0;


            double distances[4];
            distances[0]=std::sqrt(std::pow(center_x-left_bottom_x,2)+std::pow(center_y-left_bottom_y,2));
            distances[1]=std::sqrt(std::pow(center_x-left_top_x,2)+std::pow(center_y-left_top_y,2));
            distances[2]=std::sqrt(std::pow(center_x-right_top_x,2)+std::pow(center_y-right_top_y,2));
            distances[3]=std::sqrt(std::pow(center_x-right_bottom_x,2)+std::pow(center_y-right_bottom_y,2));

            int farthest_index=0;
            double max_distance=distances[0];
            for(int i=1;i<4;++i)
            {
                if(distances[i]>max_distance)
                {
                    max_distance=distances[i];
                    farthest_index=i;
                }
            }

            double farthest_x,farthest_y;

            rtree_point min_corner;
            rtree_point max_corner;
            switch(farthest_index)
            {
                case 0:
                    min_corner=rtree_point((left_bottom_x-std::get<4>(total_grid)/2.0,left_bottom_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((left_bottom_x+std::get<4>(total_grid)/2.0,left_bottom_y+std::get<4>(total_grid)/2.0));
                    break;
                case 1:
                    min_corner=rtree_point((left_top_x-std::get<4>(total_grid)/2.0,left_top_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((left_top_x+std::get<4>(total_grid)/2.0,left_top_y+std::get<4>(total_grid)/2.0));
                    break;
                case 2:
                    min_corner=rtree_point((right_top_x-std::get<4>(total_grid)/2.0,right_top_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((right_top_x+std::get<4>(total_grid)/2.0,right_top_y+std::get<4>(total_grid)/2.0));
                    break;
                case 3:
                    min_corner=rtree_point((right_bottom_x-std::get<4>(total_grid)/2.0,right_bottom_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((right_bottom_x+std::get<4>(total_grid)/2.0,right_bottom_y+std::get<4>(total_grid)/2.0));
                    break;
                default:
                    break;
            }
            return bg::model::box<rtree_point>(min_corner,max_corner);
        }

        bg::model::box<rtree_point> to_traj_far_2grid_range(const path& trajectory,const std::tuple<double,double,double,double,double> total_grid)
        {
            // 计算轨迹数据集的范围
            double min_x = std::numeric_limits<double>::max();
            double min_y = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double max_y = std::numeric_limits<double>::min();

            for (const auto& point : trajectory) 
            {
                double x = bg::get<0>(point);
                double y = bg::get<1>(point);
                if (x < min_x) min_x = x;
                if (y < min_y) min_y = y;
                if (x > max_x) max_x = x;
                if (y > max_y) max_y = y;
            }
            
            double center_x=(min_x+max_x)/2.0;
            double center_y=(min_y+max_y)/2.0;

            //找最远grid
            //左下角
            double left_bottom_x=std::get<0>(total_grid)+(std::get<4>(total_grid))/2.0;
            double left_bottom_y=std::get<1>(total_grid)+(std::get<4>(total_grid))/2.0;
            //左上角
            double left_top_x=std::get<0>(total_grid)+(std::get<4>(total_grid))/2.0;
            double left_top_y=std::get<3>(total_grid)-(std::get<4>(total_grid))/2.0;
            //右上角
            double right_top_x=std::get<2>(total_grid)-(std::get<4>(total_grid))/2.0;
            double right_top_y=std::get<3>(total_grid)-(std::get<4>(total_grid))/2.0;
            //右下角
            double right_bottom_x=std::get<2>(total_grid)-(std::get<4>(total_grid))/2.0;
            double right_bottom_y=std::get<1>(total_grid)+(std::get<4>(total_grid))/2.0;


            double distances[4];
            distances[0]=std::sqrt(std::pow(center_x-left_bottom_x,2)+std::pow(center_y-left_bottom_y,2));
            distances[1]=std::sqrt(std::pow(center_x-left_top_x,2)+std::pow(center_y-left_top_y,2));
            distances[2]=std::sqrt(std::pow(center_x-right_top_x,2)+std::pow(center_y-right_top_y,2));
            distances[3]=std::sqrt(std::pow(center_x-right_bottom_x,2)+std::pow(center_y-right_bottom_y,2));

            int farthest_index=0;
            double max_distance=distances[0];
            for(int i=1;i<4;++i)
            {
                if(distances[i]>max_distance)
                {
                    max_distance=distances[i];
                    farthest_index=i;
                }
            }

            double farthest_x,farthest_y;

            rtree_point min_corner;
            rtree_point max_corner;
            switch(farthest_index)
            {
                case 0:
                    min_corner=rtree_point((left_bottom_x-std::get<4>(total_grid)/2.0,left_bottom_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((left_bottom_x+std::get<4>(total_grid)/2.0*3,left_bottom_y+std::get<4>(total_grid)/2.0*3*3));
                    break;
                case 1:
                    min_corner=rtree_point((left_top_x-std::get<4>(total_grid)/2.0,left_top_y-std::get<4>(total_grid)/2.0*3));
                    max_corner=rtree_point((left_top_x+std::get<4>(total_grid)/2.0*3,left_top_y+std::get<4>(total_grid)/2.0));
                    break;
                case 2:
                    min_corner=rtree_point((right_top_x-std::get<4>(total_grid)/2.0*3,right_top_y-std::get<4>(total_grid)/2.0*3));
                    max_corner=rtree_point((right_top_x+std::get<4>(total_grid)/2.0,right_top_y+std::get<4>(total_grid)/2.0));
                    break;
                case 3:
                    min_corner=rtree_point((right_bottom_x-std::get<4>(total_grid)/2.0*3,right_bottom_y-std::get<4>(total_grid)/2.0));
                    max_corner=rtree_point((right_bottom_x+std::get<4>(total_grid)/2.0,right_bottom_y+std::get<4>(total_grid)/2.0*3));
                    break;
                default:
                    break;
            }
            return bg::model::box<rtree_point>(min_corner,max_corner);
        }
                        
                
                  

//不分类***************************************************************************************************************
       void grid_top_level_points_build(const string & algorithm, const std::vector<path>& paths,const rtree_t& rtree)
        {
            int nums_grid=30;
            auto result=calculate_grid_size(paths,nums_grid);
            double min_x,min_y,max_x,max_y;
            std::pair<double,double> grid_size;

            std::tie(min_x,min_y,max_x,max_y,grid_size)=result;

            std::unordered_set<int> top_level_points_set;
            
            //top层的轨迹索引
            for(int i=0;i<nums_grid;++i)
            {
                for(int j=0;j<nums_grid;++j)
                {
                    rtree_point min_corner;
                    rtree_point max_corner;
                    min_corner=rtree_point(min_x+i*(grid_size.first),min_y+j*(grid_size.second));
                    max_corner=rtree_point(min_x+(i+1)*(grid_size.first),min_y+(j+1)*(grid_size.second));
                    bg::model::box<rtree_point> current_grid_box(min_corner,max_corner);
                    std::vector<value> result_grid;
                    
                    // //1.nearest
                    // rtree.query(bgi::nearest(current_grid_box,1),std::back_inserter(result_grid));
                    
                    //2.within
                    rtree.query(bgi::within(current_grid_box),std::back_inserter(result_grid));
                    double Mul_grid=1.1;
                    while(result_grid.size()<1)
                    {   
                        double center_x = (bg::get<bg::min_corner,0>(current_grid_box)+bg::get<bg::max_corner,0>(current_grid_box)) / 2.0;
                        double center_y = (bg::get<bg::min_corner,1>(current_grid_box)+bg::get<bg::max_corner,1>(current_grid_box)) / 2.0;

                        // 计算新的最小和最大角点，增大两倍
                        double new_min_x = center_x - (center_x - bg::get<bg::min_corner,0>(current_grid_box)) * Mul_grid;
                        double new_min_y = center_y - (center_y - bg::get<bg::min_corner,1>(current_grid_box)) * Mul_grid;
                        double new_max_x = center_x + (bg::get<bg::max_corner,0>(current_grid_box) - center_x) * Mul_grid;
                        double new_max_y = center_y + (bg::get<bg::max_corner,1>(current_grid_box) - center_y) * Mul_grid;

                        // 构造并返回新的 box<rtree_point>
                        rtree_point new_min_corner(new_min_x, new_min_y);
                        rtree_point new_max_corner(new_max_x, new_max_y);
                        bg::model::box<rtree_point>current_grid_box(new_min_corner, new_max_corner);

                        rtree.query(bgi::within(current_grid_box),std::back_inserter(result_grid));
                        Mul_grid+=0.1;
                    }
                    
                    for(const auto& res:result_grid)
                    {
                        top_level_points_set.insert(res.second);
                        //std::cout<<"第"<<i+j+1<<"个网格的轨迹索引为"<<res.second<<std::endl;
                    }
                }
            }

            std::copy(top_level_points_set.begin(),top_level_points_set.end(),std::back_inserter(top_level_points));
            int top_paths_num=top_level_points.size();


            //分类
            std::vector<std::vector<std::pair<double, std::pair<int, std::pair<int,int>> > >> top_dist(top_paths_num,std::vector<std::pair<double, std::pair<int, std::pair<int,int>> > >(top_paths_num-1));
            for(int i=0;i<top_paths_num;++i)
            {
                int k=i;
                std::cout<<"顶层图:"<<"轨迹"<<top_level_points[i]<<"与其他轨迹计算相似得分"<<endl;
                for(int j=i+1;j<top_paths_num;++j)
                {                     
                    path v1=paths[top_level_points[i]];
                    path v2=paths[top_level_points[j]];
                    //std::cout<<"计算"<<top_level_points[i]<<"和"<<top_level_points[j]<<"之间的相似度"<<std::endl;
                    //分类

                    Result result=minSubTrajectory2(algorithm,v1,v2);
                    //std::cout<<"第"<<i<<"个轨迹从 "<<(result.second.first.first)<<" 到 "<<(result.second.first.second)<<"第"<<j<<"个轨迹从 "<<(result.second.second.first)<<" 到 "<<(result.second.second.second)<<" 代表性相似子轨迹得分是 "<<(-result.first)<<endl;
                    top_dist[i][k].first=result.first;
                    top_dist[i][k].second.first=top_level_points[j];
                    top_dist[i][k].second.second.first=result.second.first.first;
                    top_dist[i][k].second.second.second=result.second.first.second;

                    top_dist[k+1][i].first=result.first;
                    top_dist[k+1][i].second.first=top_level_points[i];
                    top_dist[k+1][i].second.second.first=result.second.second.first;
                    top_dist[k+1][i].second.second.second=result.second.second.second;
                    k++;
                }
            }
            //std::cout<<"top_level_graph相似度计算完成"<<std::endl;
            //std::cout<<"top_paths_nums为"<<top_paths_num<<std::endl;

            double portion_from_front=1.0;
            double portion_from_middle=0;
            double portion_from_back=0;

            //top_level层建立的边的个数
            int top_level_edge=static_cast<int>(0.8*top_paths_num);

            int count_from_front=static_cast<int>(top_level_edge*portion_from_front);
            int count_from_middle=static_cast<int>(top_level_edge*portion_from_middle);
            int count_from_back=static_cast<int>(top_level_edge*portion_from_back);

            //分配空间
            top_level_indices.resize(top_paths_num);
            top_level_distances.resize(top_paths_num);
            top_level_start_end_positions.resize(top_paths_num);
            for(int i=0;i<top_paths_num;++i)
            {
                //分类
                top_level_indices[i].second.resize(count_from_front+count_from_middle+count_from_back);
                top_level_distances[i].resize(count_from_front+count_from_middle+count_from_back);
                top_level_start_end_positions[i].resize(count_from_front+count_from_middle+count_from_back);
                //top_level_indices[i].second.resize(count_from_front+count_from_middle+count_from_back);
                //top_level_distances[i].resize(count_from_front+count_from_middle+count_from_back);
            }
      
            for(int i=0;i<top_paths_num;++i)
            {

                //for(int j=0;j<top_dist[i].size();++j)
                //{
                //std::cout << "Pair " << j << ": (" << top_dist[i][j].first << ", " << top_dist[i][j].second << ")" << std::endl;
                //}
                
                //std::cout<<"执行到here"<<std::endl;
                //std::cout<<"top_level_points[i]"<<top_level_points[i]<<std::endl;
                //top_level_indices[i].first=0;
                //std::cout<<"执行成功"<<std::endl;
                top_level_indices[i].first=top_level_points[i];
                //std::cout<<"执行到here"<<std::endl;

                std::sort(top_dist[i].begin(),top_dist[i].end());

                //std::nth_element(top_dist[i].begin(),top_dist[i].begin()+count_from_front,top_dist[i].end());
                //std::cout<<"执行到here"<<std::endl;
                //std::sort(top_dist[i].begin(),top_dist[i].begin()+count_from_front);  
                //std::cout<<"执行到here"<<std::endl;
                for(size_t j=0;j<count_from_front;++j)
                {   //std::cout<<"加入"<<top_dist[i][j].second<<std::endl;
                    //分类
                    top_level_distances[i][j]=top_dist[i][j].first;
                    top_level_indices[i].second[j]=top_dist[i][j].second.first;
                    top_level_start_end_positions[i][j].first=top_dist[i][j].second.second.first;
                    top_level_start_end_positions[i][j].second=top_dist[i][j].second.second.second;
                    //top_level_distances[i][j]=top_dist[i][j].first;
                    //top_level_indices[i].second[j]=top_dist[i][j].second;
                }
                //std::cout<<"前一部分"<<std::endl;


                int count_from_middle_start=count_from_front;
                int count_from_middle_end=top_paths_num-count_from_back;
                //std::nth_element(top_dist[i].begin()+count_from_middle_start,top_dist[i].begin()+count_from_middle_end-1,top_dist[i].end());
                for(size_t j=0;j<count_from_middle;j++)
                {
                    //std::cout<<"加入"<<top_dist[i][count_from_front+j].second<<std::endl;
                    //分类
                    top_level_distances[i][count_from_front+j]=top_dist[i][count_from_front+j].first;
                    top_level_indices[i].second[count_from_front+j]=top_dist[i][count_from_front+j].second.first;
                    top_level_start_end_positions[i][count_from_front+j].first=top_dist[i][count_from_front+j].second.second.first;
                    top_level_start_end_positions[i][count_from_front+j].second=top_dist[i][count_from_front+j].second.second.second;
                    //top_level_distances[i][count_from_front+j]=top_dist[i][count_from_front+j].first;
                    //top_level_indices[i].second[count_from_front+j]=top_dist[i][count_from_front+j].second;
                }
                //std::cout<<"中间部分"<<std::endl;

                //std::sort(top_dist[i].begin()+count_from_middle_end-1,top_dist[i].end());
                for(size_t j=0;j<count_from_back;j++)
                {
                    //std::cout<<"加入"<<top_dist[i][count_from_middle_end+j-1].second<<std::endl;
                    //分类
                    top_level_distances[i][count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].first;
                    top_level_indices[i].second[count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].second.first;
                    top_level_start_end_positions[i][count_from_front+count_from_middle+j].first=top_dist[i][count_from_middle_end+j-1].second.second.first;
                    top_level_start_end_positions[i][count_from_front+count_from_middle+j].second=top_dist[i][count_from_middle_end+j-1].second.second.second;
                    //top_level_distances[i][count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].first;
                    //top_level_indices[i].second[count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].second;
                }
                //std::cout<<"后面部分"<<std::endl;

            }
            //for(int j=0;j<top_paths_num;j++)
            //{
            //    std::cout<<top_level_points[j]<<" ";
            //}
            //std::cout<<std::endl;
            std::cout<<"top_level_graph_distances_and_indices构建完成"<<std::endl;
        }


/////////////双重过滤那一版本可以，在原来r-tree过滤的基础上增加过滤达到双重过滤
//         //构建graph*********************************************************************************
//         void naive_construction(const string & algorithm, const std::vector<path>& paths, const size_t k)
//         {
 
//             int paths_num=paths.size();

//             // 自动计算网格大小，将整个数据集平均分为1000个网格
//             //int num_grids = 10;
//             //auto total_grid = calculate_grid_size(paths, num_grids);

//             // 创建R-tree索引
//             // 插入轨迹数据到R-tree中

//             auto start_build_rtree=std::chrono::high_resolution_clock::now();
//             rtree_t rtree;

// ////            #pragma omp parallel
// ////            {     
// ////            int num_threads=omp_get_num_threads();
// ////            int thread_id=omp_get_thread_num();
// ////            int chunk_size=(paths_num+num_threads-1)/num_threads;
// ////            int start = thread_id * chunk_size;
// ////            int end = std::min((thread_id + 1) * chunk_size, paths_num);

//             ////#pragma omp parallel for
//             for (unsigned i = 0; i < paths.size(); ++i) 
//             {   
//                 const auto& current_traj=paths[i];
//                 rtree.insert(std::make_pair(to_traj_grid_range(current_traj), i));
//             }
//             auto end_build_rtree=std::chrono::high_resolution_clock::now();
//             std::chrono::duration<double> build_rtree=end_build_rtree-start_build_rtree;
//             std::cout<<"build rtree time: "<<build_rtree.count()<<" seconds\n"<<std::endl;
//             std::cout<<"rtree_build执行完成"<<std::endl;




//             std::cout<<"开始执行grid_top_level_points"<<std::endl;
//             grid_top_level_points_build(algorithm,paths,rtree);
//             std::cout<<"执行grid_top_level_points完成"<<std::endl;




//             //indices=load_from_file<int>("/home/KNN/orders_data/trajectory_new/saved_graph_index3000_front_back_test");
//             //distances=load_from_file<double>("/home/KNN/orders_data/trajectory_new/saved_graph_dist3000_front_back_test");

//             std::cout<<"开始构建底层图"<<std::endl;

// ////            #pragma omp barrier

//             auto start_build_graph_node=std::chrono::high_resolution_clock::now();
            
//             // 构建kNN图
//             //得到前0.8k个元素
//             int portion_front=0.8*k;
//             int portion_back=0.2*k;
//             //#pragma omp parallel for schedule(dynamic)
//             #pragma omp parallel for 
//             for(size_t i=0;i<paths_num;++i)
//             {
//                 //std::cout<<"第 "<<i<<" 条轨迹"<<std::endl;
                
//                 const auto& current_traj=paths[i];

//                 std::vector<value> result_grid;
//                 box_t current_traj_box=to_traj_grid_range(current_traj);

                
//                 //底层图最近邻过滤的方式********************************************************
//                 //1.intersects
//                 // rtree.query(bgi::intersects(current_traj_box),std::back_inserter(result_grid));
//                 // double Mul_grid=1.1;
//                 // while(result_grid.size()<270)
//                 // {
//                 //     current_traj_box=to_traj_Mul_grid_range(current_traj,Mul_grid);
//                 //     result_grid.clear();
//                 //     rtree.query(bgi::intersects(current_traj_box),std::back_inserter(result_grid));
//                 //     Mul_grid+=0.1;
//                 // }
//                 // std::cout<<"执行到这里"<<std::endl;
//                 //2.within
//                 // rtree.query(bgi::within(current_traj_box),std::back_inserter(result_grid));
//                 // double Mul_grid=1.1;
//                 // while(result_grid.size()<800)
//                 // {   
//                 //     //std::cout<<"搜到了"<<Mul_grid<<"倍grid"<<std::endl;
//                 //     current_traj_box=to_traj_Mul_grid_range(current_traj,Mul_grid);
//                 //     result_grid.clear();
//                 //     rtree.query(bgi::within(current_traj_box),std::back_inserter(result_grid));
//                 //     Mul_grid+=0.1;
//                 // }
//                 //std::cout<<"执行到1"<<std::endl;
//                 //3.nearest
//                 rtree.query(bgi::nearest(current_traj_box,900),std::back_inserter(result_grid));
//                 //应该用indices、distances代替;
//                 //分类
//                 std::vector<std::pair<double,std::pair<size_t,std::pair<int,int>>>> front_knn_distances;
//                 std::vector<std::pair<double,std::pair<size_t,std::pair<int,int>>>> back_knn_distances;
//                 std::set<int> selected_ids; //nearest范围之内的索引
//                 std::vector<int> remaining_ids; //nearest范围之外的索引
//                 std::pair<std::set<int>, double> filter_result_grid; //nearest范围之内经过过滤的索引和对应的最小距离
//                 std::set<int> result_ids; //nearest范围之内经过过滤的索引

//                 //过滤
//                 ///////////////////////////////////////////////////////////////////////////////////
//                 selected_ids=run_get_index(paths,result_grid,0.001);//得到nearest索引，生成反向索引
//                 //std::cout<<"result_grid.size()"<<result_grid.size()<<" selected_ids.size()"<<selected_ids.size()<<std::endl;
//                 //std::cout<<"执行到2"<<std::endl;
//                 filter_result_grid=osf_filter(current_traj,0.1,0.001,paths,result_grid,"xian",matricsType);
//                 //std::cout<<"执行到3"<<std::endl;
//                 result_ids=filter_result_grid.first;
//                 //因为result_grid不仅是索引，所以处理selected_ids、result_ids
//                 selected_ids.erase(i);
//                 result_ids.erase(i);
//                 //std::cout<<"执行到4"<<std::endl;
//                 std::cout<<"过滤后的索引数selected_ids "<<selected_ids.size()<<std::endl;
//                 std::cout<<"过滤后的索引数result_ids "<<result_ids.size()<<std::endl;
//                 if(result_ids.size()<portion_front)
//                 {
//                     int needed=portion_front-result_ids.size();
                    
//                     //过滤selected_ids中不在result_ids的元素
//                     std::vector<int> additional_ids;

//                     for(int id:selected_ids)
//                     {
//                         if(result_ids.find(id)==result_ids.end())
//                         {
//                             additional_ids.push_back(id);
//                         }
//                     }

//                     std::random_device rd;
//                     std::mt19937 gen(rd());

//                     if(additional_ids.size()<needed)
//                     {
//                         std::cout<<"Warning: Not enough additional IDs available."<<std::endl;
//                     }
//                     else 
//                     {
//                         std::shuffle(additional_ids.begin(),additional_ids.end(),gen);
//                         additional_ids.resize(needed);

//                         result_ids.insert(additional_ids.begin(),additional_ids.end());
//                     }
//                 }
//                 std::cout<<"过滤之后计算"<<result_ids.size()<<"条轨迹"<<std::endl;
//                 ///////////////////////////////////////////////////////////////////////////////////


//                 // if(i%1000==0)
//                 // {
//                 //     std::cout<<"front:  这条轨迹box范围涉及到并计算的轨迹数为 "<<result_grid.size()<<endl;
//                 // }
//                 //std::cout<<"这条轨迹第一个grid范围涉及到并计算的轨迹数为 "<<result_grid.size()<<endl;
                        
//                 selected_ids.insert(i);
                
//                 //#pragma omp parallel for schedule(dynamic)


//                 // std::cout<<"当前轨迹与其他轨迹之间的距离"<<std::endl;
//                 // int not_zero=0;
//                 // for(int tt=1;tt<paths_num;tt++)
//                 // {
//                 //     double path_distance=minSubTrajectory2(algorithm,current_traj,paths[tt]);
//                 //     if(path_distance<0)
//                 //     {
//                 //         not_zero+=1;
//                 //     }
//                 //     std::cout<<"第 1 条轨迹与第 "<<tt+1<<" 条轨迹之间的距离为: "<<path_distance<<std::endl;
//                 // }
//                 // std::cout<<"与当前轨迹相似的轨迹个数为"<<not_zero<<std::endl;

//                 int score_num=0;
//                 #pragma omp parallel for
//                 for(const auto& res:result_ids)
//                 {
//                     unsigned j=res;
//                     //selected_ids.insert(j);
//                     if(i!=j)
//                     {
//                         //double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
//                         //分类
//                         //std::pair<double, std::pair< std::pair<int,int>, std::pair<int,int> > > result; minSubTrajectory2的返回类型;
//                         //(path_distance,((start1,end1),(start2,end2)))=minSubTrajectory2(algorithm,current_traj,paths[j]);
//                         //std::cout<<"底层图front:"<<"轨迹"<<i<<"与其他轨迹"<<j<<"计算相似得分, ";
//                         Result result=minSubTrajectory2(algorithm,current_traj,paths[j]); 
//                         // if(result.first<0)
//                         // {
//                         //     score_num+=1;
//                         // }
//                         //std::cout<<"第"<<i<<"个轨迹从 "<<(result.second.first.first)<<" 到 "<<(result.second.first.second)<<"第"<<j<<"个轨迹从 "<<(result.second.second.first)<<" 到 "<<(result.second.second.second)<<" 代表性相似子轨迹得分是 "<<(-result.first)<<endl;
//                         //#pragma omp critical
//                         front_knn_distances.push_back(std::make_pair(result.first,std::make_pair(j,std::make_pair(result.second.first.first,result.second.first.second))));
//                     }
//                     else
//                     {
//                         std::cout<<"索引i是 "<<i<<" 索引j是 "<<j<<std::endl;
//                         exit(1);
//                     }
//                 }
//                 // std::cout<<"底层图front:"<<"轨迹"<<i<<"与其他轨迹相似得分大于零的数量为 "<<score_num<<std::endl;
//                 // std::cout<<std::endl<<std::endl;
//                 //std::cout<<"执行到这里"<<std::endl;

//                 std::sort(front_knn_distances.begin(),front_knn_distances.end());

//             //*********************************************************************************************************************/
//             //前面放一个集合,去掉i,这样到此直接就是remaining_ids,不用再使用find()对每个都去处理
//                 for(int t=0;t<paths_num;t++)
//                 {
//                     if(selected_ids.find(t)==selected_ids.end())
//                     {
//                         remaining_ids.push_back(t);
//                     }
//                 }
//                 //std::cout<<"执行到这里"<<std::endl;
        
        
        
//                 // if(front_knn_distances.size()>=150)
//                 // {
//                 //     for(size_t ttt=150;ttt<front_knn_distances.size();ttt++)
//                 //     {
//                 //         remaining_ids.push_back(front_knn_distances[ttt].second);
//                 //     }
//                 // }
                
//                 // std::cout<<"执行到这里"<<std::endl;
                
//                 std::random_device rd;
//                 std::mt19937 g(rd());
//                 std::shuffle(remaining_ids.begin(),remaining_ids.end(),g);
//                 std::vector<int> random_indices(remaining_ids.begin(),remaining_ids.begin()+portion_back);
//                 //std::cout<<"执行到这里"<<std::endl;
//                 //std::cout<<"执行到这"<<std::endl;
               
//                 //#pragma omp parallel for schedule(dynamic)
//                 #pragma omp parallel for
//                 for(int t=0;t<random_indices.size();t++)
//                 {
//                     unsigned j=random_indices[t];
//                     //std::cout<<"底层图back:"<<"轨迹"<<i<<"与其他轨迹"<<j<<"计算相似得分, ";
//                     Result result=minSubTrajectory2(algorithm,current_traj,paths[j]); 
//                     back_knn_distances.push_back(std::make_pair(result.first,std::make_pair(j,std::make_pair(result.second.first.first,result.second.first.second))));

//                     ////double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
//                     //#pragma omp critical
//                     ////back_knn_distances.push_back(std::make_pair(path_distance,j));
//                 }
//                 //std::cout<<"执行到这里"<<std::endl;


//                 // if(remaining_ids.size()<50)
//                 // {
//                 //     size_t num_to_get=100;
//                 //     size_t start_index = std::max(0ul, front_knn_distances.size() - num_to_get);
//                 //     for(size_t ttt=start_index;ttt<front_knn_distances.size();++ttt)
//                 //     {
//                 //         back_knn_distances.push_back(front_knn_distances[ttt]);
//                 //     }
//                 // }
//                 // else 
//                 // {
//                 //     std::random_device rd;
//                 //     std::mt19937 g(rd());
//                 //     std::shuffle(remaining_ids.begin(),remaining_ids.end(),g);
//                 //     std::vector<int> random_indices(remaining_ids.begin(),remaining_ids.begin()+portion_back);

//                 //     //std::cout<<"执行到这"<<std::endl;
                
//                 //     //#pragma omp parallel for schedule(dynamic)
//                 //     #pragma omp parallel for
//                 //     for(int t=0;t<random_indices.size();t++)
//                 //     {
//                 //         unsigned j=random_indices[t];
//                 //         double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
//                 //         //#pragma omp critical
//                 //         back_knn_distances.push_back(std::make_pair(path_distance,j));
//                 //     }
//                 //     //std::cout<<"执行到这"<<std::endl;


//                 // }
               
                

                
//                 std::sort(back_knn_distances.begin(),back_knn_distances.begin());

//                 //std::cout<<"执行到这"<<std::endl;




//                 //random nearests
//                 // std::random_device rd_;
//                 // std::mt19937 g_(rd_());
//                 // std::shuffle(g_);
//                 // std::vector<int> random_nearest();

//             //    int res_num=front_knn_distances.size();
//             //    std::cout<<"within "<<res_num<<" 条轨迹与这个data轨迹之间的距离: "<<std::endl;
//             //    for(size_t j=0;j<res_num;++j)
//             //    {
//             //        std::cout<<"这个data轨迹与within "<<j+1<<" 索引号为 "<<front_knn_distances[j].second<<" 轨迹之间的距离 "<<front_knn_distances[j].first<<std::endl;
//             //    }
//             //    std::cout<<std::endl;
//             //    int res_num1=0;
//                 //std::cout<<"第"<<i<<"条轨迹的存储内容"<<std::endl;
//                 for(size_t j=0;j<portion_front;++j)
//                 {
//                     //分类
//                     distances[i][j]=front_knn_distances[j].first;
//                     indices[i][j]=front_knn_distances[j].second.first;
//                     start_end_positions[i][j].first=front_knn_distances[j].second.second.first;
//                     start_end_positions[i][j].second=front_knn_distances[j].second.second.second;
//                     //std::cout<<"第"<<j<<"个邻居与我的start "<<start_end_positions[i][j].first<<" to end "<<start_end_positions[i][j].second<<"相似"<<std::endl;

//                     //distances[i][j]=front_knn_distances[j].first;
//                     //indices[i][j]=front_knn_distances[j].second;
//             //       std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这个data轨迹之间的距离: "<<distances[i][j]<<std::endl;
//             //       res_num1++;
//                 }
//                 //std::cout<<std::endl;
//             //     std::cout<<"执行到这里"<<std::endl;
//             //     for(size_t j=0;j<20;++j)
//             //     {
//             //         distances[i][60+j]=front_knn_distances[60+4*j].first;
//             //         indices[i][60+j]=front_knn_distances[60+4*j].second;
//             // //        std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这个data轨迹之间的距离: "<<distances[i][60+j]<<std::endl;
//             // //       res_num1++;
//             //     }
//             //     std::cout<<"执行到这里"<<std::endl;




//                 // for(size_t j=0;j<portion_front;++j)
//                 // {
//                 //     distances[i][j]=front_knn_distances[j].first;
//                 //     indices[i][j]=front_knn_distances[j].second;
//                 // }
//                 //std::cout<<"执行到这"<<std::endl;
//                 for(size_t j=0;j<portion_back;++j)
//                 {
//                     distances[i][portion_front+j]=back_knn_distances[j].first;
//                     indices[i][portion_front+j]=back_knn_distances[j].second.first;
//                     start_end_positions[i][portion_front+j].first=back_knn_distances[j].second.second.first;
//                     start_end_positions[i][portion_front+j].second=back_knn_distances[j].second.second.second;
//                     //distances[i][portion_front+j]=back_knn_distances[j].first;
//                     //indices[i][portion_front+j]=back_knn_distances[j].second;
//             //       std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这条data轨迹之间的距离为: "<<distances[i][portion_front+j]<<std::endl;
//             //       res_num1++;
//                 }
//                 ////std::cout<<"执行到这里"<<std::endl;
//                //exit(1);
//                 //std::cout<<"graph_node_build执行完成"<<std::endl;
//                 //std::cout<<"graph_node_build执行完成"<<std::endl;
// ////            }                
//             }

//             //indices=load_from_file<int>("/home/KNN/orders_data/trajectory_new_chengdu/saved_graph_index5000_front_back_test_final_final");
//             //distances=load_from_file<double>("/home/KNN/orders_data/trajectory_new_chengdu/saved_graph_dist5000_front_back_test_final_final");
            

//             //3.过滤+类似最开始的顶层图
//             // std::cout<<"开始创建top-level-graph"<<std::endl;
//             // auto start_build_top_level_graph=std::chrono::high_resolution_clock::now();
//             // //front_knn_distances的大小一定会大于2*k/3*k或者是自己设置的nearest的值，可以根据这个和想要消掉的邻居数设置drop_neighbors_num;
//             // int drop_neighbors=0.9*k;
//             // std::vector<char> top_level_points_visitpool(paths_num,false);
//             // std::unordered_set<int> top_level_points_set;
//             // for(int t=0;t<paths_num;t++)
//             // {
//             //     top_level_points_set.insert(t);
//             // }

//             // int add_top_level_point=rand()%paths_num;
//             // top_level_points.push_back(add_top_level_point);

//             // top_level_points_visitpool[add_top_level_point]=true;
//             // top_level_points_set.erase(add_top_level_point);
            
//             // while(!top_level_points_set.empty())
//             // {
//             //     //std::cout<<"add_top_level_point "<<add_top_level_point<<std::endl;
//             //     for(int j=0;j<drop_neighbors;++j)
//             //     {
//             //         top_level_points_visitpool[indices[add_top_level_point][j]]=true;
//             //         top_level_points_set.erase(indices[add_top_level_point][j]);
//             //     }
//             //     auto it=std::next(top_level_points_set.begin(),std::rand()%top_level_points_set.size());
//             //     add_top_level_point=*it;
//             //     top_level_points.push_back(add_top_level_point);
                    
//             //     top_level_points_visitpool[add_top_level_point]=true;
//             //     top_level_points_set.erase(add_top_level_point);
//             //     //std::cout<<"又选了一个add_top_level_point"<<std::endl;
//             // }


//             // for(int t=0;t<top_level_points.size();t++)
//             // {
//             //     std::cout<<top_level_points[t]<<" ";
//             // }
//             // std::cout<<std::endl;

//             // int top_paths_num=top_level_points.size();
//             // std::cout<<"top_level_points_set的size "<<top_paths_num<<std::endl;
//             // top_level_indices.resize(top_paths_num);
//             // top_level_distances.resize(top_paths_num);
//             // for(int t=0;t<top_paths_num;t++)
//             // {
//             //     top_level_indices[t].first=top_level_points[t];
//             //     top_level_indices[t].second.resize(top_paths_num-1);
//             //     top_level_distances[t].resize(top_paths_num-1);
//             // }

//             // std::sort(top_level_points.begin(),top_level_points.end());

//             // for(int t=0;t<top_paths_num;t++)
//             // {                
//             //     for(int tt=t+1;tt<top_paths_num;tt++)
//             //     {
//             //         //std::cout<<"执行到t="<<t<<" tt="<<tt<<std::endl;
//             //         double top_paths_distance=minSubTrajectory2(algorithm,paths[top_level_points[t]],paths[top_level_points[tt]]);
//             //         top_level_distances[t][tt-1]=top_paths_distance;
//             //         top_level_indices[t].second[tt-1]=top_level_points[tt];
//             //         //std::cout<<"执行到此"<<std::endl;

//             //         top_level_distances[tt][t]=top_paths_distance;
//             //         //std::cout<<"执行到此"<<std::endl;
//             //         top_level_indices[tt].second[t]=top_level_points[t];
//             //         //std::cout<<"执行到此"<<std::endl;
//             //     }

//             // }
//             // auto end_build_top_level_graph=std::chrono::high_resolution_clock::now();
//             // std::chrono::duration<double> build_top_level_graph=end_build_top_level_graph-start_build_top_level_graph;
//             // std::cout<<"build top_level_graph time: "<<build_top_level_graph.count()<<" seconds\n"<<std::endl;
//             // std::cout<<"top_level_graph_build构建完成"<<std::endl;
//             auto end_build_graph_node=std::chrono::high_resolution_clock::now();
//             std::chrono::duration<double> build_graph_node=end_build_graph_node-start_build_graph_node;
//             std::cout<<"build graph_node time: "<<build_graph_node.count()<<" seconds\n"<<std::endl;
//             std::cout<<"graph_node_build执行完成"<<std::endl;
//         }
//第三版并发，没有运行
// void naive_construction(const std::string &algorithm, const std::vector<path> &paths, const size_t k) {
//     std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
//     int paths_num = paths.size();
//     auto start_build_rtree = std::chrono::high_resolution_clock::now();
//     rtree_t rtree;

//     // Build R-tree in parallel
//     #pragma omp parallel
//     {
//         std::vector<std::pair<box_t, int>> local_inserts;
        
//         #pragma omp for
//         for (unsigned i = 0; i < paths_num; ++i) {
//             const auto &current_traj = paths[i];
//             local_inserts.emplace_back(to_traj_grid_range(current_traj), i);
//         }

//         // Only one thread should insert into the R-tree
//         #pragma omp critical
//         {
//             for (const auto &insert : local_inserts) {
//                 rtree.insert(insert);
//             }
//         }
//     }
    
//     auto end_build_rtree = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> build_rtree = end_build_rtree - start_build_rtree;
//     std::cout << "build rtree time: " << build_rtree.count() << " seconds\n" << std::endl;

//     std::cout << "开始执行grid_top_level_points" << std::endl;
//     grid_top_level_points_build(algorithm, paths, rtree);
//     std::cout << "执行grid_top_level_points完成" << std::endl;

//     std::cout << "开始构建底层图" << std::endl;
//     auto start_build_graph_node = std::chrono::high_resolution_clock::now();
//     int portion_front = static_cast<int>(0.8 * k);
//     int portion_back = static_cast<int>(0.2 * k);

//     // Prepare data structures
//     std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> front_knn_distances(paths_num);
//     std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> back_knn_distances(paths_num);

//     // Build reverse index once
//     std::vector<std::set<int>> global_selected_ids = run_get_index(paths, {}, 0.001);  // 假设这里的 {} 代表所有路径
//     // Main parallel loop
//     #pragma omp parallel
//     {
//         std::vector<std::set<int>> selected_ids(paths_num);

//         #pragma omp for
//         for (size_t i = 0; i < paths_num; ++i) {
//             std::cout << "1" << std::endl;
//             const auto &current_traj = paths[i];
//             std::vector<value> result_grid;
//             box_t current_traj_box = to_traj_grid_range(current_traj);

//             rtree.query(bgi::nearest(current_traj_box, 900), std::back_inserter(result_grid));

//             auto filter_result_grid = osf_filter(current_traj, 0.1, 0.001, paths, result_grid, "xian", matricsType);
//             auto result_ids = filter_result_grid.first;

//             selected_ids[i] = global_selected_ids[i]; // 使用预先构建的反向索引
//             selected_ids[i].erase(i);
//             result_ids.erase(i);

//             std::cout << "2" << std::endl;
//             // Ensure sufficient IDs for front KNN
//             if (result_ids.size() < portion_front) {
//                 int needed = portion_front - result_ids.size();
//                 std::vector<int> additional_ids;

//                 for (int id : selected_ids[i]) {
//                     if (result_ids.find(id) == result_ids.end()) {
//                         additional_ids.push_back(id);
//                     }
//                 }

//                 std::random_device rd;
//                 std::mt19937 gen(rd());
//                 std::shuffle(additional_ids.begin(), additional_ids.end(), gen);

//                 if (additional_ids.size() >= needed) {
//                     additional_ids.resize(needed);
//                     result_ids.insert(additional_ids.begin(), additional_ids.end());
//                 } else {
//                     std::cout << "Warning: Not enough additional IDs available." << std::endl;
//                 }
//             }

//             std::cout << "3" << std::endl;
//             selected_ids[i].insert(i);
//             // Calculate front KNN distances
//             std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>> local_front_knn;
//             for (const auto &res : result_ids) {
//                 unsigned j = res;
//                 if (i != j) {
//                     Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
//                     local_front_knn.emplace_back(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second)));
//                 }
//             }

//             std::cout << "4" << std::endl;
//             #pragma omp critical
//             {
//                 front_knn_distances[i].insert(front_knn_distances[i].end(), local_front_knn.begin(), local_front_knn.end());
//             }

//             // Back KNN calculations
//             std::vector<int> remaining_ids;
//             for (int t = 0; t < paths_num; t++) {
//                 if (selected_ids[i].find(t) == selected_ids[i].end()) {
//                     remaining_ids.push_back(t);
//                 }
//             }

//             std::random_device rd;
//             std::mt19937 g(rd());
//             std::shuffle(remaining_ids.begin(), remaining_ids.end(), g);
//             std::vector<int> random_indices(remaining_ids.begin(), remaining_ids.begin() + std::min(portion_back, static_cast<int>(remaining_ids.size())));

//             std::cout << "5" << std::endl;
//             std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>> local_back_knn;
//             for (int t : random_indices) {
//                 unsigned j = t;
//                 Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
//                 local_back_knn.emplace_back(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second)));
//             }

//             std::cout << "6" << std::endl;
//             #pragma omp critical
//             {
//                 back_knn_distances[i].insert(back_knn_distances[i].end(), local_back_knn.begin(), local_back_knn.end());
//             }
//         }
//     }

//     // Post-processing results
//     for (size_t i = 0; i < paths_num; ++i) {
//         std::sort(front_knn_distances[i].begin(), front_knn_distances[i].end());
//         for (size_t j = 0; j < portion_front; ++j) {
//             distances[i][j] = front_knn_distances[i][j].first;
//             indices[i][j] = front_knn_distances[i][j].second.first;
//             start_end_positions[i][j].first = front_knn_distances[i][j].second.second.first;
//             start_end_positions[i][j].second = front_knn_distances[i][j].second.second.second;
//         }

//         std::sort(back_knn_distances[i].begin(), back_knn_distances[i].end());
//         for (size_t j = 0; j < portion_back; ++j) {
//             distances[i][portion_front + j] = back_knn_distances[i][j].first;
//             indices[i][portion_front + j] = back_knn_distances[i][j].second.first;
//             start_end_positions[i][portion_front + j].first = back_knn_distances[i][j].second.second.first;
//             start_end_positions[i][portion_front + j].second = back_knn_distances[i][j].second.second.second;
//         }
//     }

//     auto end_build_graph_node = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> build_graph_node = end_build_graph_node - start_build_graph_node;
//     std::cout << "build graph_node time: " << build_graph_node.count() << " seconds\n" << std::endl;
//     std::cout << "graph_node_build执行完成" << std::endl;
// }
//第二版并发有时会出错，可以运行
void naive_construction(const std::string &algorithm, const std::vector<path> &paths, const size_t k) {
    omp_set_num_threads(50);
    std::cout<<"Number of threads: "<<omp_get_max_threads()<<std::endl;
    int paths_num = paths.size();
    auto start_build_rtree = std::chrono::high_resolution_clock::now();
    rtree_t rtree;

    //构建rtree的时候，利用local_inserts来存储所有轨迹的矩形，然后插入到rtree中，使用local_inserts实现并行
    // Build R-tree in parallel
    #pragma omp parallel
    {
        std::vector<std::pair<box_t, int>> local_inserts;
        
        #pragma omp for
        for (unsigned i = 0; i < paths_num; ++i) {
            const auto &current_traj = paths[i];
            local_inserts.emplace_back(to_traj_grid_range(current_traj), i);
        }

        // Only one thread should insert into the R-tree
        #pragma omp critical
        {
            for (const auto &insert : local_inserts) {
                rtree.insert(insert);
            }
        }
    }
    
    auto end_build_rtree = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_rtree = end_build_rtree - start_build_rtree;
    std::cout << "build rtree time: " << build_rtree.count() << " seconds\n" << std::endl;

    std::cout << "开始执行grid_top_level_points" << std::endl;
    grid_top_level_points_build(algorithm, paths, rtree);
    std::cout << "执行grid_top_level_points完成" << std::endl;

    std::cout << "开始构建底层图" << std::endl;
    auto start_build_graph_node = std::chrono::high_resolution_clock::now();
    int portion_front = static_cast<int>(1.0 * k);
    int portion_back = static_cast<int>(0* k);

    //所有轨迹的邻居都存储在front_knn_distances,back_knn_distances中，实现并行化
    // Prepare data structures
    std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> front_knn_distances(paths_num);
    std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> back_knn_distances(paths_num);

    //run_get_index(paths,0.001);
    // Main parallel loop
    //std::vector<std::set<int>> selected_ids(paths_num);
    #pragma omp parallel
    {
        //std::vector<value> result_grid;
        //std::set<int> selected_ids;
        //std::map<std::pair<int, int>, std::set<int>> invert_index;

        //做了这一处修改,上面的selected_ids、invert_index内移到for循环内可以加快，但是local_front_knn,local_back_knn内移到for循环内时间不变480s没有加快性能从0.5到0.3，去除只留front_knn_distances，back_knn_distances时间还会有所增加600s性能到0.18；
        std::vector<std::vector<std::pair<double, std::pair<int,std::pair<double,double>>>>> local_front_knn(paths_num);
        std::vector<std::vector<std::pair<double, std::pair<int,std::pair<double,double>>>>> local_back_knn(paths_num);

        #pragma omp for
        for (size_t i = 0; i < paths_num; ++i) {
            //std::cout<<"1"<<std::endl;
            const auto &current_traj = paths[i];
            //std::vector<value> result_grid;

            //做了这一处修改
            std::vector<value> result_grid;
            std::set<int> selected_ids;
            std::map<std::pair<int, int>, std::set<int>> invert_index;


            box_t current_traj_box = to_traj_grid_range(current_traj);
            //std::cout<<"1**"<<std::endl;
            //*******可以将result_grid定义为和selected_ids一样的path_num个大小
            rtree.query(bgi::nearest(current_traj_box, 7500), std::back_inserter(result_grid));
            //std::cout<<"2//"<<std::endl;
            // for(const auto& result:result_grid)
            // {
            //     selected_id[i].insert(result.second);
            // }
            //********selected_ids[i]就是result_grid.second,得到的就是第一层过滤后的索引
            //std::cout<<"3&&"<<std::endl;
            std::tie(selected_ids, invert_index)=run_get_index(paths,result_grid,0.001);
            //std::cout<<"3&&"<<std::endl;
            //********nearest过滤得到的索引放到run_get_index中去构建反向索引，这个反向索引也就是第一层过滤后的轨迹作为osf过滤的基本，osf从这里面过滤会留下更少的轨迹，加快计算
            auto filter_result_grid = osf_filter(current_traj, invert_index, 0.1, 0.001, paths, result_grid, "xian", matricsType);
            auto result_ids = filter_result_grid.first;
            //std::cout<<"4##"<<std::endl;
            selected_ids.erase(i);
            result_ids.erase(i);
            //std::cout<<"5^^"<<std::endl;
            //std::cout<<"2"<<std::endl;
            // Ensure sufficient IDs for front KNN
            if (result_ids.size() < portion_front) {
                int needed = portion_front - result_ids.size();
                std::vector<int> additional_ids;

                for (int id : selected_ids) {
                    if (result_ids.find(id) == result_ids.end()) {
                        additional_ids.push_back(id);
                    }
                }

                std::random_device rd;
                std::mt19937 gen(rd());
                std::shuffle(additional_ids.begin(), additional_ids.end(), gen);

                if (additional_ids.size() >= needed) {
                    std::shuffle(additional_ids.begin(),additional_ids.end(),gen);
                    additional_ids.resize(needed);
                    result_ids.insert(additional_ids.begin(), additional_ids.end());
                } else {
                    std::cout << "Warning: Not enough additional IDs available." << std::endl;
                }
            }

            //std::cout<<"3"<<std::endl;
            selected_ids.insert(i);
            // Calculate front KNN distances
            //std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>> local_front_knn;
            for (const auto &res : result_ids) {
                unsigned j = res;
                if (i != j) {
                    Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
                    local_front_knn[i].emplace_back(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second)));
                }
            }

            //std::cout<<"4"<<std::endl;
            //#pragma omp critical
            //{
            //    front_knn_distances[i].insert(front_knn_distances[i].end(), local_front_knn.begin(), local_front_knn.end());
            //}

            // Back KNN calculations
            std::vector<int> remaining_ids;
            for (int t = 0; t < paths_num; t++) {
                if (selected_ids.find(t) == selected_ids.end()) {
                    remaining_ids.push_back(t);
                }
            }

            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(remaining_ids.begin(), remaining_ids.end(), g);
            std::vector<int> random_indices(remaining_ids.begin(), remaining_ids.begin() + std::min(portion_back, static_cast<int>(remaining_ids.size())));

            //std::cout<<"5"<<std::endl;
            //std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>> local_back_knn;
            for (int t : random_indices) {
                unsigned j = t;
                Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
                local_back_knn[i].emplace_back(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second)));
            }

            //std::cout<<"6"<<std::endl;
            //#pragma omp critical
            //{
            //    back_knn_distances[i].insert(back_knn_distances[i].end(), local_back_knn.begin(), local_back_knn.end());
            //}

            //std::cout<<"4"<<std::endl;
            #pragma omp critical
            {
                front_knn_distances[i].insert(front_knn_distances[i].end(), local_front_knn[i].begin(), local_front_knn[i].end());
                back_knn_distances[i].insert(back_knn_distances[i].end(), local_back_knn[i].begin(), local_back_knn[i].end());
            }
            if(i%1000==0)
            {   
                std::cout<<"执行了一千条轨迹"<<std::endl;
                std::cout<<"这条轨迹计算front时的轨迹条数"<<result_ids.size()<<std::endl;
                std::cout<<"这条轨迹计算back时的轨迹条数"<<random_indices.size()<<std::endl;
                std::cout<<"这条轨迹计算"<<(result_ids.size()+random_indices.size())<<"条轨迹"<<std::endl;
            }
        }
    }

    // Post-processing results
    for (size_t i = 0; i < paths_num; ++i) {
        std::sort(front_knn_distances[i].begin(), front_knn_distances[i].end());
        for (size_t j = 0; j < portion_front; ++j) {
            distances[i][j] = front_knn_distances[i][j].first;
            indices[i][j] = front_knn_distances[i][j].second.first;
            start_end_positions[i][j].first = front_knn_distances[i][j].second.second.first;
            start_end_positions[i][j].second = front_knn_distances[i][j].second.second.second;
        }

        std::sort(back_knn_distances[i].begin(), back_knn_distances[i].end());
        for (size_t j = 0; j < portion_back; ++j) {
            distances[i][portion_front + j] = back_knn_distances[i][j].first;
            indices[i][portion_front + j] = back_knn_distances[i][j].second.first;
            start_end_positions[i][portion_front + j].first = back_knn_distances[i][j].second.second.first;
            start_end_positions[i][portion_front + j].second = back_knn_distances[i][j].second.second.second;
        }
    }

    auto end_build_graph_node = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_graph_node = end_build_graph_node - start_build_graph_node;
    std::cout << "build graph_node time: " << build_graph_node.count() << " seconds\n" << std::endl;
    std::cout << "graph_node_build执行完成" << std::endl;
}
//第一版并发，没有运行是不完整的代码
// void naive_construction(const std::string &algorithm, const std::vector<path> &paths, const size_t k) {
//     int paths_num = paths.size();
//     auto start_build_rtree = std::chrono::high_resolution_clock::now();
//     rtree_t rtree;

//     // Build R-tree in parallel
//     #pragma omp parallel for
//     for (unsigned i = 0; i < paths_num; ++i) {
//         const auto &current_traj = paths[i];
//         rtree.insert(std::make_pair(to_traj_grid_range(current_traj), i));
//     }
    
//     auto end_build_rtree = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> build_rtree = end_build_rtree - start_build_rtree;
//     std::cout << "build rtree time: " << build_rtree.count() << " seconds\n" << std::endl;

//     std::cout << "开始执行grid_top_level_points" << std::endl;
//     grid_top_level_points_build(algorithm, paths, rtree);
//     std::cout << "执行grid_top_level_points完成" << std::endl;

//     std::cout << "开始构建底层图" << std::endl;
//     auto start_build_graph_node = std::chrono::high_resolution_clock::now();
//     int portion_front = static_cast<int>(0.8 * k);
//     int portion_back = static_cast<int>(0.2 * k);

//     // Prepare data structures
//     std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> front_knn_distances(paths_num);
//     std::vector<std::vector<std::pair<double, std::pair<int, std::pair<double, double>>>>> back_knn_distances(paths_num);

//     // Main parallel loop
//     #pragma omp parallel
//     {
//         std::vector<std::unordered_set<int>> selected_ids(paths_num);

//         #pragma omp for
//         for (size_t i = 0; i < paths_num; ++i) {
//             const auto &current_traj = paths[i];
//             std::vector<value> result_grid;
//             box_t current_traj_box = to_traj_grid_range(current_traj);

//             rtree.query(bgi::nearest(current_traj_box,900), std::back_inserter(result_grid));

//             auto filter_result_grid = osf_filter(current_traj, 0.1, 0.001, paths, result_grid, "xian", matricsType);
//             auto result_ids = filter_result_grid.first;

//             selected_ids[i].insert(i);

//             // Ensure sufficient IDs for front KNN
//             if (result_ids.size() < portion_front) {
//                 int needed = portion_front - result_ids.size();
//                 std::vector<int> additional_ids;

//                 for (int id : selected_ids[i]) {
//                     if (result_ids.find(id) == result_ids.end()) {
//                         additional_ids.push_back(id);
//                     }
//                 }

//                 std::random_device rd;
//                 std::mt19937 gen(rd());
//                 std::shuffle(additional_ids.begin(), additional_ids.end(), gen);
                
//                 std::cout<<"additional_ids.size() "<<additional_ids.size()<<std::endl;
//                 std::cout<<"needed "<<needed<<std::endl;

//                 if (additional_ids.size() >= needed) {
//                     additional_ids.resize(needed);
//                     result_ids.insert(additional_ids.begin(), additional_ids.end());
//                 } else {
//                     std::cout << "Warning: Not enough additional IDs available." << std::endl;
//                 }
//             }

//             // Calculate front KNN distances in parallel
//             #pragma omp parallel for
//             for (const auto &res : result_ids) {
//                 unsigned j = res;
//                 if (i != j) {
//                     Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
//                     #pragma omp critical
//                     {
//                         front_knn_distances[i].push_back(std::make_pair(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second))));
//                     }
//                 }
//             }

//             // Back KNN calculations
//             std::vector<int> remaining_ids;
//             for (int t = 0; t < paths_num; t++) {
//                 if (selected_ids[i].find(t) == selected_ids[i].end()) {
//                     remaining_ids.push_back(t);
//                 }
//             }

//             std::random_device rd;
//             std::mt19937 g(rd());
//             std::shuffle(remaining_ids.begin(), remaining_ids.end(), g);
//             std::vector<int> random_indices(remaining_ids.begin(), remaining_ids.begin() + portion_back);

//             #pragma omp parallel for
//             for (int t : random_indices) {
//                 unsigned j = t;
//                 Result result = minSubTrajectory2(algorithm, current_traj, paths[j]);
//                 #pragma omp critical
//                 {
//                     back_knn_distances[i].push_back(std::make_pair(result.first, std::make_pair(j, std::make_pair(result.second.first.first, result.second.first.second))));
//                 }
//             }
//         }
//     }

//     // Post-processing results
//     for (size_t i = 0; i < paths_num; ++i) {
//         std::sort(front_knn_distances[i].begin(), front_knn_distances[i].end());
//         for (size_t j = 0; j < portion_front; ++j) {
//             distances[i][j] = front_knn_distances[i][j].first;
//             indices[i][j] = front_knn_distances[i][j].second.first;
//             start_end_positions[i][j].first = front_knn_distances[i][j].second.second.first;
//             start_end_positions[i][j].second = front_knn_distances[i][j].second.second.second;
//         }

//         std::sort(back_knn_distances[i].begin(), back_knn_distances[i].end());
//         for (size_t j = 0; j < portion_back; ++j) {
//             distances[i][portion_front + j] = back_knn_distances[i][j].first;
//             indices[i][portion_front + j] = back_knn_distances[i][j].second.first;
//             start_end_positions[i][portion_front + j].first = back_knn_distances[i][j].second.second.first;
//             start_end_positions[i][portion_front + j].second = back_knn_distances[i][j].second.second.second;
//         }
//     }

//     auto end_build_graph_node = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> build_graph_node = end_build_graph_node - start_build_graph_node;
//     std::cout << "build graph_node time: " << build_graph_node.count() << " seconds\n" << std::endl;
//     std::cout << "graph_node_build执行完成" << std::endl;
// }
//分类***************************************************************************************************************




//不分类***************************************************************************************************************
//        void grid_top_level_points_build(const string & algorithm, const std::vector<path>& paths,const rtree_t& rtree)
//         {
//             int nums_grid=30;
//             auto result=calculate_grid_size(paths,nums_grid);
//             double min_x,min_y,max_x,max_y;
//             std::pair<double,double> grid_size;

//             std::tie(min_x,min_y,max_x,max_y,grid_size)=result;

//             std::unordered_set<int> top_level_points_set;
            
//             //top层的轨迹索引
//             for(int i=0;i<nums_grid;++i)
//             {
//                 for(int j=0;j<nums_grid;++j)
//                 {
//                     rtree_point min_corner;
//                     rtree_point max_corner;
//                     min_corner=rtree_point(min_x+i*(grid_size.first),min_y+j*(grid_size.second));
//                     max_corner=rtree_point(min_x+(i+1)*(grid_size.first),min_y+(j+1)*(grid_size.second));
//                     bg::model::box<rtree_point> current_grid_box(min_corner,max_corner);
//                     std::vector<value> result_grid;
                    
//                     //1.nearest
//                     //rtree.query(bgi::nearest(current_grid_box,1),std::back_inserter(result_grid));
                    
//                     //2.within
//                     rtree.query(bgi::within(current_grid_box),std::back_inserter(result_grid));
//                     double Mul_grid=1.1;
//                     while(result_grid.size()<1)
//                     {   
//                         double center_x = (bg::get<bg::min_corner,0>(current_grid_box)+bg::get<bg::max_corner,0>(current_grid_box)) / 2.0;
//                         double center_y = (bg::get<bg::min_corner,1>(current_grid_box)+bg::get<bg::max_corner,1>(current_grid_box)) / 2.0;

//                         // 计算新的最小和最大角点，增大两倍
//                         double new_min_x = center_x - (center_x - bg::get<bg::min_corner,0>(current_grid_box)) * Mul_grid;
//                         double new_min_y = center_y - (center_y - bg::get<bg::min_corner,1>(current_grid_box)) * Mul_grid;
//                         double new_max_x = center_x + (bg::get<bg::max_corner,0>(current_grid_box) - center_x) * Mul_grid;
//                         double new_max_y = center_y + (bg::get<bg::max_corner,1>(current_grid_box) - center_y) * Mul_grid;

//                         // 构造并返回新的 box<rtree_point>
//                         rtree_point new_min_corner(new_min_x, new_min_y);
//                         rtree_point new_max_corner(new_max_x, new_max_y);
//                         bg::model::box<rtree_point>current_grid_box(new_min_corner, new_max_corner);

//                         rtree.query(bgi::within(current_grid_box),std::back_inserter(result_grid));
//                         Mul_grid+=0.1;
//                     }
                    
//                     for(const auto& res:result_grid)
//                     {
//                         top_level_points_set.insert(res.second);
//                         //std::cout<<"第"<<i+j+1<<"个网格的轨迹索引为"<<res.second<<std::endl;
//                     }
//                 }
//             }

//             std::copy(top_level_points_set.begin(),top_level_points_set.end(),std::back_inserter(top_level_points));
//             int top_paths_num=top_level_points.size();

//             std::vector<std::vector<std::pair<double,int>>> top_dist(top_paths_num,std::vector<std::pair<double,int>>(top_paths_num-1));
//             for(int i=0;i<top_paths_num;++i)
//             {
//                 int k=i;
//                 for(int j=i+1;j<top_paths_num;++j)
//                 {                     
//                     path v1=paths[top_level_points[i]];
//                     path v2=paths[top_level_points[j]];
//                     //std::cout<<"计算"<<top_level_points[i]<<"和"<<top_level_points[j]<<"之间的相似度"<<std::endl;
//                     double top_path_distance=minSubTrajectory2(algorithm,v1,v2);
//                     top_dist[i][k].first=top_path_distance;
//                     top_dist[i][k].second=top_level_points[j];

//                     top_dist[k+1][i].first=top_path_distance;
//                     top_dist[k+1][i].second=top_level_points[i];
//                     k++;
//                 }
//             }
//             //std::cout<<"top_level_graph相似度计算完成"<<std::endl;
//             //std::cout<<"top_paths_nums为"<<top_paths_num<<std::endl;

//             double portion_from_front=0.1;
//             double portion_from_middle=0.1;
//             double portion_from_back=0.8;

//             //top_level层建立的边的个数
//             int top_level_edge=static_cast<int>(0.6*top_paths_num);

//             int count_from_front=static_cast<int>(top_level_edge*portion_from_front);
//             int count_from_middle=static_cast<int>(top_level_edge*portion_from_middle);
//             int count_from_back=static_cast<int>(top_level_edge*portion_from_back);

//             //分配空间
//             top_level_indices.resize(top_paths_num);
//             top_level_distances.resize(top_paths_num);
//             for(int i=0;i<top_paths_num;++i)
//             {
//                 top_level_indices[i].second.resize(count_from_front+count_from_middle+count_from_back);
//                 top_level_distances[i].resize(count_from_front+count_from_middle+count_from_back);
//             }
      
//             for(int i=0;i<top_paths_num;++i)
//             {

//                 //for(int j=0;j<top_dist[i].size();++j)
//                 //{
//                 //std::cout << "Pair " << j << ": (" << top_dist[i][j].first << ", " << top_dist[i][j].second << ")" << std::endl;
//                 //}
                
//                 //std::cout<<"执行到here"<<std::endl;
//                 //std::cout<<"top_level_points[i]"<<top_level_points[i]<<std::endl;
//                 //top_level_indices[i].first=0;
//                 //std::cout<<"执行成功"<<std::endl;
//                 top_level_indices[i].first=top_level_points[i];
//                 //std::cout<<"执行到here"<<std::endl;

//                 std::sort(top_dist[i].begin(),top_dist[i].end());

//                 //std::nth_element(top_dist[i].begin(),top_dist[i].begin()+count_from_front,top_dist[i].end());
//                 //std::cout<<"执行到here"<<std::endl;
//                 //std::sort(top_dist[i].begin(),top_dist[i].begin()+count_from_front);  
//                 //std::cout<<"执行到here"<<std::endl;
//                 for(size_t j=0;j<count_from_front;++j)
//                 {   //std::cout<<"加入"<<top_dist[i][j].second<<std::endl;
//                     top_level_distances[i][j]=top_dist[i][j].first;
//                     top_level_indices[i].second[j]=top_dist[i][j].second;
//                 }
//                 //std::cout<<"前一部分"<<std::endl;


//                 int count_from_middle_start=count_from_front;
//                 int count_from_middle_end=top_paths_num-count_from_back;
//                 //std::nth_element(top_dist[i].begin()+count_from_middle_start,top_dist[i].begin()+count_from_middle_end-1,top_dist[i].end());
//                 for(size_t j=0;j<count_from_middle;j++)
//                 {
//                     //std::cout<<"加入"<<top_dist[i][count_from_front+j].second<<std::endl;
//                     top_level_distances[i][count_from_front+j]=top_dist[i][count_from_front+j].first;
//                     top_level_indices[i].second[count_from_front+j]=top_dist[i][count_from_front+j].second;
//                 }
//                 //std::cout<<"中间部分"<<std::endl;

//                 //std::sort(top_dist[i].begin()+count_from_middle_end-1,top_dist[i].end());
//                 for(size_t j=0;j<count_from_back;j++)
//                 {
//                     //std::cout<<"加入"<<top_dist[i][count_from_middle_end+j-1].second<<std::endl;
//                     top_level_distances[i][count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].first;
//                     top_level_indices[i].second[count_from_front+count_from_middle+j]=top_dist[i][count_from_middle_end+j-1].second;
//                 }
//                 //std::cout<<"后面部分"<<std::endl;

//             }
//             //for(int j=0;j<top_paths_num;j++)
//             //{
//             //    std::cout<<top_level_points[j]<<" ";
//             //}
//             //std::cout<<std::endl;
//             std::cout<<"top_level_graph_distances_and_indices构建完成"<<std::endl;
//         }




//         //构建graph*********************************************************************************
//         void naive_construction(const string & algorithm, const std::vector<path>& paths, const size_t k)
//         {
 
//             int paths_num=paths.size();

//             // 自动计算网格大小，将整个数据集平均分为1000个网格
//             //int num_grids = 10;
//             //auto total_grid = calculate_grid_size(paths, num_grids);

//             // 创建R-tree索引
//             // 插入轨迹数据到R-tree中

//             auto start_build_rtree=std::chrono::high_resolution_clock::now();
//             rtree_t rtree;

// ////            #pragma omp parallel
// ////            {     
// ////            int num_threads=omp_get_num_threads();
// ////            int thread_id=omp_get_thread_num();
// ////            int chunk_size=(paths_num+num_threads-1)/num_threads;
// ////            int start = thread_id * chunk_size;
// ////            int end = std::min((thread_id + 1) * chunk_size, paths_num);

//             ////#pragma omp parallel for
//             for (unsigned i = 0; i < paths.size(); ++i) 
//             {   
//                 const auto& current_traj=paths[i];
//                 rtree.insert(std::make_pair(to_traj_grid_range(current_traj), i));
//             }
//             auto end_build_rtree=std::chrono::high_resolution_clock::now();
//             std::chrono::duration<double> build_rtree=end_build_rtree-start_build_rtree;
//             std::cout<<"build rtree time: "<<build_rtree.count()<<" seconds\n"<<std::endl;
//             std::cout<<"rtree_build执行完成"<<std::endl;




//             std::cout<<"开始执行grid_top_level_points"<<std::endl;
//             grid_top_level_points_build(algorithm,paths,rtree);
//             std::cout<<"执行grid_top_level_points完成"<<std::endl;




//             //indices=load_from_file<int>("/home/KNN/orders_data/trajectory_new/saved_graph_index3000_front_back_test");
//             //distances=load_from_file<double>("/home/KNN/orders_data/trajectory_new/saved_graph_dist3000_front_back_test");

//             std::cout<<"开始构建底层图"<<std::endl;

// ////            #pragma omp barrier

//             auto start_build_graph_node=std::chrono::high_resolution_clock::now();
            
//             // 构建kNN图
//             //得到前0.8k个元素
//             int portion_front=0.9*k;
//             int portion_back=0.1*k;
//             //#pragma omp parallel for schedule(dynamic)
//             #pragma omp parallel for 
//             for(size_t i=0;i<paths_num;++i)
//             {
//                 //std::cout<<"第 "<<i<<" 条轨迹"<<std::endl;

//                 const auto& current_traj=paths[i];

//                 std::vector<value> result_grid;
//                 box_t current_traj_box=to_traj_grid_range(current_traj);

                
//                 //底层图最近邻过滤的方式********************************************************
//                 //1.intersects
//                 // rtree.query(bgi::intersects(current_traj_box),std::back_inserter(result_grid));
//                 // double Mul_grid=1.1;
//                 // while(result_grid.size()<270)
//                 // {
//                 //     current_traj_box=to_traj_Mul_grid_range(current_traj,Mul_grid);
//                 //     result_grid.clear();
//                 //     rtree.query(bgi::intersects(current_traj_box),std::back_inserter(result_grid));
//                 //     Mul_grid+=0.1;
//                 // }
//                 // std::cout<<"执行到这里"<<std::endl;
//                 //2.within
//                 // rtree.query(bgi::within(current_traj_box),std::back_inserter(result_grid));
//                 // double Mul_grid=1.1;
//                 // while(result_grid.size()<800)
//                 // {   
//                 //     //std::cout<<"搜到了"<<Mul_grid<<"倍grid"<<std::endl;
//                 //     current_traj_box=to_traj_Mul_grid_range(current_traj,Mul_grid);
//                 //     result_grid.clear();
//                 //     rtree.query(bgi::within(current_traj_box),std::back_inserter(result_grid));
//                 //     Mul_grid+=0.1;
//                 // }

//                 //3.nearest
//                 rtree.query(bgi::nearest(current_traj_box,900),std::back_inserter(result_grid));


//                 //应该用indices、distances代替;
//                 std::vector<std::pair<double,size_t>> front_knn_distances;
//                 std::vector<std::pair<double,size_t>> back_knn_distances;
//                 std::set<int> selected_ids;
//                 std::vector<int> remaining_ids;

//                 if(i%100==0)
//                 {
//                     std::cout<<"front:  这条轨迹box范围涉及到并计算的轨迹数为 "<<result_grid.size()<<endl;
//                 }
//                 //std::cout<<"这条轨迹第一个grid范围涉及到并计算的轨迹数为 "<<result_grid.size()<<endl;
                        
//                 selected_ids.insert(i);
                
//                 //#pragma omp parallel for schedule(dynamic)


//                 // std::cout<<"当前轨迹与其他轨迹之间的距离"<<std::endl;
//                 // int not_zero=0;
//                 // for(int tt=1;tt<paths_num;tt++)
//                 // {
//                 //     double path_distance=minSubTrajectory2(algorithm,current_traj,paths[tt]);
//                 //     if(path_distance<0)
//                 //     {
//                 //         not_zero+=1;
//                 //     }
//                 //     std::cout<<"第 1 条轨迹与第 "<<tt+1<<" 条轨迹之间的距离为: "<<path_distance<<std::endl;
//                 // }
//                 // std::cout<<"与当前轨迹相似的轨迹个数为"<<not_zero<<std::endl;

//                 #pragma omp parallel for
//                 for(const auto& res:result_grid)
//                 {
//                     unsigned j=res.second;
//                     selected_ids.insert(j);
//                     if(i!=j)
//                     {
//                         double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
                        
//                         //#pragma omp critical
//                         front_knn_distances.push_back(std::make_pair(result.first,std::make_pair(j,std::make_pair(result.second.first.first,result.second.first.second))));
//                     }
//                 }
//                 //std::cout<<"执行到这里"<<std::endl;

//                 std::sort(front_knn_distances.begin(),front_knn_distances.end());


//         ////        for(int t=0;t<paths_num;t++)
//         ////        {
//         ////            if(selected_ids.find(t)==selected_ids.end())
//         ////            {
//         ////                remaining_ids.push_back(t);
//         ////            }
//         ////        }
//         ////        //std::cout<<"执行到这里"<<std::endl;
//         ////
//         ////
//         ////
//         ////        if(front_knn_distances.size()>=150)
//         ////        {
//         ////            for(size_t ttt=150;ttt<front_knn_distances.size();ttt++)
//         ////            {
//         ////                remaining_ids.push_back(front_knn_distances[ttt].second);
//         ////            }
//         ////        }
//         ////        
//         ////        std::cout<<"执行到这里"<<std::endl;
//         ////        
//         ////        std::random_device rd;
//         ////        std::mt19937 g(rd());
//         ////        std::shuffle(remaining_ids.begin(),remaining_ids.end(),g);
//         ////        std::vector<int> random_indices(remaining_ids.begin(),remaining_ids.begin()+portion_back);
//         ////        std::cout<<"执行到这里"<<std::endl;
//         ////        //std::cout<<"执行到这"<<std::endl;
//         ////       
//         ////        //#pragma omp parallel for schedule(dynamic)
//         ////        #pragma omp parallel for
//         ////        for(int t=0;t<random_indices.size();t++)
//         ////        {
//         ////            unsigned j=random_indices[t];
//         ////            double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
//         ////            //#pragma omp critical
//         ////            back_knn_distances.push_back(std::make_pair(path_distance,j));
//         ////        }
//         ////        std::cout<<"执行到这里"<<std::endl;


//                 // if(remaining_ids.size()<50)
//                 // {
//                 //     size_t num_to_get=100;
//                 //     size_t start_index = std::max(0ul, front_knn_distances.size() - num_to_get);
//                 //     for(size_t ttt=start_index;ttt<front_knn_distances.size();++ttt)
//                 //     {
//                 //         back_knn_distances.push_back(front_knn_distances[ttt]);
//                 //     }
//                 // }
//                 // else 
//                 // {
//                 //     std::random_device rd;
//                 //     std::mt19937 g(rd());
//                 //     std::shuffle(remaining_ids.begin(),remaining_ids.end(),g);
//                 //     std::vector<int> random_indices(remaining_ids.begin(),remaining_ids.begin()+portion_back);

//                 //     //std::cout<<"执行到这"<<std::endl;
                
//                 //     //#pragma omp parallel for schedule(dynamic)
//                 //     #pragma omp parallel for
//                 //     for(int t=0;t<random_indices.size();t++)
//                 //     {
//                 //         unsigned j=random_indices[t];
//                 //         double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
//                 //         //#pragma omp critical
//                 //         back_knn_distances.push_back(std::make_pair(path_distance,j));
//                 //     }
//                 //     //std::cout<<"执行到这"<<std::endl;


//                 // }
               
                

                
//         ////        std::sort(back_knn_distances.begin(),back_knn_distances.begin());

//                 //std::cout<<"执行到这"<<std::endl;




//                 //random nearests
//                 // std::random_device rd_;
//                 // std::mt19937 g_(rd_());
//                 // std::shuffle(g_);
//                 // std::vector<int> random_nearest();

//             //    int res_num=front_knn_distances.size();
//             //    std::cout<<"within "<<res_num<<" 条轨迹与这个data轨迹之间的距离: "<<std::endl;
//             //    for(size_t j=0;j<res_num;++j)
//             //    {
//             //        std::cout<<"这个data轨迹与within "<<j+1<<" 索引号为 "<<front_knn_distances[j].second<<" 轨迹之间的距离 "<<front_knn_distances[j].first<<std::endl;
//             //    }
//             //    std::cout<<std::endl;
//             //    int res_num1=0;
//                 for(size_t j=0;j<150;++j)
//                 {
//                     distances[i][j]=front_knn_distances[j].first;
//                     indices[i][j]=front_knn_distances[j].second;
//             //       std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这个data轨迹之间的距离: "<<distances[i][j]<<std::endl;
//             //       res_num1++;
//                 }
//             //     std::cout<<"执行到这里"<<std::endl;
//             //     for(size_t j=0;j<20;++j)
//             //     {
//             //         distances[i][60+j]=front_knn_distances[60+4*j].first;
//             //         indices[i][60+j]=front_knn_distances[60+4*j].second;
//             // //        std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这个data轨迹之间的距离: "<<distances[i][60+j]<<std::endl;
//             // //       res_num1++;
//             //     }
//             //     std::cout<<"执行到这里"<<std::endl;




//                 // for(size_t j=0;j<portion_front;++j)
//                 // {
//                 //     distances[i][j]=front_knn_distances[j].first;
//                 //     indices[i][j]=front_knn_distances[j].second;
//                 // }
//                 //std::cout<<"执行到这"<<std::endl;
//         ////        for(size_t j=0;j<portion_back;++j)
//         ////        {
//         ////            distances[i][portion_front+j]=back_knn_distances[j].first;
//         ////            indices[i][portion_front+j]=back_knn_distances[j].second;
//             //       std::cout<<"加入第"<<res_num1+1<<"条轨迹边,索引号为 "<<front_knn_distances[j].second<<" 与这条data轨迹之间的距离为: "<<distances[i][portion_front+j]<<std::endl;
//             //       res_num1++;
//         ////        }
//         ////        std::cout<<"执行到这里"<<std::endl;
//                //exit(1);
//             ////    std::cout<<"graph_node_build执行完成"<<std::endl;
//             ////    std::cout<<"graph_node_build执行完成"<<std::endl;
// ////            }                
//             }

//             //indices=load_from_file<int>("/home/KNN/orders_data/trajectory_new_chengdu/saved_graph_index5000_front_back_test_final_final");
//             //distances=load_from_file<double>("/home/KNN/orders_data/trajectory_new_chengdu/saved_graph_dist5000_front_back_test_final_final");
            

//             //3.过滤+类似最开始的顶层图
//             // std::cout<<"开始创建top-level-graph"<<std::endl;
//             // auto start_build_top_level_graph=std::chrono::high_resolution_clock::now();
//             // //front_knn_distances的大小一定会大于2*k/3*k或者是自己设置的nearest的值，可以根据这个和想要消掉的邻居数设置drop_neighbors_num;
//             // int drop_neighbors=0.9*k;
//             // std::vector<char> top_level_points_visitpool(paths_num,false);
//             // std::unordered_set<int> top_level_points_set;
//             // for(int t=0;t<paths_num;t++)
//             // {
//             //     top_level_points_set.insert(t);
//             // }

//             // int add_top_level_point=rand()%paths_num;
//             // top_level_points.push_back(add_top_level_point);

//             // top_level_points_visitpool[add_top_level_point]=true;
//             // top_level_points_set.erase(add_top_level_point);
            
//             // while(!top_level_points_set.empty())
//             // {
//             //     //std::cout<<"add_top_level_point "<<add_top_level_point<<std::endl;
//             //     for(int j=0;j<drop_neighbors;++j)
//             //     {
//             //         top_level_points_visitpool[indices[add_top_level_point][j]]=true;
//             //         top_level_points_set.erase(indices[add_top_level_point][j]);
//             //     }
//             //     auto it=std::next(top_level_points_set.begin(),std::rand()%top_level_points_set.size());
//             //     add_top_level_point=*it;
//             //     top_level_points.push_back(add_top_level_point);
                    
//             //     top_level_points_visitpool[add_top_level_point]=true;
//             //     top_level_points_set.erase(add_top_level_point);
//             //     //std::cout<<"又选了一个add_top_level_point"<<std::endl;
//             // }


//             // for(int t=0;t<top_level_points.size();t++)
//             // {
//             //     std::cout<<top_level_points[t]<<" ";
//             // }
//             // std::cout<<std::endl;

//             // int top_paths_num=top_level_points.size();
//             // std::cout<<"top_level_points_set的size "<<top_paths_num<<std::endl;
//             // top_level_indices.resize(top_paths_num);
//             // top_level_distances.resize(top_paths_num);
//             // for(int t=0;t<top_paths_num;t++)
//             // {
//             //     top_level_indices[t].first=top_level_points[t];
//             //     top_level_indices[t].second.resize(top_paths_num-1);
//             //     top_level_distances[t].resize(top_paths_num-1);
//             // }

//             // std::sort(top_level_points.begin(),top_level_points.end());

//             // for(int t=0;t<top_paths_num;t++)
//             // {                
//             //     for(int tt=t+1;tt<top_paths_num;tt++)
//             //     {
//             //         //std::cout<<"执行到t="<<t<<" tt="<<tt<<std::endl;
//             //         double top_paths_distance=minSubTrajectory2(algorithm,paths[top_level_points[t]],paths[top_level_points[tt]]);
//             //         top_level_distances[t][tt-1]=top_paths_distance;
//             //         top_level_indices[t].second[tt-1]=top_level_points[tt];
//             //         //std::cout<<"执行到此"<<std::endl;

//             //         top_level_distances[tt][t]=top_paths_distance;
//             //         //std::cout<<"执行到此"<<std::endl;
//             //         top_level_indices[tt].second[t]=top_level_points[t];
//             //         //std::cout<<"执行到此"<<std::endl;
//             //     }

//             // }
//             // auto end_build_top_level_graph=std::chrono::high_resolution_clock::now();
//             // std::chrono::duration<double> build_top_level_graph=end_build_top_level_graph-start_build_top_level_graph;
//             // std::cout<<"build top_level_graph time: "<<build_top_level_graph.count()<<" seconds\n"<<std::endl;
//             // std::cout<<"top_level_graph_build构建完成"<<std::endl;
                









//             auto end_build_graph_node=std::chrono::high_resolution_clock::now();
//             std::chrono::duration<double> build_graph_node=end_build_graph_node-start_build_graph_node;
//             std::cout<<"build graph_node time: "<<build_graph_node.count()<<" seconds\n"<<std::endl;
//             std::cout<<"graph_node_build执行完成"<<std::endl;
//         }

//不分类***************************************************************************************************************








////                // Calculate number of nearest neighbors to consider
////                size_t front_num_nearest = std::max(static_cast<size_t>(portion_front), static_cast<size_t>(1)); // Ensure at least 1 neighbor
////                front_num_nearest = std::min(front_num_nearest, front_knn_distances.size() + 1); // Cap at available neighbors plus one
////
////                // Ensure we have at least 0.8 * k nearest neighbors
////                if(front_knn_distances.size() < front_num_nearest) 
////                {
////                    // Query more distant trajectories
////                    result_grid.clear(); // Clear the previous results
////                    //rtree_point grid_index=to_grid_range(centroid,grid_size*2);
////                    current_traj_box=to_traj_2grid_range(current_traj);
////
////                    rtree.query(bgi::intersects(current_traj_box), std::back_inserter(result_grid));
////
////                    if(i%10==0)
////                    {
////                        std::cout<<"front:  这条轨迹第二个grid范围涉及到的轨迹数为 "<<result_grid.size()<<endl;
////                    }
////
////            ////        std::cout<<"这条轨迹第二个grid范围涉及到的轨迹数为 "<<result_grid.size()<<endl;
////                    bool found_new_neighbor = false;
////                    for (const auto& res : result_grid) 
////                    {
////                        unsigned j=res.second;
////                        if (j != i && std::find_if(front_knn_distances.begin(), front_knn_distances.end(),
////                                            [&](const std::pair<double, size_t>& p) { return p.second == j; }) == front_knn_distances.end()) 
////                        {
////                            double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
////                            front_knn_distances.push_back(std::make_pair(path_distance,j));
////                            found_new_neighbor = true;
////                            if(front_knn_distances.size()>=portion_front)
////                                break; // Exit the loop after finding one new neighbor
////                        }
////                    }
////                }
////
////                if(i%10==0)
////                {
////                    std::cout<<"front:  这条轨迹总共计算的轨迹数为 "<<front_knn_distances.size()<<endl;
////                }
////            ////    std::cout<<"这条轨迹总共计算的轨迹数为 "<<front_knn_distances.size()<<endl;
////
////                // 对kNN图按照距离排序，并输出结果
////                std::sort(front_knn_distances.begin(), front_knn_distances.end());
////
////                // Now ensure we have exactly front_num_nearest neighbors
////                if (front_knn_distances.size() > portion_front) 
////                {
////                    front_knn_distances.resize(portion_front);
////                }
////
////
////            //得到后0.2k个元素
////                int portion_back=0.2*k;
////                std::vector<value> result_far_grid;
////                box_t current_traj_far_box=to_traj_far_grid_range(current_traj,total_grid);
////
////                rtree.query(bgi::intersects(current_traj_far_box),std::back_inserter(result_far_grid));
////
////                //应该用indices、distances代替;
////                std::vector<std::pair<double,size_t>> back_knn_distances;
////                if(i%10==0)
////                {
////                    std::cout<<"back:  这条轨迹第一个grid范围涉及到并计算的轨迹数为 "<<result_far_grid.size()<<endl;
////                }
////            ////    std::cout<<"这条轨迹第一个grid范围涉及到并计算的轨迹数为 "<<result_far_grid.size()<<endl;
////                for(const auto& far_res:result_far_grid)
////                {
////                    unsigned j=far_res.second;
////                    if(i!=j)
////                    {
////                        double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
////                        back_knn_distances.push_back(std::make_pair(path_distance,j));
////                    }
////                }
////        
////                // Calculate number of nearest neighbors to consider
////                size_t back_num_nearest = std::max(static_cast<size_t>(portion_back), static_cast<size_t>(1)); // Ensure at least 1 neighbor
////                back_num_nearest = std::min(back_num_nearest, back_knn_distances.size() + 1); // Cap at available neighbors plus one
////
////                // Ensure we have at least 0.8 * k nearest neighbors
////                if(back_knn_distances.size() < back_num_nearest) 
////                {
////                    // Query more distant trajectories
////                    result_far_grid.clear(); // Clear the previous results
////                    //rtree_point grid_index=to_grid_range(centroid,grid_size*2);
////                    current_traj_far_box=to_traj_far_2grid_range(current_traj,total_grid);
////
////                    rtree.query(bgi::intersects(current_traj_far_box), std::back_inserter(result_far_grid));
////
////                    if(i%10==0)
////                    {
////                        std::cout<<"back:  这条轨迹第二个grid范围涉及到的轨迹数为 "<<result_far_grid.size()<<endl;
////                    }
////
////            ////        std::cout<<"这条轨迹第二个grid范围涉及到的轨迹数为 "<<result_far_grid.size()<<endl;
////                    bool found_new_neighbor = false;
////                    for (const auto& far_res : result_far_grid) 
////                    {
////                        unsigned j=far_res.second;
////                        if (j != i && std::find_if(back_knn_distances.begin(), back_knn_distances.end(),
////                                            [&](const std::pair<double, size_t>& p) { return p.second == j; }) == back_knn_distances.end()) 
////                        {
////                            double path_distance=minSubTrajectory2(algorithm,current_traj,paths[j]);
////                            back_knn_distances.push_back(std::make_pair(path_distance,j));
////                            found_new_neighbor = true;
////                            if(back_knn_distances.size()>=portion_back)
////                                break; // Exit the loop after finding one new neighbor
////                        }
////                    }
////                }
////
////                if(i%10==0)
////                {
////                    std::cout<<"back:  这条轨迹总共计算的轨迹数为 "<<back_knn_distances.size()<<endl;
////                }
////            ////    std::cout<<"这条轨迹总共计算的轨迹数为 "<<knn_distances.size()<<endl;
////
////                // 对kNN图按照距离排序，并输出结果
////                std::sort(back_knn_distances.begin(), back_knn_distances.end());
////
////                // Now ensure we have exactly back_num_nearest neighbors
////
////                //替换成distances、indices
////                distances[i].resize(k);
////                indices[i].resize(k);
////                for(size_t j=0;j<portion_front;++j)
////                {
////                    distances[i][j]=front_knn_distances[j].first;
////                    indices[i][j]=front_knn_distances[j].second;
////                }
////
////                auto it=back_knn_distances.end()-portion_back;
////                size_t j=0.8*k;
////                for(; it != back_knn_distances.end() && j <= k; ++it, ++j)
////                {
////                    distances[i][j]=back_knn_distances[j].first;
////                    indices[i][j]=back_knn_distances[j].second;
////                }
////
////                if(i%10==0)
////                {
////                    std::cout<<"这条轨迹创建 "<<distances[i].size()<<" 条边"<<std::endl;
////                }
////            ////    std::cout<<"这条轨迹创建 "<<knn_distances.size()<<" 条边"<<std::endl;
////            }
////        }




    ////    void nth_index_element(const std::vector<double>& dist, const int length, std::vector<double>& distances_k, std::vector<int>& indices_k, const size_t n)
    ////    {
            //construct a vector with index
            //std::vector<std::pair<double, int>> v;
            //for(int i=0;i<length;++i)
            //    v.push_back(std::pair<double, int>(dist[i], i));

            //since it contain the point itself and find n+1 minimum elements
            //std::nth_element(v.begin(), v.begin()+n+1, v.end());
            //std::sort(v.begin(), v.begin()+n+1);
            //eliminate the first element(itself)
            //for(size_t i=0;i<n;++i)
            //{
            //    distances_k[i] = v[i+1]->first;
            //    indices_k[i] = v[i+1]->second;
            //}
            //std::cout<<"nth_index_element执行完毕"<<endl;


            //取length中的前0.8*n和后0.2*n;
            //取length中的前0.6*n,中间的0.2*n和后0.2*n;

    ////        std::vector<std::pair<double,int>> v;
    ////        for(int i=0;i<length;++i)
    ////            v.push_back(std::pair<double,int>(dist[i],i));

            //double portion_from_front=0.8;
    ////        double portion_from_front=0.6;
    ////        double portion_from_middle=0.2;
    ////        double portion_from_back=0.2;

    ////        int count_from_front=static_cast<int>(n*portion_from_front);
    ////        int count_from_middle=static_cast<int>(n*portion_from_middle);
    ////        int count_from_back=static_cast<int>(n*portion_from_back);
            
    ////        std::nth_element(v.begin(),v.begin()+count_from_front+1,v.end());
    ////        std::sort(v.begin(),v.begin()+count_from_front+1);
    ////        for(size_t i=0;i<count_from_front;++i)
    ////        {
    ////            distances_k[i]=v[i+1].first;
    ////            indices_k[i]=v[i+1].second;
    ////        }


    ////        int count_from_middle_start=(length-count_from_front-count_from_back)/2-(count_from_middle/2)+count_from_front;
    ////        int count_from_middle_end=(length-count_from_front-count_from_back)/2+(count_from_middle/2)+count_from_front;
    ////        std::nth_element(v.begin(),v.begin()+count_from_middle_start+1,v.end());
    ////        std::nth_element(v.begin(),v.begin()+count_from_middle_end+1,v.end());
    ////        std::sort(v.begin()+count_from_middle_start+1,v.begin()+count_from_middle_end+1);
    ////        for(size_t i=0;i<count_from_middle_end-count_from_middle_start;++i)
    ////        {
    ////            distances_k[count_from_front+i]=v[count_from_middle_start+i+1].first;
    ////            indices_k[count_from_front+i]=v[count_from_middle_start+i+1].second;
    ////        }


    ////        std::nth_element(v.begin(),v.begin()+length-count_from_back,v.end());
    ////        std::sort(v.begin()+length-count_from_back,v.end());
            //82分
            //for(size_t i=0;i<count_from_back;++i)
            //{  
            //    distances_k[count_from_front+i]=v[length-count_from_back+i].first;
            //   indices_k[count_from_front+i]=v[length-count_from_back+i].second;
            //}
            //622分
    ////        for(size_t i=0;i<count_from_back;++i)
    ////        {
    ////            distances_k[count_from_front+count_from_middle+i]=v[length-count_from_back+i].first;
    ////            indices_k[count_from_front+count_from_middle+i]=v[length-count_from_back+i].second;
    ////        }

    ////    }

        /*
         * @brief use a naive construction way to build a knn nearest neighbor graph
         * @param data: the coordinate of the point in matrix type with shape [point_num, dim]
         */




    ////    void naive_construction(const string & algorithm, const std::vector<path>& paths, const size_t k)
    ////    {
    ////        int paths_num = paths.size();

            //自己添加的内容
            //int paths_num=3000;

    ////        std::vector<std::vector<double>> dist(paths_num,std::vector<double>(paths_num));


            //compute the distance between each two points
    ////        for(int i=0;i<paths_num;++i)
    ////        {
    ////            int k=0;
    ////            for(int j=i+1;j<paths_num;++j)
    ////            {
    ////                k++;
    ////                path v1 = paths[i];
    ////                path v2 = paths[j];
                    //std::cout<<"naive_construction中for循环里paths[i]/paths[j]的长度"<<endl;
                    //std::cout<<paths[i].size()<<endl;
                    //std::cout<<paths[j].size()<<endl;
                    //std::cout<<"naive_construction中for循环里v1/v2的长度"<<endl;
                    
                    //std::cout<<"v1的size "<<v1.size()<<" "<<"v2的size "<<v2.size()<<endl;
                    //std::cout<<v1.size()<<endl;
                    //std::cout<<v2.size()<<endl;
                    
                    //std::cout<<"执行minSubTrajectory"<<endl;

                    //effifientAlgorithm时间复杂度太大无法接受
                    //double path_distance=execute_data_to_data("efficientAlgorithm",v1,v2);
                    //dist[i][j]=dist[j][i]=path_distance;


                    //data-data原始版
                    //double path_distance=minSubTrajectory1(algorithm, v1, v2);
                    
                    //data-data改进版
    ////                double path_distance=minSubTrajectory2(algorithm, v1, v2);
                    
                    //std::cout<<"执行结束minSubTrajectory"<<endl;
    ////                dist[i][j] = dist[j][i] = path_distance;
    ////            }
                //std::cout<<"k的值为 "<<k<<endl;
    ////        }
    ////        save_dist_to_file(dist,"/home/KNN/orders_data/trajectory_new/saved_graph_index30000_622_top_dist_300");

            //std::cout<<"执行nth_index_element"<<endl;
    ////        for(int i=0;i<paths_num;++i)
    ////        {
    ////            nth_index_element(dist[i], paths_num, distances[i], indices[i], k);
    ////        }


            //保存创建的graph
            //save_to_file<int>(indices, "/home/KNN/orders_data/trajectory_new/saved_graph_index3000_622_top_pre_300");
            //save_to_file<double>(distances, "/home/KNN/orders_data/trajectory_new/saved_graph_dist3000_622_top_pre_300");


            //下载dist,indices,distance
        ////    dist=load_2d_double_vector_from_file("/home/KNN/orders_data/trajectory_new/saved_graph_index3000_622_top_dist_300");
        ////    indices=load_from_file<int>("/home/KNN/orders_data/trajectory_new/saved_graph_index3000_622_top_300");
        ////    distances=load_from_file<double>("/home/KNN/orders_data/trajectory_new/saved_graph_dist3000_622_top_300");

        ////    top_level_points=load_top_level_points_from_file("/home/KNN/orders_data/trajectory_new/top_level_graph_points_300");
        ////    top_level_indices=load_top_level_indices_from_file("/home/KNN/orders_data/trajectory_new/top_level_graph_indices_300");



    ////        std::cout<<"开始创建top-level-grpah"<<endl;
            // //top-level-graph选择上层graph的节点集合
            // //加入点的邻居不再考虑加入
    ////        int drop_neighbors_num=k;
            // //得到上层graph中的节点集合top_level_points
    ////        this->top_level_points_build(dist,paths_num,drop_neighbors_num);
    ////        std::cout<<"top_level_points_build构建完成"<<endl;

    ////        std::sort(top_level_points.begin(),top_level_points.end());

    ////        int top_level_paths_num=top_level_points.size();
    ////        std::vector<std::pair<int,std::vector<std::pair<double,int>>>> top_level_points_dist(top_level_paths_num);
            
    ////        for(auto& element1:top_level_points_dist)
    ////        {
    ////            element1.second.resize(top_level_paths_num);
    ////        }
            
             //得到上层graph的节点之间的距离dist
    ////        for(auto top_level_point_row:top_level_points)
    ////        {
    ////            std::pair<int,std::vector<std::pair<double,int>>> element2;
    ////           element2.first=top_level_point_row;

    ////            for(auto top_level_point_line:top_level_points)
    ////            {
    ////                element2.second.push_back({dist[top_level_point_row][top_level_point_line],top_level_point_line});
    ////            }
    ////            top_level_points_dist.push_back(element2);
    ////        }
    ////        std::cout<<"执行到此1"<<endl;
            //对行按照节点索引进行排序
    ////         std::sort(top_level_points_dist.begin(),top_level_points_dist.end());

    ////        for(auto& element_:top_level_points_dist)
    ////        {
    ////            std::sort(element_.second.begin(),element_.second.end());
                //保留上层graph中每个节点的距离后30%的邻居来构建上层graph
    ////            std::vector<int> percent30_back;
    ////            int num_to_keep=static_cast<int>(0.3*element_.second.size());
    ////            for(auto rit=element_.second.rbegin();rit!=(element_.second.rbegin()+num_to_keep);++rit)
    ////            {
    ////                percent30_back.push_back(rit->second);
    ////            }
                //std::cout<<"添加的element_.first "<<element_.first<<endl;
                //std::cout<<"percent30_back ";
                //for(auto n:percent30_back)
                //{
                //    std::cout<<n<<" ";
                //}
                //std::cout<<endl;
    ////            top_level_indices.emplace_back(element_.first,percent30_back);
    ////        }
    ////        std::cout<<"top_level_indices构建完成"<<endl;
            // // std::cout<<"top_level_indices的size "<<top_level_indices.size()<<" "<<"top_level_points集合"<<std::endl;
            // // for(const auto& pair_vector:top_level_indices)
            // // {
            // //     std::cout<<"top_level_point "<<pair_vector.first<<" 它的邻居节点 ";
            // //     for(int value_vector:pair_vector.second)
            // //     {
            // //         std::cout<<value_vector<<" ";
            // //     }
            // //     std::cout<<std::endl;     
            // // }


    ////        std::cout<<"保存创建的top_level_graph"<<endl;
    ////    }

    public:

        /*Default constructor*/
        Knn_Graph(){}

        ~Knn_Graph()
        {
            //if(indices.ptr())
            //    delete[] indices.ptr();
            //if(distances.ptr())
            //    delete[] distances.ptr();
        }

        /*
         *@brief use the data to build a k nearest graph
         *@param data: the points data in shape [point_num, vec_len]
         *@param method: the method to build the graph
         *@param k: the param k of knn graph
         */
        void build_graph(const string & algorithm, const std::vector<path>& paths, const size_t k, BUILD_GRAPH_METHOD method)
        {
            size_t path_num = paths.size();
            indices.resize(path_num,std::vector<int>(k));
            distances.resize(path_num,std::vector<double>(k));
            //分类
            start_end_positions.resize(path_num,std::vector<std::pair<int,int>>(k));

            
            this->k = k;
            //std::cout<<"执行naive_construction"<<endl;


            // //(DTW、EDR、ERP)+(CMA、POS、PSS)、DTSM计算一对轨迹相似的时间
            // auto start_compute_trajectories_pair=std::chrono::high_resolution_clock::now();
            // int paths_num=paths.size();
            // int t=0;
            // for(size_t i=1;i<=100;i+=10)
            // {
            //     //Result result=minSubTrajectory2(algorithm,paths[i],paths[i+1]);
            //     t++;
            //     std::cout<<"计算第"<<t<<"条";

            //     double res=100000000000.0;
            //     const path& current_trajectory=paths[i];
            //     const path& target_trajectory=paths[i+1];
            //     double sum=0;
            //     for(size_t start=0;start<current_trajectory.size()-1;start++)
            //     {
            //         for(size_t end=start+1;end<current_trajectory.size();end++)
            //         {
            //             sum++;
            //             path sub_path(current_trajectory.begin()+start, current_trajectory.begin()+end);
            //             subResult r=execute("efficientAlgorithm",sub_path,target_trajectory);
            //             if(res>r.second)
            //             {
            //                 res=r.second;
            //             }
            //         }
            //     }
            //     std::cout<<sum<<std::endl;
            // }
            // auto end_compute_trajectories_pair=std::chrono::high_resolution_clock::now();
            // std::chrono::duration<double> elapsed_compute_trajectories_pair=end_compute_trajectories_pair-start_compute_trajectories_pair;
            // std::cout<<"Elpased time: "<<elapsed_compute_trajectories_pair.count()<<" seconds\n"<<endl;
            // exit(1);


            if(method==NAIVE)
            {
                naive_construction(algorithm, paths, k);
            }
        }

        /*
         * @brief load the graph from a saved graph file (index and distance)
         * @param index_file_name: the file path of the index file
         * @param dist_file_path: the file path of the distance file
         * @exception define in the io.h
         */

        //分类
        void load_graph(const std::string& graph_index_path,const std::string& graph_dist_path,const std::string& graph_start_end_position_path,const std::string& top_level_graph_points_path,const std::string& top_level_graph_indices_path,const std::string& top_level_graph_distances_path,const std::string& top_level_graph_start_end_positions_path)
        {
            try{
                indices=load_from_file<int>(graph_index_path);
                // std::cout<<indices[0].size()<<std::endl;
                // for(int i=0;i<indices.size();i++)
                // {
                //     std::cout<<i<<" 这一个点的邻居索引"<<" ";
                //     for(int j=0;j<indices[i].size();j++)
                //     {
                //         std::cout<<indices[i][j]<<" ";
                //     }
                //     std::cout<<std::endl;
                // }
                // exit(1);
                distances=load_from_file<double>(graph_dist_path);
                start_end_positions=load_start_end_positions_from_file(graph_start_end_position_path);
                std::cout<<"执行到此8"<<endl;
                top_level_points=load_top_level_points_from_file(top_level_graph_points_path);
                std::cout<<"执行到此9"<<endl;
                top_level_indices=load_top_level_indices_from_file(top_level_graph_indices_path);
                top_level_distances=load_from_file<double>(top_level_graph_distances_path);
                top_level_start_end_positions=load_top_level_start_end_positions_from_file(top_level_graph_start_end_positions_path);
                std::cout<<"执行到此10"<<endl;

                k=indices[0].size();
            }catch(std::exception& e)
            {
                throw e;
            }
        }
        void save_graph(const std::string& index_file_path, const std::string& dist_file_path, const std::string& start_end_position_file_path, const std::string& top_level_points_file_path, const std::string& top_level_indices_file_path, const std::string& top_level_distances_file_path, const std::string& top_level_start_end_positions_file_path)
        {
            try{
                save_to_file<int>(indices,index_file_path);
                save_to_file<double>(distances,dist_file_path);
                save_start_end_positions_to_file(start_end_positions,start_end_position_file_path);
                std::cout<<"执行到此2"<<endl;
                save_top_level_points_to_file(top_level_points,top_level_points_file_path);
                std::cout<<"执行到此3"<<endl;
                save_top_level_indices(top_level_indices,top_level_indices_file_path);
                save_to_file<double>(top_level_distances,top_level_distances_file_path);
                std::cout<<"执行到此4"<<endl;
                save_top_level_start_end_positions_to_file(top_level_start_end_positions,top_level_start_end_positions_file_path);
                }catch(std::exception& e)
                {
                    throw e;
                }
        }
        // void load_graph(const std::string& graph_index_path,const std::string& graph_dist_path,const std::string& top_level_graph_points_path,const std::string& top_level_graph_indices_path,const std::string& top_level_graph_distances_path)
        // {
        //     try{
        //         indices = load_from_file<int>(graph_index_path);
        //         std::cout<<indices[0].size()<<std::endl;
        //         distances = load_from_file<double>(graph_dist_path);
        //         std::cout<<"执行到此8"<<endl;
        //         //下载top-level-graph节点和对应节点的邻居索引
        //         top_level_points = load_top_level_points_from_file(top_level_graph_points_path);
        //         std::cout<<"执行到此9"<<endl;
        //         top_level_indices = load_top_level_indices_from_file(top_level_graph_indices_path);
        //         top_level_distances = load_from_file<double>(top_level_graph_distances_path);
        //         std::cout<<"执行到此10"<<endl;

        //         //if(indices.rows != distances.rows || indices.cols != distances.cols)
        //         //{
        //         //    throw GnnsException("Saved graph file error\n");
        //         //}
        //         k = indices[0].size();
        //     }catch (std::exception& e){
        //         throw e;
        //     }
        // }

        // /*
        //  * @brief save the graph into the disk (index and distance)
        //  * @param index_file_name: the file path of the index file
        //  * @param dist_file_path: the file path of the distance file
        //  * @exception define in the io.h
        //  */
        // void save_graph(const std::string& index_file_path, const std::string& dist_file_path, const std::string& top_level_points_file_path, const std::string& top_level_indices_file_path, const std::string& top_level_distances_file_path)
        // {
        //     try{
        //         save_to_file<int>(indices, index_file_path);
        //         save_to_file<double>(distances, dist_file_path);
        //         std::cout<<"执行到此2"<<endl;
        //         //保存top-level-graph中的节点top_level_points和每个对应节点的邻居索引top_level_indices
        //         save_top_level_points_to_file(top_level_points,top_level_points_file_path);
        //         std::cout<<"执行到此3"<<endl;
        //         save_top_level_indices(top_level_indices,top_level_indices_file_path);
        //         save_to_file<double>(top_level_distances,top_level_distances_file_path);
        //         std::cout<<"执行到此4"<<endl;
        //     }catch (std::exception& e){
        //         throw e;
        //     }
        // }

        void get_neighbors(const int search_index,
            std::vector<int> indices,
            std::vector<float> dists,
             int graph_search_expand=-1)
        {
            if(graph_search_expand == -1 || graph_search_expand > k)
            {
                graph_search_expand = k;
            }

            indices.resize(graph_search_expand);
            dists.resize(graph_search_expand);
            for(int i=0;i<graph_search_expand;++i)
            {
                indices[i] = this->indices[search_index][i];
                dists[i] = this->distances[search_index][i];
            }
        }

        std::vector<int> get_neighbors(const int search_index, int graph_search_expand=-1)
        {
            if(graph_search_expand == -1 || graph_search_expand > k)
            {
                graph_search_expand = k;
            }

            std::vector<int> neighbors;
            neighbors.resize(graph_search_expand);
            for(int i=0;i<graph_search_expand;++i)
            {
                neighbors[i] = this->indices[search_index][i];
            }

            return neighbors;
        }

        std::vector<std::pair<int,int>> get_start_end_positions(const int search_index, int graph_search_expand=-1)
        {
            if(graph_search_expand == -1 || graph_search_expand > k)
            {
                graph_search_expand = k;
            }

            std::vector<std::pair<int,int>> start_end_position;
            start_end_position.resize(graph_search_expand);
            for(int i=0;i<graph_search_expand;++i)
            {

                start_end_position[i].first = this->start_end_positions[search_index][i].first;
                start_end_position[i].second = this->start_end_positions[search_index][i].second;
            }
            return start_end_position;
        }

        
        //分别得到top-level-graph的节点数量，和对应index的节点序号
        int get_top_level_points_num()
        {
            return top_level_points.size();
        }
        int get_top_level_point_number_of_index(int top_level_v_it_index)
        {
            return top_level_points[top_level_v_it_index];
        }
        //得到top-level-graph中节点的邻居
        std::vector<int> get_top_level_neighbors(const int search_index)
        {
            std::vector<int> neighbors;
            for(const auto& top_level_point_pair:top_level_indices)
            {
                if(top_level_point_pair.first==search_index)
                {
                    neighbors.resize(top_level_point_pair.second.size());
                    for(int i=0;i<top_level_point_pair.second.size();i++)
                    {
                        neighbors[i]=top_level_point_pair.second[i];
                    }
                    //std::cout<<"添加邻居成功 "<<std::endl;
                    break;
                }
            }
            return neighbors;
        }

        std::vector<std::pair<int,int>> get_top_level_start_end_positions(const int search_index, int graph_search_expand=-1)
        {
            if(graph_search_expand == -1 || graph_search_expand > k)
            {
                graph_search_expand = k;
            }

            std::vector<std::pair<int,int>> start_end_position;
            start_end_position.resize(graph_search_expand);
            for(int i=0;i<graph_search_expand;++i)
            {

                start_end_position[i].first = this->start_end_positions[search_index][i].first;
                start_end_position[i].second = this->start_end_positions[search_index][i].second;
            }
            return start_end_position;
        }

        //得到top-level-graph中节点的索引号top_level_points,对应节点的邻居索引
        std::vector<int>& get_top_level_points()
        {
            return top_level_points;
        }
        std::vector<std::pair<int,std::vector<int>>>&  get_top_level_indices()
        {
            return top_level_indices;
        }
       std::vector<std::vector<std::pair<int,int>>>& get_top_level_start_end_positions()
       {
            return top_level_start_end_positions;
       }
    private:

        //distance of the two points
        std::vector<std::vector<double>> distances;

        //nearst neighbor index
        std::vector<std::vector<int>> indices;

        std::vector<std::vector<std::pair<int,int>>> start_end_positions;


        std::vector<std::vector<double>> reverse_distances;
        std::vector<std::vector<int>> reverse_indices;

        //k nearst neighbor in the graph
        size_t k;

        
        //top-level-graph
        //top-level-graph上层graph中的节点集合,节点索的序号
        std::vector<int> top_level_points;
        //top-level-graph上层graph中对应节点的邻居索引
        std::vector<std::pair<int,std::vector<int>>> top_level_indices;
        std::vector<std::vector<double>> top_level_distances;

        std::vector<std::vector<std::pair<int,int>>> top_level_start_end_positions;
 
        //Distance distance;
    };

//}

#endif //GNNS_KNN_GRAPH_H