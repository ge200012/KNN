#ifndef GNNS_GNNS_INDEX_H
#define GNNS_GNNS_INDEX_H

#include <set>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include "knn_graph.h"
#include "../util/params.h"

#include "mostSimilar.h"

//自己添加的内容
#include <unordered_set>

//namespace gnns
//{
    class Gnns_Params : public Index_Params
    {
    public:
        //Gnns_Params(size_t Graph_k = 1000, BUILD_GRAPH_METHOD method = NAIVE)
        Gnns_Params(size_t Graph_k = 300, BUILD_GRAPH_METHOD method = NAIVE)
        {
            algorithm = "GNNS";
            this->Graph_k = Graph_k;
            this->method = method;
        }
    public:
        size_t Graph_k;
        BUILD_GRAPH_METHOD method;

    };

    template<typename Distance>
    class Gnns_Index
    {
        typedef typename Distance::ElementType ElementType;
        typedef typename Distance::DistanceType DistanceType;
    
    private:

        //自己添加的内容
        //Knn_Graph
        Knn_Graph<Distance> graph;
        //paths
        std::vector<path> paths;
        //path number
        size_t paths_num;
        //param k of knn Graph
        size_t k;
        //build graph method
        BUILD_GRAPH_METHOD method;
        //Distance 
        Distance distance;
    private:
        void setDataset(const std::vector<path>& datas)
        {
            this->paths_num=datas.size();
            for(int i=0;i<paths_num;i++)
            {
                paths.push_back(datas[i]);
            }
        }
    public:
        Gnns_Index(){}

        //自己添加的内容
        Gnns_Index(int paths_num,Gnns_Params params=Gnns_Params())
        {
            this->paths_num=paths_num;
            k=params.Graph_k;
            this->method=params.method;
        }
        Gnns_Index(const std::vector<path>& datas,Gnns_Params params=Gnns_Params())
        {
            setDataset(datas);
            k=params.Graph_k;
            this->method=params.method;
        }


        //分类
        void build_index(const std::string& algorithm,const std::string& graph_index_path,const std::string& graph_dist_path,const std::string& graph_start_end_position_path,const std::string& top_level_graph_points_saved_path,
                        const std::string& top_level_graph_indices_saved_path, const std::string& top_level_graph_distances_saved_path, const std::string& top_level_graph_start_end_position_saved_path, bool rebuild=false)
        {
            if(rebuild==true)
            {
                std::cout<<"The saved graph is not exist and building a graph, this may cost lots of time..."<<std::endl;
                graph.build_graph(algorithm,paths,k,method);
                graph.save_graph(graph_index_path,graph_dist_path,graph_start_end_position_path,top_level_graph_points_saved_path,top_level_graph_indices_saved_path,top_level_graph_distances_saved_path,top_level_graph_start_end_position_saved_path);
            }
            else 
            {
                try{
                    std::cout<<"The saved graph is exist, downloading the graph..."<<std::endl;
                    graph.load_graph(graph_index_path,graph_dist_path,graph_start_end_position_path,top_level_graph_points_saved_path,top_level_graph_indices_saved_path,top_level_graph_distances_saved_path,top_level_graph_start_end_position_saved_path);
                }catch(std::exception& e)
                {
                    std::cout<<"download file error, next building a graph"<<std::endl;
                    std::cout<<"The saved graph is not exist and building a graph, this may cost lots of time..."<<std::endl;
                    graph.build_graph(algorithm,paths,k,method);
                    graph.save_graph(graph_index_path,graph_dist_path,graph_start_end_position_path,top_level_graph_points_saved_path,top_level_graph_indices_saved_path,top_level_graph_distances_saved_path,top_level_graph_start_end_position_saved_path);
                }
            }
        }
        //分类
        void knn_search(const std::vector<path>& queryies,
            std::vector<std::vector<int>>& indices,
            std::vector<std::vector<float>>& dists,
            size_t knn,
            const Search_Params& params)
        {
            if(params.E>k)
            {
                std::cout<<"WARNINGS: The search expand E in param exceeds the k (param of the build knn graph)"<<std::endl;
            }
            srand((unsigned)time(0));
            for(int i=0;i<queryies.size();++i)
            {   //std::cout<<"执行到这1"<<std::endl;
                find_neighbors(queryies[i],indices[i],dists[i],knn,params);
            }
        }
        //分类 
    private:
        void find_neighbors(const path& query,
        std::vector<int>& index,
        std::vector<float>& dist,
        size_t knn,
        const Search_Params& params)
        {
            size_t R=params.R;
            std::set<std::pair<float,int>> dist_and_index;

            std::vector<std::pair<std::pair<bool,double>,std::pair<int,int>>> visitpool(paths.size());
            int sum_calculate=0;
            while(R--)
            {
                std::vector<int> top_level_points=graph.get_top_level_points();
                // std::cout<<"top_level_points的索引号 ";
                // for(int i=0;i<top_level_points.size();i++)
                // {
                //     std::cout<<top_level_points[i]<<" ";
                // }
                //std::cout<<std::endl;
                //std::vector<std::pair<int,std::vector<int>>> top_level_indices=graph.get_top_level_indices();
                //std::vector<std::vector<std::pair<int,int>>> top_level_start_end_positions=graph.get_top_level_start_end_positions();

                int top_level_points_num=top_level_points.size();
                size_t top_level_v_it_index=rand()%top_level_points_num;
                size_t top_level_v_it=top_level_points[top_level_v_it_index];
                std::cout<<"顶层图进入点 "<<top_level_v_it<<std::endl;
                subResult r=execute("efficientAlgorithm",query,paths[top_level_v_it]);
                //std::cout<<"执行到这2"<<std::endl;

                //存储的是和query的distances和index;
                std::set<std::pair<double,int>> top_level_dist_and_index;

                //存储的是index、bool、distance、start、end;
                std::vector<std::pair<int,std::pair<std::pair<bool,double>,std::pair<int,int>>>> top_level_visitpool;

                for(const auto& it:top_level_points)
                {
                    top_level_visitpool.push_back(std::make_pair(it,std::make_pair(std::make_pair(false,0),std::make_pair(-1,-1))));
                }
                //std::cout<<"执行到这里1"<<std::endl;
                for(auto& it:top_level_visitpool)
                {
                    if(it.first==top_level_v_it)
                    {
                        it.second.first.first=true;
                        it.second.first.second=r.second;
                        it.second.second.first=r.first.first;
                        it.second.second.second=r.first.second;
                        break;
                    }
                }
                //std::cout<<"执行到这里2"<<std::endl;
                double top_level_min_dist=-1;
                if(r.second<top_level_min_dist||top_level_min_dist==-1)
                {
                    top_level_min_dist=r.second;
                }
                //std::cout<<"执行到这3"<<std::endl;
                while(true)
                {   
                    std::vector<int> top_level_neighbors=graph.get_top_level_neighbors(top_level_v_it);
                    std::vector<std::pair<int,int>> top_level_start_end_positions=graph.get_top_level_start_end_positions(top_level_v_it);
                    //存储distances和(start,end);

                    //std::cout<<"top_level_start_end_position的内容"<<std::endl;
                    //for(int j=0;j<top_level_start_end_positions.size();j++)
                    //{
                    //    std::cout<<"start "<<top_level_start_end_positions[j].first<<" to end "<<top_level_start_end_positions[j].second<<std::endl;
                    //}
                    //std::cout<<std::endl;

                    std::vector<std::pair<double,std::pair<int,int>>> top_level_dist_to_query(top_level_neighbors.size());

                    int top_calculate_neighbor_node_num=0;
                    //std::cout<<"top_level_neighbors.size() "<<top_level_neighbors.size()<<std::endl;
                    //std::cout<<"遍历邻居 ";
                    //for(int j=0;j<top_level_neighbors.size();++j)
                    //{
                    //    std::cout<<top_level_neighbors[j]<<" ";
                    //}
                    //std::cout<<std::endl;
                    for(int i=0;i<top_level_neighbors.size();++i)
                    {
                        int top_level_neighbor_number=top_level_neighbors[i];
                        //std::cout<<top_level_neighbor_number<<" ";
                        for(auto& it:top_level_visitpool)
                        {
                            if(it.first==top_level_neighbor_number)
                            {
                                if(it.second.first.first==true)
                                {
                                    top_level_dist_to_query[i].first=it.second.first.second;
                                    top_level_dist_to_query[i].second.first=it.second.second.first;
                                    top_level_dist_to_query[i].second.second=it.second.second.second;
                                    //top_calculate_neighbor_node_num++;
                                }
                                else 
                                {
                                    subResult r=execute("efficientAlgorithm",query,paths[top_level_neighbor_number]);
                                    top_level_dist_to_query[i].first=r.second;
                                    top_level_dist_to_query[i].second.first=r.first.first;
                                    top_level_dist_to_query[i].second.second=r.first.second;
                                    //top_calculate_neighbor_node_num++;

                                    it.second.first.first=true;
                                    it.second.first.second=r.second;
                                    it.second.second.first=r.first.first;
                                    it.second.second.second=r.first.second;
                                    sum_calculate+=1;
                                }
                                //std::cout<<"执行到这里3"<<std::endl;
                                break;
                            }
                        }
                        top_level_dist_and_index.insert(std::pair<float,int>(top_level_dist_to_query[i].first,top_level_neighbor_number));
                    }
                    //std::cout<<"执行到这里4"<<std::endl;
                    int top_level_min_index=(std::minmax_element(top_level_dist_to_query.begin(),top_level_dist_to_query.end())).first-top_level_dist_to_query.begin();
                    if(top_level_dist_to_query[top_level_min_index].first<top_level_min_dist||top_level_min_dist==-1)
                    {
                        top_level_v_it=top_level_neighbors[top_level_min_index];
                        top_level_min_dist=top_level_dist_to_query[top_level_min_index].first;
                        subResult r=std::make_pair(std::make_pair(top_level_dist_to_query[top_level_min_index].second.first,top_level_dist_to_query[top_level_min_index].second.second),top_level_dist_to_query[top_level_min_index].first);
                        //std::cout<<"执行到这里5"<<std::endl;
                    }
                    else 
                    {
                        break;
                    }
                }
                //std::cout<<"执行到这4"<<std::endl;
                size_t v_it=top_level_v_it;

                //size_t v_it=rand()%30000;
                std::cout<<"底层图进入点v_it "<<v_it<<std::endl;
                double min_dist=-1;
                while(true)
                {
                    std::vector<int> neighbors=graph.get_neighbors(v_it,params.E);


                    //dist_and_index.clear();


                    std::vector<std::pair<int,int>> start_end_positions=graph.get_start_end_positions(v_it,params.E);
                    std::vector<std::pair<double,std::pair<int,int>>> dist_to_query(neighbors.size());

                    // std::cout<<"start_end_position的内容"<<std::endl;
                    // for(int j=0;j<start_end_positions.size();j++)
                    // {
                    //     std::cout<<"start "<<start_end_positions[j].first<<" to end "<<start_end_positions[j].second<<std::endl;
                    // }
                    // std::cout<<std::endl;

                    int calculate_neighbor_node_num=0;
                    //std::cout<<"筛选和计算的邻居 "<<std::endl;
                    //std::cout<<"还能到这1"<<std::endl;
                    for(int i=0;i<neighbors.size();++i)
                    {
                        int neighbor_index=neighbors[i];
                        // if(start_end_positions[i].second==-1&&start_end_positions[i].first==-1)
                        // {
                        //     //std::cout<<"筛掉"<<neighbor_index<<" query similar from "<<r.first.first<<" to "<<r.first.second<<" data similar from "<<start_end_positions[i].first<<" to "<<start_end_positions[i].second<<std::endl;
                        //     continue;
                        // }

                        // else 
                        {   //std::cout<<"还能到这6"<<std::endl;
                            //std::cout<<"neighbor_index "<<neighbor_index<<std::endl;
                            //std::cout<<"计算"<<neighbor_index<<" query similar from "<<r.first.first<<" to "<<r.first.second<<" data similar from "<<start_end_positions[i].first<<" to "<<start_end_positions[i].second<<std::endl;
                            if(visitpool[neighbor_index].first.first==true)
                            {
                                //std::cout<<"还能到这22"<<std::endl;
                                //std::cout<<"假的计算"<<neighbor_index<<" query similar from "<<r.first.first<<" to "<<r.first.second<<" data similar from "<<start_end_positions[i].first<<" to "<<start_end_positions[i].second<<std::endl;
                                dist_to_query[calculate_neighbor_node_num].first=visitpool[neighbor_index].first.second;
                                dist_to_query[calculate_neighbor_node_num].second.first=visitpool[neighbor_index].second.first;
                                dist_to_query[calculate_neighbor_node_num].second.second=visitpool[neighbor_index].second.second;
                                calculate_neighbor_node_num++;
                                //std::cout<<"还能到这2"<<std::endl;
                            }
                            else
                            {
                                //std::cout<<"还能到这33"<<std::endl;
                                //std::cout<<"真的计算"<<neighbor_index<<" query similar from "<<r.first.first<<" to "<<r.first.second<<" data similar from "<<start_end_positions[i].first<<" to "<<start_end_positions[i].second<<std::endl;
                                subResult r=execute("efficientAlgorithm",query,paths[neighbor_index]);

                                //std::cout<<"efficientAlgorithm计算的query与data之间的r query similar from "<<r.first.first<<" to "<<r.first.second<<std::endl;
                                dist_to_query[calculate_neighbor_node_num].first=r.second;
                                dist_to_query[calculate_neighbor_node_num].second.first=r.first.first;
                                dist_to_query[calculate_neighbor_node_num].second.second=r.first.second;
                                calculate_neighbor_node_num++;

                                visitpool[neighbor_index].first.first=true;
                                visitpool[neighbor_index].first.second=r.second;
                                visitpool[neighbor_index].second.first=r.first.first;
                                visitpool[neighbor_index].second.second=r.first.second;
                                sum_calculate+=1;
                                //std::cout<<"还能到这3"<<std::endl;
                            }
                            dist_and_index.insert(std::pair<double,int>(dist_to_query[i].first,neighbor_index));
                            //std::cout<<"还能到这4"<<std::endl;
                        }
                    }
                    //std::cout<<std::endl;
                    //std::cout<<"还能到这5"<<std::endl;
                    int min_index=(std::minmax_element(dist_to_query.begin(),dist_to_query.end())).first-dist_to_query.begin();
                    // std::cout<<"v_it "<<v_it<<" 的邻居 ";
                    // for(int t=0;t<neighbors.size();t++)
                    // {
                    //     std::cout<<neighbors[t]<<" ";
                    // }
                    // std::cout<<std::endl;
                    if(dist_to_query[min_index].first<min_dist||min_dist==-1)
                    {
                        v_it=neighbors[min_index];
                        std::cout<<"下一个底层图搜索点v_it "<<v_it<<std::endl;
                        min_dist=dist_to_query[min_index].first;
                        subResult r=std::make_pair(std::make_pair(dist_to_query[min_index].second.first,dist_to_query[min_index].second.second),dist_to_query[min_index].first);
                    }
                    else 
                    {
                        break;
                    }
                }
                if(R==0&&dist_and_index.size()<knn)
                {
                    std::cout<<"R++"<<std::endl;
                    R+=1;
                }
            }
            size_t k=0;
            //std::cout<<"预测top-k "<<std::endl;
            for(auto it=dist_and_index.begin();it!=dist_and_index.end();++it)
            {
                dist[k]=it->first;
                index[k]=it->second;
                //std::cout<<index[k]<<" ";
                k++;
                if(k==knn)
                {
                    break;
                }
            }
            //std::cout<<std::endl;
            std::cout<<"计算这个query的top-k总共计算了 "<<sum_calculate<<" 条轨迹"<<std::endl;
        }
};


//         void build_index(const std::string & algorithm,const std::string& graph_index_path,const std::string& graph_dist_path,const std::string& top_level_graph_points_path,const std::string& top_level_graph_indices_path,const std::string& top_level_graph_distances_path,bool rebuild=false)
//         {
//             if(rebuild==true)
//             {
//                 std::cout<<"The saved graph is not exist and building a graph, this may cost lots of time..."<<std::endl;
//                 graph.build_graph(algorithm,paths,k,method);
//                 graph.save_graph(graph_index_path,graph_dist_path,top_level_graph_points_path,top_level_graph_indices_path,top_level_graph_distances_path);
//             }
//             else 
//             {
//                 try{
//                     std::cout<<"The saved graph is exist, downloading the graph..."<<std::endl;
//                     graph.load_graph(graph_index_path,graph_dist_path,top_level_graph_points_path,top_level_graph_indices_path,top_level_graph_distances_path);
//                 }catch(std::exception& e )
//                 {
//                     std::cout<<"download file error, next building a graph"<<std::endl;
//                     std::cout<<"The saved graph is not exist and building a graph, this may cost lots of time..."<<std::endl;
//                     graph.build_graph(algorithm,paths,k,method);
//                     graph.save_graph(graph_index_path,graph_dist_path,top_level_graph_points_path,top_level_graph_indices_path,top_level_graph_distances_path);
//                 }
//             }
//         }
//         void knn_search(const std::vector<path>& queryies,
//             std::vector<std::vector<int>>& indices,
//             std::vector<std::vector<float>>& dists,
//             size_t knn,
//             const Search_Params& params)
//         {
//             if(params.E>k)
//             {
//                 std::cout<<"WARNINGS: The search expand E in param exceeds the k (param of the build knn graph)"<<std::endl;
//             }
//             srand((unsigned)time(0));

//             for(int i=0;i<queryies.size();++i)
//             {
//                 //std::cout<<"第 "<<i<<" 条query"<<endl;
//                 find_neighbors(queryies[i],indices[i],dists[i],knn,params);
//             }
//         }
//     private:
//         void find_neighbors(const path& query,
//         std::vector<int>& index,
//         std::vector<float>& dist,
//         size_t knn,
//         const Search_Params& params)

//         //2.过滤+top_selectedgrid过滤网格顶层图搜索**************************************************************************
//                 {
//             size_t R=params.R;
//             std::set<std::pair<float,int>>dist_and_index;

//             std::vector<std::pair<bool,double>> visitpool(paths.size());
//             int sum_calculate=0;
//             while(R--)
//             {   

//                 //top-level-graph上层graph所用数据结构
//                 std::vector<int> top_level_points=graph.get_top_level_points();
//                 //std::cout<<"top_level_points的size "<<top_level_points.size()<<" "<<"top_level_points集合"<<std::endl;
//                 //for(int i=0;i<top_level_points.size();i++)
//                 //{
//                 //    std::cout<<top_level_points[i]<<" ";
//                 //}
//                 //std::cout<<std::endl;

//                 std::vector<std::pair<int,std::vector<int>>> top_level_indices=graph.get_top_level_indices();
//                 // std::cout<<"top_level_indices的size "<<top_level_indices.size()<<" "<<"top_level_points集合"<<std::endl;
//                 // for(const auto& pair_vector:top_level_indices)
//                 // {
//                 //     std::cout<<"top_level_point "<<pair_vector.first<<" 它的邻居节点 ";
//                 //     for(int value_vector:pair_vector.second)
//                 //     {
//                 //         std::cout<<value_vector<<" ";
//                 //     }
//                 //     std::cout<<std::endl;
                    
//                 // }
//                 int top_level_points_num=top_level_points.size();
//                 size_t top_level_v_it_index=rand()%top_level_points_num;
//                 size_t top_level_v_it=top_level_points[top_level_v_it_index];
//                 std::cout<<"顶层图进入点 "<<top_level_v_it<<std::endl;

//                 std::set<std::pair<double,int>> top_level_dist_and_index;
//                 std::vector<std::pair<int,std::pair<bool,double>>> top_level_visitpool;
//                 //std::cout<<"执行到此5"<<std::endl;
//                 for(const auto& it:top_level_points)
//                 {
//                     top_level_visitpool.push_back(std::make_pair(it,std::make_pair(false,0)));
//                     //std::cout<<"执行for 循环"<<std::endl;
//                 }
//                 //std::cout<<"top_level_visitpool "<<top_level_visitpool.size()<<std::endl;

//                 double top_level_min_dist=-1;
//                 while(true)
//                 {
//                     //存储邻居的索引号和邻居到query的距离
//                     std::vector<int> top_level_neighbors=graph.get_top_level_neighbors(top_level_v_it);
//                     //std::cout<<"top_level_v_it "<<top_level_v_it<<std::endl;
//                     std::vector<double> top_level_dist_to_query(top_level_neighbors.size());
//                     //std::cout<<"top_level_neighbors的size "<<top_level_neighbors.size()<<std::endl;
//                     //std::cout<<"top_level_v_it的邻居 ";
//                     //for(int i=0;i<top_level_neighbors.size();i++)
//                     //{
//                     //    std::cout<<top_level_neighbors[i]<<" ";
//                     //}
//                     //std::cout<<std::endl;

//                     for(int i=0;i<top_level_neighbors.size();++i)
//                     {
//                         //top_level_neighbor_number表示节点在全体轨迹节点中的序号
//                         int top_level_neighbor_number=top_level_neighbors[i];
//                         for(auto& it:top_level_visitpool)
//                         {
//                             if(it.first==top_level_neighbor_number)
//                             {
//                                 if(it.second.first==true)
//                                 {
//                                     top_level_dist_to_query[i]=it.second.second;
//                                     //std::cout<<"进入if "<<std::endl;
//                                 }
//                                 else 
//                                 {
//                                     subResult r=execute("efficientAlgorithm",query,paths[top_level_neighbor_number]);
//                                     top_level_dist_to_query[i]=r.second;
//                                     it.second.first=true;
//                                     it.second.second=r.second;
//                                     //std::cout<<"进入else "<<std::endl;
//                                     sum_calculate+=1;
//                                 }
//                                 break;
//                             }
//                             //std::cout<<"不是这个节点"<<std::endl;   
//                         }
//                         //std::cout<<"执行这个for 循环"<<std::endl;
//                         top_level_dist_and_index.insert(std::pair<float,int>(top_level_dist_to_query[i],top_level_neighbor_number));
//                         //std::cout<<"插入成功"<<std::endl;
//                     }
//                     //std::cout<<"for 循环执行结束"<<std::endl;

//                     //std::cout<<"top_level_dist_to_query "<<std::endl;
//                     //for(int i=0;i<top_level_dist_to_query.size();i++)
//                     //{
//                     //    std::cout<<top_level_dist_to_query[i]<<" ";
//                     //}
//                     //std::cout<<std::endl;

//                     int top_level_min_index=(std::minmax_element(top_level_dist_to_query.begin(),top_level_dist_to_query.end())).first-top_level_dist_to_query.begin();
//                     //std::cout<<"top_level_min_index "<<top_level_min_index<<std::endl;

//                     if(top_level_dist_to_query[top_level_min_index]<top_level_min_dist||top_level_min_dist==-1)
//                     {
//                         //std::cout<<"到达if 换下一个top_level_v_it"<<std::endl;
                        
//                         //for(int i=0;i<top_level_neighbors.size();i++)
//                         //{
//                         //    std::cout<<top_level_neighbors[i]<<" ";
//                         //}
//                         //std::cout<<std::endl;

//                         top_level_v_it=top_level_neighbors[top_level_min_index];
//                         top_level_min_dist=top_level_dist_to_query[top_level_min_index];
                        
//                     }
//                     else 
//                     {
//                         break;
//                     }
//                 }
//                 //std::cout<<"执行到此6"<<std::endl;
//                 size_t v_it=top_level_v_it;
//                 //size_t v_it=rand()%paths_num;
//                 std::cout<<"底层图进入点v_it "<<v_it<<std::endl;
//                 double min_dist=-1;


//                 while(true)
//                 {
//                     std::vector<int> neighbors=graph.get_neighbors(v_it,params.E);


//                     //dist_and_index的清除,只存储当前搜索点的邻居内容;目前做的唯一修改
//                     //dist_and_index.clear();



//                     std::vector<double> dist_to_query(neighbors.size());
//                     for(int i=0;i<neighbors.size();++i)
//                     {
//                         int neighbor_index=neighbors[i];
//                         if(visitpool[neighbor_index].first==true)
//                         {
//                             dist_to_query[i]=visitpool[neighbor_index].second;
//                         }
//                         else 
//                         {
//                             subResult r=execute("efficientAlgorithm",query,paths[neighbor_index]);
//                             dist_to_query[i]=r.second;
//                             visitpool[neighbor_index].first=true;
//                             visitpool[neighbor_index].second=r.second;
//                             sum_calculate+=1;
//                         }
//                         dist_and_index.insert(std::pair<double,int>(dist_to_query[i],neighbor_index));
//                     }
//                     int min_index=(std::minmax_element(dist_to_query.begin(),dist_to_query.end())).first-dist_to_query.begin();
//                     std::cout<<"v_it "<<v_it<<" 的邻居 ";
//                     for(int t=0;t<neighbors.size();t++)
//                     {
//                         std::cout<<neighbors[t]<<" ";
//                     }
//                     std::cout<<std::endl;
//                     if(dist_to_query[min_index]<min_dist||min_dist==-1)
//                     {
//                         // std::cout<<"v_it "<<v_it<<" 的邻居 ";
//                         // for(int t=0;t<neighbors.size();t++)
//                         // {
//                         //     std::cout<<neighbors[t]<<" ";
//                         // }
//                         // std::cout<<std::endl;
//                         // std::cout<<"dist_and_index中的内容"<<std::endl;
//                         // for(const auto& element:dist_and_index)
//                         // {
//                         //     std::cout<<"from the neighbor to query's Distance: "<<element.first<<", the neighbor's Index: "<<element.second<<std::endl;
//                         // }
//                         v_it=neighbors[min_index]; 
//                         std::cout<<"下一个底层图搜索点v_it "<<v_it<<std::endl;
//                         min_dist=dist_to_query[min_index];
//                     }
//                     else 
//                     {   

//                         // //计算一下反向表
//                         // if()
//                         // {}
//                         // else 
//                         // {
//                         //     break;
//                         // }

//                         // std::cout<<"底层图最优点"<<v_it<<std::endl;
//                         // std::cout<<"底层最优点 "<<v_it<<" 的邻居 ";
//                         // for(int t=0;t<neighbors.size();t++)
//                         // {
//                         //     std::cout<<neighbors[t]<<" ";
//                         // }
//                         // std::cout<<std::endl;
//                         // std::cout<<"dist_and_index中的内容"<<std::endl;
//                         // for(const auto& element:dist_and_index)
//                         // {
//                         //     std::cout<<"from the neighbor to query's Distance: "<<element.first<<", the neighbor's Index: "<<element.second<<std::endl;
//                         // }
//                         break;
//                     }
//                 }

                
//                 if(R==0&&dist_and_index.size()<knn)
//                 {
//                     std::cout<<"R++"<<std::endl;
//                     R+=1;
//                 }
//                 //std::cout<<"执行到此7"<<std::endl;
//             }
//             size_t k=0;
//             std::cout<<"预测top-k "<<std::endl;
//             for(auto it=dist_and_index.begin();it!=dist_and_index.end();++it)
//             {
//                 dist[k]=it->first;
//                 index[k]=it->second;
//                 std::cout<<index[k]<<" ";
//                 k++;
//                 if(k==knn)
//                     break;
//             }
//             std::cout<<std::endl;
//             std::cout<<"计算这个query的top-k总共计算了 "<<sum_calculate<<" 条轨迹"<<std::endl;
//         }
//     };


















// //1.top_dropnearester_selectedrandom最开始的顶层图搜索***********************************************************************
// ////        {
// ////            size_t R=params.R;
// ////            std::set<std::pair<float,int>>dist_and_index;
// ////            std::vector<std::pair<bool,double>> visitpool(paths.size());
// ////            while(R--)
// ////            {   

//                 //top-level-graph上层graph所用数据结构
// ////                std::vector<int> top_level_points=graph.get_top_level_points();
//                 //std::cout<<"top_level_points的size "<<top_level_points.size()<<" "<<"top_level_points集合"<<std::endl;
//                 //for(int i=0;i<top_level_points.size();i++)
//                 //{
//                 //    std::cout<<top_level_points[i]<<" ";
//                 //}
//                 //std::cout<<std::endl;

// ////                std::vector<std::pair<int,std::vector<int>>> top_level_indices=graph.get_top_level_indices();
//                 // std::cout<<"top_level_indices的size "<<top_level_indices.size()<<" "<<"top_level_points集合"<<std::endl;
//                 // for(const auto& pair_vector:top_level_indices)
//                 // {
//                 //     std::cout<<"top_level_point "<<pair_vector.first<<" 它的邻居节点 ";
//                 //     for(int value_vector:pair_vector.second)
//                 //     {
//                 //         std::cout<<value_vector<<" ";
//                 //     }
//                 //     std::cout<<std::endl;
                    
//                 // }
// ////                int top_level_points_num=top_level_points.size();
// ////                size_t top_level_v_it_index=rand()%top_level_points_num;
// ////                size_t top_level_v_it=top_level_points[top_level_v_it_index];

// ////                std::set<std::pair<float,int>> top_level_dist_and_index;
// ////                std::vector<std::pair<int,std::pair<bool,double>>> top_level_visitpool;
// ////                std::cout<<"执行到此5"<<std::endl;
// ////                for(const auto& it:top_level_points)
// ////                {
// ////                    top_level_visitpool.push_back(std::make_pair(it,std::make_pair(false,0)));
//                     //std::cout<<"执行for 循环"<<std::endl;
// ////                }
//                 //std::cout<<"top_level_visitpool "<<top_level_visitpool.size()<<std::endl;

// ////                float top_level_min_dist=-1;
// ////                while(true)
// ////                {
//                     //存储邻居的索引号和邻居到query的距离
// ////                    std::vector<int> top_level_neighbors=graph.get_top_level_neighbors(top_level_v_it);
//                     //std::cout<<"top_level_v_it "<<top_level_v_it<<std::endl;
// ////                    std::vector<float> top_level_dist_to_query(top_level_neighbors.size());
//                     //std::cout<<"top_level_neighbors的size "<<top_level_neighbors.size()<<std::endl;
//                     //std::cout<<"top_level_v_it的邻居 ";
//                     //for(int i=0;i<top_level_neighbors.size();i++)
//                     //{
//                     //    std::cout<<top_level_neighbors[i]<<" ";
//                     //}
//                     //std::cout<<std::endl;

// ////                    for(int i=0;i<top_level_neighbors.size();++i)
// ////                    {
//                         //top_level_neighbor_number表示节点在全体轨迹节点中的序号
// ////                        int top_level_neighbor_number=top_level_neighbors[i];
// ////                        for(auto& it:top_level_visitpool)
// ////                        {
// ////                            if(it.first==top_level_neighbor_number)
// ////                            {
// ////                                if(it.second.first==true)
// ////                                {
// ////                                    top_level_dist_to_query[i]=it.second.second;
//                                     //std::cout<<"进入if "<<std::endl;
// ////                                }
// ////                                else 
// ////                                {
// ////                                    subResult r=execute("efficientAlgorithm",query,paths[top_level_neighbor_number]);
// ////                                    top_level_dist_to_query[i]=r.second;
// ////                                    it.second.first=true;
// ////                                    it.second.second=r.second;
//                                     //std::cout<<"进入else "<<std::endl;
// ////                                }
// ////                                break;
// ////                            }
//                             //std::cout<<"不是这个节点"<<std::endl;   
// ////                        }
//                         //std::cout<<"执行这个for 循环"<<std::endl;
// ////                        top_level_dist_and_index.insert(std::pair<float,int>(top_level_dist_to_query[i],top_level_neighbor_number));
//                         //std::cout<<"插入成功"<<std::endl;
// ////                    }
//                     //std::cout<<"for 循环执行结束"<<std::endl;

//                     //std::cout<<"top_level_dist_to_query "<<std::endl;
//                     //for(int i=0;i<top_level_dist_to_query.size();i++)
//                     //{
//                     //    std::cout<<top_level_dist_to_query[i]<<" ";
//                     //}
// ////                    std::cout<<std::endl;

// ////                    int top_level_min_index=(std::minmax_element(top_level_dist_to_query.begin(),top_level_dist_to_query.end())).first-top_level_dist_to_query.begin();
//                     //std::cout<<"top_level_min_index "<<top_level_min_index<<std::endl;

// ////                    if(top_level_dist_to_query[top_level_min_index]<top_level_min_dist||top_level_min_dist==-1)
// ////                    {
//                         //std::cout<<"到达if 换下一个top_level_v_it"<<std::endl;
                        
//                         //for(int i=0;i<top_level_neighbors.size();i++)
//                         //{
//                         //    std::cout<<top_level_neighbors[i]<<" ";
//                         //}
//                         //std::cout<<std::endl;

// ////                        top_level_v_it=top_level_neighbors[top_level_min_index];
// ////                        top_level_min_dist=top_level_dist_to_query[top_level_min_index];
                        
// ////                    }
// ////                    else 
// ////                    {
// ////                        break;
// ////                    }
// ////                }
// ////                std::cout<<"执行到此6"<<std::endl;
// ////                size_t v_it=top_level_v_it;
// ////                float min_dist=-1;
// ////                while(true)
// ////                {
// ////                    std::vector<int> neighbors=graph.get_neighbors(v_it,params.E);
// ////                    std::vector<float> dist_to_query(neighbors.size());
// ////                    for(int i=0;i<neighbors.size();++i)
// ////                    {
// ////                        int neighbor_index=neighbors[i];
// ////                        if(visitpool[neighbor_index].first==true)
// ////                        {
// ////                            dist_to_query[i]=visitpool[neighbor_index].second;
// ////                        }
// ////                        else 
// ////                        {
// ////                            subResult r=execute("efficientAlgorithm",query,paths[neighbor_index]);
// ////                            dist_to_query[i]=r.second;
// ////                            visitpool[neighbor_index].first=true;
// ////                            visitpool[neighbor_index].second=r.second;
// ////                        }
// ////                        dist_and_index.insert(std::pair<float,int>(dist_to_query[i],neighbor_index));
// ////                    }
// ////                    int min_index=(std::minmax_element(dist_to_query.begin(),dist_to_query.end())).first-dist_to_query.begin();
// ////                    if(dist_to_query[min_index]<min_dist||min_dist==-1)
// ////                    {
// ////                        v_it=neighbors[min_index];
// ////                        min_dist=dist_to_query[min_index];
// ////                    }
// ////                    else 
// ////                    {
// ////                        break;
// ////                    }
// ////                }
// ////                if(R==0&&dist_and_index.size()<knn)
// ////                {
// ////                    std::cout<<"R++"<<std::endl;
// ////                    R+=1;
// ////                }
//                 //std::cout<<"执行到此7"<<std::endl;
// ////            }
// ////            size_t k=0;
// ////            for(auto it=dist_and_index.begin();it!=dist_and_index.end();++it)
// ////            {
// ////                dist[k]=it->first;
// ////                index[k]=it->second;
// ////                k++;
// ////                if(k==knn)
// ////                    break;
// ////            }
// ////        }
// ////    };



// //0.最开始的无顶层图搜索********************************************************

//         //用四个注释符注释
//         ////{
//         ////    size_t R=params.R;
//         ////    std::set<std::pair<float,int>>dist_and_index;
//         ////    std::vector<std::pair<bool,double>> visitpool(paths.size());
//         ////    //std::cout<<"搜索路径起始点";
//         ////    while(R--)
//         ////    {
//         ////        size_t v_it=rand()%paths_num;
//         ////        //std::cout<<v_it;
//         ////        
//         ////        //size_t v_it;
//         ////        //查看query与data的距离;
//         ////        
//         ////        //while(true)
//         ////        //{
//         ////        //    int v_it1,v_it2,v_it3;
//         ////        //    double distance1,distance2,distance3;
//         ////        //    v_it1=rand()%paths_num;
//         ////        //    v_it2=rand()%paths_num;
//         ////        //    v_it3=rand()%paths_num;
//         ////
//         ////        //    subResult r1=execute("efficientAlgorithm",query,paths[v_it1]);
//         ////        //    distance1=r1.second;
//         ////        //    subResult r2=execute("efficientAlgorithm",query,paths[v_it2]);
//         ////        //    distance2=r2.second;
//         ////        //    subResult r3=execute("efficientAlgorithm",query,paths[v_it3]);
//         ////        //    distance3=r3.second;
//         ////        //    double min_distance=std::min({distance1,distance2,distance3});
//         ////        //    double max_distance=std::max({distance1,distance2,distance3});
//         ////
//         ////        //    if(min_distance<55)
//         ////        //    {
//         ////        //        if(distance1==min_distance)
//         ////        //            v_it=v_it1;
//         ////        //        else if(distance2=min_distance)
//         ////        //            v_it=v_it2;
//         ////        //        else 
//         ////        //            v_it=v_it3;
//         ////        //        break;
//         ////        //    }
//         ////        //    else 
//         ////        //    {
//         ////                //std::vector<int>neighbors1=graph.get_neighbors(v_it1,100);
//         ////                //std::vector<int>neighbors2=graph.get_neighbors(v_it2,100);
//         ////                //std::vector<int>neighbors3=graph.get_neighbors(v_it3,100);
//         ////                //std::vector<int>complement;
//         ////                //for(int i=1;i<3000;++i)
//         ////                //{
//         ////                //    if(i != v_it1 && i != v_it2 && i != v_it3 && std::find(neighbors1.begin(), neighbors1.end(), i) == neighbors1.end() &&
//         ////                //        std::find(neighbors2.begin(), neighbors2.end(), i) == neighbors2.end() &&
//         ////                //        std::find(neighbors3.begin(), neighbors3.end(), i) == neighbors3.end())
//         ////                //    {
//         ////                //        complement.push_back(i);
//         ////                //    }
//         ////                //}
//         ////                //size_t complement_index=rand()%(complement.size());
//         ////                //v_it=complement[complement_index];
//         ////                //break;
//         ////
//         ////        //        int max_distance_v_it;
//         ////        //        if(distance1==max_distance)
//         ////        //            max_distance_v_it=v_it1;
//         ////        //        else if(distance2=max_distance)
//         ////        //            max_distance_v_it=v_it2;
//         ////        //        else 
//         ////        //            max_distance_v_it=v_it3;
//         ////        //        std::vector<int>one_hop_neigbhors=graph.get_neighbors(max_distance_v_it,100);
//         ////        //        std::unordered_set<int>two_hop_neighbors;
//         ////        //        for(int neighbor:one_hop_neigbhors)
//         ////        //        {
//         ////        //            std::vector<int> neighbors_of_neighbor=graph.get_neighbors(neighbor,100);
//         ////        //            for(int n:neighbors_of_neighbor)
//         ////        //            {
//         ////        //                two_hop_neighbors.insert(n);
//         ////        //            }
//         ////        //        }
//         ////        //        std::vector<int> complement;
//         ////        //        for (int i = 1; i < 3000; ++i) 
//         ////        //        {
//         ////        //            if (i != v_it1 && i != v_it2 && i != v_it3 && two_hop_neighbors.find(i) == two_hop_neighbors.end()) 
//         ////        //            {
//         ////        //                complement.push_back(i);
//         ////        //            }
//         ////        //        }
//         ////        //        size_t complement_index=rand()%(complement.size());
//         ////        //        v_it=complement[complement_index];
//         ////        //        break;
//         ////
//         ////
//         ////        //    }
//         ////        //}
//         ////        float min_dist=-1;
//         ////        
//         ////        
//         ////        while(true)
//         ////        {
//         ////            //记录v_it最近的params.E个邻居的index
//         ////            std::vector<int> neighbors=graph.get_neighbors(v_it,params.E);
//         ////            std::vector<float> dist_to_query(neighbors.size());
//         ////            for(int i=0;i<neighbors.size();++i)
//         ////            {
//         ////                int neighbor_index=neighbors[i];
//         ////                
//         ////                if(visitpool[neighbor_index].first==true)
//         ////                {
//         ////                    dist_to_query[i]=visitpool[neighbor_index].second;
//         ////                }
//         ////                else 
//         ////                {
//         ////                    //自己添加的内容
//         ////                    //1.第一种
//         ////                    //qurey-data使用efficient
//         ////                    subResult r=execute("efficientAlgorithm",query,paths[neighbor_index]);
//         ////                
//         ////                    //2.第二种
//         ////                    //query-data使用exactS
//         ////                    //std::cout<<"query的size "<<query.size()<<endl;
//         ////                    //std::cout<<"neighbor的size "<<paths[neighbor_index].size()<<endl;
//         ////
//         ////                    //subResult r=execute("exactS",query,paths[neighbor_index]);
//         ////
//         ////                    dist_to_query[i]=r.second;
//         ////                    //std::cout<<"query与data的distance ";
//         ////                    //std::cout<<r.second<<endl;
//         ////
//         ////
//         ////                    //3.1第三种原始版
//         ////                    //query-data原始
//         ////                    //double dist=minSubTrajectory1query(query,paths[neighbor_index]);
//         ////
//         ////                    //3.2第三种改进版
//         ////                    //query-data改进
//         ////                    //double dist=minSubTrajectory2query(query,paths[neighbor_index]);
//         ////
//         ////                    //dist_to_query[i]=dist;
//         ////
//         ////                 
//         ////                    //dist_to_query[i]=distance(paths[neighbor_index],query);
//         ////
//         ////                    visitpool[neighbor_index].first=true;
//         ////                    visitpool[neighbor_index].second=r.second;
//         ////                }
//         ////                dist_and_index.insert(std::pair<float,int>(dist_to_query[i],neighbor_index));
//         ////            }
//         ////
//         ////            int min_index=(std::minmax_element(dist_to_query.begin(),dist_to_query.end())).first-dist_to_query.begin();
//         ////
//         ////
//         ////            if(dist_to_query[min_index]<min_dist||min_dist==-1)
//         ////            {
//         ////                v_it=neighbors[min_index];
//         ////                //std::cout<<" "<<v_it;
//         ////                min_dist=dist_to_query[min_index];
//         ////            }
//         ////            else 
//         ////            {
//         ////                break;
//         ////            }
//         ////        }
//         ////        if(R==0&&dist_and_index.size()<knn)
//         ////        {
//         ////            std::cout<<"R++"<<std::endl;
//         ////            R+=1;
//         ////        }
//         ////    }
//         ////    size_t k=0;
//         ////    //std::cout<<endl;
//         ////    //std::cout<<" 最终最近点的top-k邻居";
//         ////    for(auto it=dist_and_index.begin();it!=dist_and_index.end();++it)
//         ////    {
//         ////        dist[k]=it->first;
//         ////        index[k]=it->second;
//         ////        //std::cout<<" "<<index[k];
//         ////        k++;
//         ////        if(k==knn)
//         ////            break;
//         ////    }
//         ////    //std::cout<<endl;
//         ////}
//     //};
// //}



#endif //GNNS_GNNS_INDEX_H
