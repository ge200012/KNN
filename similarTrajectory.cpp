#include "similarTrajectory.h"
#include <cmath>
//自己添加的内容

//namespace gnns
//{


//自己添加的内容
int edr1(rtree_point p1, rtree_point p2, double e)
 {
    double distance = sqrt((p1.get<0>() - p2.get<0>()) * (p1.get<0>() - p2.get<0>()) + (p1.get<1>() - p2.get<1>()) * (p1.get<1>() - p2.get<1>()));
    if (e < distance) return 0;
    return 1;
}
//data-data原始

double minSubTrajectory1(const string & algorithm, const path& path1, const path& path2)
{
    std::vector<std::vector<int>> array_2d(path1.size(),std::vector<int>(path2.size()));
    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            if(edr1(path1[i],path2[j],0.0005)==1)
            {
                array_2d[i][j]=1;
            }
            else 
            {
                array_2d[i][j]=-1;
            }
        }
    }

    //输出创建的数组
    //for(int i=0;i<2;i++)
    //{
    //    for(int j=0;j<path2.size();j++)
    //    {
    //        std::cout<<array_2d[i][j];
    //    }
    //    std::cout<<endl;
    //}
    
    bool breakFlag=false;

    int SimilarScore=-10000;
    int start1=-1;
    int start2=-1;
    int end1=-1;
    int end2=-1;

    bool NoSimilar=true;

    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            
            if(array_2d[i][j]==-1)
            {
                continue;
            }
            NoSimilar=false;
            //std::cout<<"i "<<"j "<<i+1<<" "<<j+1<<endl;
            int total=2*array_2d[i][j];
            int subtrajectotal=total;
            int start_1=i;
            int start_2=j;
            int end_1=i;
            int end_2=j;

            bool breakFlag=false;
            for(int i1=i;i1<path1.size()-1;)
            {
                int j1=j;
                //for(int j1=j;j1<path2.size()-1;)
                for(;j1<path2.size()-1;)
                {
                    //判断朝哪个方向
                    //std::cout<<"执行到这里0"<<endl;
                    if(array_2d[i1+1][j1+1]==1)
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        //std::cout<<"执行到这里1"<<endl;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]==-1)
                    {
                        if(array_2d[i1][j1+1]>=array_2d[i1+1][j1])
                        {
                            total+=array_2d[i1][j1+1];
                            j1+=1;
                            //std::cout<<"执行到这里2"<<endl;
                            if(total>subtrajectotal)
                            {
                                subtrajectotal=total;
                                end_1=i1;
                                end_2=j1;
                            }
                        }
                        else 
                        {
                            total+=array_2d[i1+1][j1];
                            i1+=1;
                            //std::cout<<"执行到这里3"<<endl;
                            if(total>subtrajectotal)
                            {
                                subtrajectotal=total;
                                end_1=i1;
                                end_2=j1;
                            }
                        }
                    }
                    else 
                    {
                        //std::cout<<"执行到这里4"<<endl;
                        total+=-1;
                        j1+=1;
                    }

                    //std::cout<<"total "<<total<<endl;
                    //std::cout<<"i1 "<<"j1 "<<i1+1<<" "<<j1+1<<endl;
                    if(total<=0||i1+1>=path1.size()||j1+1>=path2.size())
                    {
                        breakFlag=true;
                        break;
                    }
                }
                if(j1>=path2.size()-1)
                {
                    breakFlag=true;
                    break;
                }
                if(breakFlag==true)
                {
                    break;
                }
            }
            //std::cout<<"此次子轨迹得分 "<<subtrajectotal<<endl;
            if(subtrajectotal>SimilarScore)
            {
                SimilarScore=subtrajectotal;
                start1=start_1;
                start2=start_2;
                end1=end_1;
                end2=end_2;
                //std::cout<<"第一个子轨迹从 "<<start1+1<<" 到 "<<end1+1<<"第二个子轨迹从 "<<start2+1<<" 到 "<<end2+1<<" 代表性相似子轨迹得分是 "<<SimilarScore<<endl;
            }
            //std:cout<<"subtrajectotal "<<subtrajectotal<<endl;
        }
    }

    if(NoSimilar==true)
        return 0;
    else
        return -SimilarScore;
}


//data-data改进

Result minSubTrajectory2(const string & algorithm, const path& path1, const path& path2)
{
    std::vector<std::vector<int>> array_2d(path1.size(),std::vector<int>(path2.size()));
    int array_sum=0;
    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            if(edr1(path1[i],path2[j],0.0005)==1)
            {
                array_2d[i][j]=1;
                array_sum+=1;
            }
            else 
            {
                array_2d[i][j]=-1;
            }
        }
    }
    //if(path1.size()==126 and path2.size()==111)
    //{
    //    std::cout<<"{";
    //    for(int i=0;i<path1.size();i++)
    //    {
    //        std::cout<<"{";
    //        for(int j=0;j<path2.size();j++)
    //        {
    //            if(j<path2.size()-1)
    //                std::cout<<array_2d[i][j]<<",";
    //           else 
    //                std::cout<<array_2d[i][j];
    //        }
    //        std::cout<<"},"<<endl;
    //    }
    //    std::cout<<"};"<<endl;
    //}

    //输出创建的数组
    //for(int i=0;i<2;i++)
    //{
    //    for(int j=0;j<path2.size();j++)
    //    {
    //        std::cout<<array_2d[i][j];
    //    }
    //    std::cout<<endl;
    //}
    
    bool breakFlag=false;

    int SimilarScore=-10000;
    int start1=-1;
    int start2=-1;
    int end1=-1;
    int end2=-1;

    bool NoSimilar=true;


    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            
            if(array_2d[i][j]==-1)
            {
                continue;
            }
            NoSimilar=false;
            //std::cout<<"i "<<"j "<<i+1<<" "<<j+1<<endl;
            int total=2*array_2d[i][j];
            int subtrajectotal=total;
            int start_1=i;
            int start_2=j;
            int end_1=i;
            int end_2=j;

            bool breakFlag=false;
            for(int i1=i;i1<path1.size()-1;)
            {
                int j1=j;
                //for(int j1=j;j1<path2.size()-1;)
                for(;j1<path2.size()-1;)
                {
                    //判断朝哪个方向
                    if(array_2d[i1][j1]==1&&array_2d[i1+1][j1+1]==1)
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]==1&&array_2d[i1][j1+1]>=array_2d[i1+1][j1])
                    {
                        total-=1;
                        j1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]==1&&array_2d[i1+1][j1]>array_2d[i1][j1+1])
                    {
                        total-=1;
                        i1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]==-1&&array_2d[i1][j1+1]==1)
                    {
                        total+=3*array_2d[i1][j1+1];
                        j1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]==-1&&array_2d[i1+1][j1]==1)
                    {
                        total+=3*array_2d[i1+1][j1];
                        i1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else if(array_2d[i1][j1]!=array_2d[i1+1][j1+1])
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        //std::cout<<"run this code"<<endl;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    else
                    {
                        total+=array_2d[i1][j1+1];
                        j1+=1;
                        if(total>subtrajectotal)
                        {
                            subtrajectotal=total;
                            end_1=i1;
                            end_2=j1;
                        }
                    }
                    //std::cout<<"total "<<total<<endl;
                    //std::cout<<"i1 "<<"j1 "<<i1+1<<" "<<j1+1<<endl;
                    if(total<=0||i1+1>=path1.size()||j1+1>=path2.size())
                    {
                        breakFlag=true;
                        break;
                    }
                }
                if(j1>=path2.size()-1)
                {
                    breakFlag=true;
                    break;
                }
                if(breakFlag==true)
                {
                    break;
                }        
            }
            //std::cout<<"此次子轨迹得分 "<<subtrajectotal<<endl;
            if(subtrajectotal>SimilarScore)
            {
                SimilarScore=subtrajectotal;
                start1=start_1;
                start2=start_2;
                end1=end_1;
                end2=end_2;
                //std::cout<<"第一个子轨迹从 "<<start1+1<<" 到 "<<end1+1<<"第二个子轨迹从 "<<start2+1<<" 到 "<<end2+1<<" 代表性相似子轨迹得分是 "<<SimilarScore<<endl;
            }
            //std:cout<<"subtrajectotal "<<subtrajectotal<<endl;
        }
    }
    if(NoSimilar==true)
    {
        double path=0;
        //std::cout<<"相似得分为: "<<path<<endl;
        std::pair<int,int> first_start_end=std::make_pair(start1,end1);
        std::pair<int,int> second_start_end=std::make_pair(start2,end2);
        Result result=std::make_pair(path,std::make_pair(first_start_end,second_start_end));
        return result;
    }
    else
    {
        double path=-SimilarScore;
        //std::cout<<"相似得分为: "<<(-path)<<endl;
        std::pair<int,int> first_start_end=std::make_pair(start1,end1);
        std::pair<int,int> second_start_end=std::make_pair(start2,end2);
        Result result=std::make_pair(path,std::make_pair(first_start_end,second_start_end));
        return result;
    }
}




//query-data原始

double minSubTrajectory1query(const path& path1, const path& path2)
{
    std::vector<std::vector<int>> array_2d(path1.size(),std::vector<int>(path2.size()));
    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            if(edr1(path1[i],path2[j],0.0005)==1)
            {
                array_2d[i][j]=1;
            }
            else 
            {
                array_2d[i][j]=-1;
            }
        }
    }

    //输出创建的数组
    //for(int i=0;i<2;i++)
    //{
    //    for(int j=0;j<path2.size();j++)
    //    {
    //        std::cout<<array_2d[i][j];
    //    }
    //    std::cout<<endl;
    //}
    
    bool breakFlag=false;

    int SimilarScore=-10000;
    int start1=-1;
    int start2=-1;
    int end1=-1;
    int end2=-1;

    bool NoSimilar=true;

    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            NoSimilar=false;
            //std::cout<<"i "<<"j "<<i+1<<" "<<j+1<<endl;
            int total=2*array_2d[i][j];
            int subtrajectotal=total;
            int start_1=i;
            int start_2=j;
            int end_1=i;
            int end_2=j;

            int breakFlag=0;
            for(int i1=i;i1<path1.size()-1;)
            {
                int j1=j;
                //for(int j1=j;j1<path2.size()-1;)
                for(;j1<path2.size()-1;)
                {
                    //判断朝哪个方向
                    //std::cout<<"执行到这里0"<<endl;
                    if(array_2d[i1+1][j1+1]==1)
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        //std::cout<<"执行到这里1"<<endl;
                        subtrajectotal=total;
                        end_2=j1;
                    }
                    else if(array_2d[i1][j1]==-1)
                    {
                        if(array_2d[i1][j1+1]>=array_2d[i1+1][j1])
                        {
                            total+=array_2d[i1][j1+1];
                            j1+=1;
                            //std::cout<<"执行到这里2"<<endl;
                            subtrajectotal=total;
                            end_2=j1;
                        }
                        else 
                        {
                            total+=array_2d[i1+1][j1];
                            i1+=1;
                            //std::cout<<"执行到这里3"<<endl;
                            subtrajectotal=total;
                            end_2=j1;
                        }
                    }
                    else 
                    {
                        //std::cout<<"执行到这里4"<<endl;
                        total+=-1;
                        j1+=1;
                    }

                    //std::cout<<"total "<<total<<endl;
                    //std::cout<<"i1 "<<"j1 "<<i1+1<<" "<<j1+1<<endl;
                    if(i1+1>=path1.size())
                    {
                        breakFlag=1;
                        break;
                    }
                    if(j1+1>=path2.size())
                    {
                        breakFlag=2;
                        break;
                    }
                }
                if(j1+1>=path2.size())
                {
                    breakFlag=2;
                    break;
                }
                if(breakFlag==1)
                {
                    break;
                }
                if(breakFlag==2)
                {
                    for(;i1<path1.size()-1;)
                    {
                        total+=array_2d[i1+1][j1];
                        i1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                    }
                    break;
                }
            }
            //std::cout<<"此次子轨迹得分 "<<subtrajectotal<<endl;
            if(subtrajectotal>SimilarScore)
            {
                SimilarScore=subtrajectotal;
                start1=start_1;
                start2=start_2;
                end2=end_2;
                //std::cout<<"第一个子轨迹从 "<<start1+1<<" 到 "<<end1+1<<"第二个子轨迹从 "<<start2+1<<" 到 "<<end2+1<<" 代表性相似子轨迹得分是 "<<SimilarScore<<endl;
            }
            //std:cout<<"subtrajectotal "<<subtrajectotal<<endl;
        }
    }

    if(NoSimilar==true)
        return 0;
    else
        return -SimilarScore;
}


//query-data改进
double minSubTrajectory2query(const path& path1, const path& path2)
{
    std::vector<std::vector<int>> array_2d(path1.size(),std::vector<int>(path2.size()));
    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            if(edr1(path1[i],path2[j],0.0005)==1)
            {
                array_2d[i][j]=1;
            }
            else 
            {
                array_2d[i][j]=-1;
            }
        }
    }

    //输出创建的数组
    //for(int i=0;i<2;i++)
    //{
    //    for(int j=0;j<path2.size();j++)
    //    {
    //        std::cout<<array_2d[i][j];
    //    }
    //    std::cout<<endl;
    //}
    
    bool breakFlag=false;

    int SimilarScore=-10000;
    int start1=-1;
    int start2=-1;
    int end1=-1;
    int end2=-1;

    bool NoSimilar=true;

    for(int i=0;i<path1.size();i++)
    {
        for(int j=0;j<path2.size();j++)
        {
            NoSimilar=false;
            //std::cout<<"i "<<"j "<<i+1<<" "<<j+1<<endl;
            int total=2*array_2d[i][j];
            int subtrajectotal=total;
            int start_1=i;
            int start_2=j;
            int end_1=i;
            int end_2=j;

            int breakFlag=0;
            for(int i1=i;i1<path1.size()-1;)
            {
                int j1=j;
                //for(int j1=j;j1<path2.size()-1;)
                for(;j1<path2.size()-1;)
                {
                    //判断朝哪个方向
                    //std::cout<<"执行到这里0"<<endl;
                    if(array_2d[i1][j1]==1&&array_2d[i1+1][j1+1]==1)
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"right and down"<<endl;
                    }
                    else if(array_2d[i1][j1]==1&&array_2d[i1][j1+1]>=array_2d[i1+1][j1])
                    {
                        total-=1;
                        j1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"right"<<endl;
                    }
                    else if(array_2d[i1][j1]==1&&array_2d[i1+1][j1]>array_2d[i1][j1+1])
                    {
                        total-=1;
                        i1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"down"<<endl;
                    }
                    else if(array_2d[i1][j1]==-1&&array_2d[i1][j1+1]==1)
                    {
                        total+=3*array_2d[i1][j1+1];
                        j1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"right"<<endl;
                    }
                    else if(array_2d[i][j1]==-1&&array_2d[i1+1][j1]==1)
                    {
                        total+=3*array_2d[i1+1][j1];
                        i1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"down"<<endl;
                    }
                    else if(array_2d[i1][j1]!=array_2d[i1+1][j1+1])
                    {
                        total+=2*array_2d[i1+1][j1+1];
                        i1+=1;
                        j1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"right and down"<<endl;
                    }
                    else
                    {
                        total+=array_2d[i1][j1+1];
                        j1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                        //std::cout<<"right"<<endl;
                    }

                    //std::cout<<"total "<<total<<endl;
                    //std::cout<<"i1 "<<"j1 "<<i1+1<<" "<<j1+1<<endl;
                    if(i1+1>=path1.size())
                    {
                        breakFlag=1;
                        break;
                    }
                    if(j1+1>=path2.size())
                    {
                        breakFlag=2;
                        break;
                    }
                }
                if(j1+1>=path2.size())
                {
                    breakFlag=2;
                    break;
                }
                if(breakFlag==1)
                {
                    break;
                }
                if(breakFlag==2)
                {
                    for(;i1<path1.size()-1;)
                    {
                        total+=array_2d[i1+1][j1];
                        i1+=1;
                        subtrajectotal=total;
                        end_2=j1;
                    }
                    break;
                }
            }
            //std::cout<<"此次子轨迹得分 "<<subtrajectotal<<endl;
            if(subtrajectotal>SimilarScore)
            {
                SimilarScore=subtrajectotal;
                start1=start_1;
                start2=start_2;
                end2=end_2;
                //std::cout<<"第一个子轨迹从 "<<start1+1<<" 到 "<<end1+1<<"第二个子轨迹从 "<<start2+1<<" 到 "<<end2+1<<" 代表性相似子轨迹得分是 "<<SimilarScore<<endl;
            }
            //std:cout<<"subtrajectotal "<<subtrajectotal<<endl;
        }
    }

    if(NoSimilar==true)
        return 0;
    else
        return -SimilarScore;
}

double minSubTrajectory(const string & algorithm, const path& path1, const path& path2) {
    double minDistance = 1000000;
    //std::cout<<path1.size()<<endl;
    //std::cout<<path2.size()<<endl;

    for(int ii=0; ii<path1.size()-1;++ii)
    {
        for(int jj=ii+1;jj<path1.size();++jj)
        {
            //std::cout<<"执行一条path1的子轨迹"<<endl;
            path subpath1(path1.begin()+ii,path1.begin()+jj+1);
            
            for (int i = 0; i < path2.size() - 1; ++i) 
            {
                for (int j = i + 1; j < path2.size(); ++j) 
                {   
                    //std::cout<<"执行wedDistance"<<endl;
                    minDistance = min(wedDistance(subpath1, path2, i, j), minDistance);
                }
            }

        }
    }
    return minDistance;
}

double efficientAlgorithmWED_data_to_data(path path1, path path2) {
    double final_r=100000;
    for(int path1_i=0;path1_i<path1.size();path1_i++)
    {
        for(int path1_j=path1_i+1;path1_j<path1.size()+1;path1_j++)
        {
            auto path1_start=path1.begin()+path1_i;
            auto path1_end=path1.begin()+path1_j;

            int lensubPath1=path1_end-path1_start;
            int lenPath2=path2.size();
            auto allCost=new double[path2.size()];
            auto allCostTmp=new double[path2.size()];
            auto starts=new int[path2.size()];
            auto startsTmp=new int[path2.size()];
            double empty[lensubPath1];
            for(int i=0;i<lensubPath1;++i)
            {
                if(i==0)
                {
                    empty[i]=distance(nullPoint,*(path1_start+i));
                }
                else 
                {
                    empty[i]=empty[i-1]+distance(nullPoint,*(path1_start+i));
                }
            }
            for(int i=0;i<lensubPath1;++i)
            {
                for(int j=0;j<lenPath2;++j)
                {
                    if(i==0)
                    {
                        starts[j]=j;
                        allCost[j]=distance(*(path1_start),path2[j]);
                    }
                    else 
                    {
                        if(j==0)
                        {
                            allCost[0]=min(allCostTmp[0]+distance(nullPoint,*(path1_start+i)),empty[i-1]+distance(*(path1_start+i),path2[0]));
                            starts[0]=0;
                            continue;
                        }
                        double c1=allCostTmp[j]+distance(nullPoint,*(path1_start+i));
                        double c2=allCostTmp[j-1]+distance(*(path1_start+i),path2[j]);
                        double c3=allCost[j-1]+distance(nullPoint,path2[j-1])-distance(*(path1_start+i),path2[j-1])+(c2-allCostTmp[j-1]);
                        if(c1>=c2)
                        {
                            allCost[j]=c2;
                            starts[j]=startsTmp[j-1];
                        }
                        else 
                        {
                            allCost[j]=c1;
                            starts[j]=startsTmp[j];
                        }
                        if(c3<allCost[j])
                        {
                            starts[j]=starts[j-1];
                            allCost[j]=c3;
                        }
                    }
                }
                swap(allCost,allCostTmp);
                swap(starts, startsTmp);
            }
            int start,end=0;
            double res=allCostTmp[0];
            for(int j=0;j<lenPath2;++j)
            {
                if(res>allCostTmp[j])
                {
                    res=min(res,allCostTmp[j]);
                    end=j;
                }
            }
            start=startsTmp[end];
            delete[] allCostTmp;
            delete[] allCost;
            delete[] starts;
            delete[] startsTmp;
            if(final_r>res)
            {
                final_r=res;
            }
        }
    }
    return final_r;
}


subResult efficientAlgorithmWED(path path1, path path2) {
    int lenPath1 = path1.size();
    int lenPath2 = path2.size();
    auto allCost = new double[path2.size()];
    auto allCostTmp = new double[path2.size()];
    auto starts = new int[path2.size()];
    auto startsTmp = new int[path2.size()];
    double empty[lenPath1];
    for (int i = 0; i < lenPath1; ++i) {
        if (i == 0) {
            empty[i] = distance(nullPoint, path1[i]);
        } else {
            empty[i] = empty[i - 1] + distance(nullPoint, path1[i]);
        }
    }
    for (int i = 0; i < lenPath1; ++i) {
        //        cout << "from:  ";
        for (int j = 0; j < lenPath2; ++j) {
            if (i == 0) {
                starts[j] = j;
                allCost[j] = distance(path1[0], path2[j]);
            } else{
                //                int x = 0;
                if (j == 0) {
                    allCost[0] = min(allCostTmp[0] + distance(nullPoint, path1[i]), empty[i-1] + distance(path1[i], path2[0]));
                    starts[0] = 0;
                    continue;
                }
                double c1 = allCostTmp[j] + distance(nullPoint, path1[i]);
                double c2 = allCostTmp[j - 1] + distance(path1[i], path2[j]);
                double c3 = allCost[j - 1] + distance(nullPoint, path2[j - 1]) - distance(path1[i], path2[j-1])
                        + (c2 - allCostTmp[j-1]);
                if (c1 >= c2) {
                    allCost[j] = c2;
                    starts[j] = startsTmp[j - 1];
                    //                    x = 0;
                } else {
                    allCost[j] = c1;
                    starts[j] = startsTmp[j];
                    //                    x = 1;
                }
                if (c3 < allCost[j]) {
                    starts[j] = starts[j-1];
                    allCost[j] = c3;
                    //                    x = -1;
                }
                //                cout << x << " ";

                /*
                 * 8945128451248657
                 * 9481659841637894616841131216574121965023657864513246364124863125628965130
                 */

            }
        }
        //        cout << endl;
        //        for (int j = 0; j < path2.size(); ++j) {
        //            cout << *(allCost + j) << " ";
        //        }
        //        cout << endl;
        swap(allCost,allCostTmp);
        swap(starts, startsTmp);
    }
    int start, end = 0;
    double res = allCostTmp[0];
    for (int j = 0; j < lenPath2; ++j) {
        if (res > allCostTmp[j]) {
            res = min(res, allCostTmp[j]);
            end = j;
        }
    }
    start = startsTmp[end];
    delete[] allCostTmp;
    delete[] allCost;
    delete[] starts;
    delete[] startsTmp;
    subResult r;
    r.first.first = start;
    r.first.second = end;
    r.second = res;
    return r;
}


double efficientAlgorithmDTW_data_to_data(path path1, path path2) {
    double final_r=100000;
    for(int path1_i=0;path1_i<path1.size();path1_i++)
    {
        for(int path1_j=path1_i+1;path1_j<path1.size()+1;path1_j+1)
        {
            auto path1_start=path1.begin()+path1_i;
            auto path1_end=path1.begin()+path1_j;

            int lensubPath1=path1_end-path1_start;
            int lenPath2=path2.size();
            auto allCost=new double[path2.size()];
            auto allCostTmp=new double[path2.size()];
            auto starts=new int[path2.size()];
            auto startsTmp=new int[path2.size()];
            double empty[lensubPath1];
            for(int i=0;i<lensubPath1;++i)
            {
                if(i==0)
                {
                    empty[i]=pointDistance(path2[0],*(path1_start+i));
                }
                else 
                {
                    empty[i]=empty[i-1]+pointDistance(path2[0],*(path1_start+i));
                }
            }
            for(int i=0;i<lensubPath1;++i)
            {
                for(int j=0;j<lenPath2;++j)
                {
                    if(i==0)
                    {
                        starts[j]=j;
                        allCost[j]=pointDistance(*(path1_start),path2[j]);
                    }
                    else 
                    {
                        if(j==0)
                        {
                            starts[0]=0;
                            allCost[0]=min(allCostTmp[0]+pointDistance(path2[0],*(path1_start+i)),empty[i-1]+pointDistance(*(path1_start+i),path2[0]));
                            continue;
                        }
                        if(allCostTmp[j]>=allCostTmp[j-1])
                        {
                            starts[j]=startsTmp[j-1];
                        }
                        else 
                        {
                            starts[j]=startsTmp[j];
                        }
                        if(allCostTmp[j]>allCost[j-1]&&allCostTmp[j-1]>allCost[j-1])
                        {
                            starts[j]=starts[j-1];
                        }
                        allCost[j]=min(allCostTmp[j],min(allCostTmp[j-1],allCost[j-1]))+pointDistance(*(path1_start+i),path2[j]);
                    }
                }
                swap(allCost,allCostTmp);
                swap(starts,startsTmp);
            }
            int start;
            int end=0;
            double res=allCostTmp[0];
            for(int j=0;j<lenPath2;++j)
            {
                if(res>allCostTmp[j])
                {
                    res=min(res,allCostTmp[j]);
                    end=j;
                }
            }
            start=startsTmp[end];
            delete[] allCostTmp;
            delete[] allCost;
            delete[] starts;
            delete[] startsTmp;
            if(final_r>res)
            {
                final_r=res;
            }
        }
    }
    return final_r;
}


subResult efficientAlgorithmDTW(path path1, path path2) {
    int lenPath1 = path1.size();
    int lenPath2 = path2.size();
    auto allCost = new double[path2.size()];
    auto allCostTmp = new double[path2.size()];
    auto starts = new int[path2.size()];
    auto startsTmp = new int[path2.size()];
    double empty[lenPath1];
    for (int i = 0; i < lenPath1; ++i) {
        if (i == 0) {
            empty[i] = pointDistance(path2[0], path1[i]);
        } else {
            empty[i] = empty[i - 1] + pointDistance(path2[0], path1[i]);
        }
    }
    for (int i = 0; i < lenPath1; ++i) {
        for (int j = 0; j < lenPath2; ++j) {
            if (i == 0) {
                starts[j] = j;
                allCost[j] = pointDistance(path1[0], path2[j]);
            } else {
                if (j == 0) {
                    starts[0] = 0;
                    allCost[0] = min(allCostTmp[0] + pointDistance(path2[0], path1[i]), empty[i-1] + pointDistance(path1[i], path2[0]));
                    continue;
                }
                if (allCostTmp[j] >= allCostTmp[j - 1]) {
                    starts[j] = startsTmp[j - 1];
                } else {
                    starts[j] = startsTmp[j];
                }
                if (allCostTmp[j] > allCost[j-1] && allCostTmp[j - 1] > allCost[j-1]) {
                    starts[j] = starts[j - 1];
                }
                allCost[j] = min(allCostTmp[j], min(allCostTmp[j - 1], allCost[j-1])) + pointDistance(path1[i], path2[j]);
            }
        }
        swap(allCost,allCostTmp);
        swap(starts, startsTmp);
    }
    int start;
    int end = 0;
    double res = allCostTmp[0];
    for (int j = 0; j < lenPath2; ++j) {
        if (res > allCostTmp[j]) {
            res = min(res, allCostTmp[j]);
            end = j;
        }
    }
    start = startsTmp[end];
    delete[] allCostTmp;
    delete[] allCost;
    delete[] starts;
    delete[] startsTmp;
    subResult r;
    r.first.first = start;
    r.first.second = end;
    r.second = res;
    return r;
}


double efficientAlgorithmFC_data_to_data(path path1, path path2)
{
    double final_r=100000;
    for(int path1_i=0;path1_i<path1.size();path1_i++)
    {
        for(int path1_j=path1_i+1;path1_j<path1.size()+1;path1_j++)
        {
            auto path1_start=path1.begin()+path1_i;
            auto path1_end=path1.begin()+path1_j;

            int lensubPath1=path1_end-path1_start;
            int lenPath2=path2.size();
            auto allCost=new double[path2.size()];
            auto allCostTmp=new double[path2.size()];
            auto starts=new int[path2.size()];
            auto startsTmp=new int[path2.size()];
            double empty[lensubPath1];
            for(int i=0;i<lensubPath1;++i)
            {
                if(i==0)
                {
                    empty[i]=pointDistance(path2[0],*(path1_start+i));
                }
                else 
                {
                    empty[i]=max(empty[i-1],pointDistance(path2[0],*(path1_start+i)));
                }
            }
            for(int i=0;i<lensubPath1;++i)
            {
                for(int j=0;j<lenPath2;++j)
                {
                    if(i==0)
                    {
                        starts[j]=j;
                        allCost[j]=pointDistance(*(path1_start),path2[j]);
                    }
                    else 
                    {
                        if(j==0)
                        {
                            starts[0]=0;
                            allCost[0]=min(max(allCostTmp[0],pointDistance(path2[0],*(path1_start+i))),max(empty[i-1],pointDistance(*(path1_start+i),path2[0])));
                            continue;
                        }
                        if(allCostTmp[j]>=allCostTmp[j-1])
                        {
                            starts[j]=startsTmp[j];
                        }
                        else
                        {
                            starts[j]=startsTmp[j];
                        }
                        if(allCostTmp[j]>allCost[j-1]&&allCostTmp[j-1]>allCost[j-1])
                        {
                            starts[j]=starts[j-1];
                        }
                        allCost[j]=max(min(allCostTmp[j],min(allCostTmp[j-1],allCost[j-1])),pointDistance(path1[i],path2[j]));
                    }
                }
                swap(allCost,allCostTmp);
                swap(starts,startsTmp);
            }
            int start;
            int end=0;
            double res=allCostTmp[0];
            for(int j=0;j<lenPath2;++j)
            {
                if(res>allCostTmp[j])
                {
                    res=min(res,allCostTmp[j]);
                    end=j;
                }
            }
            start=startsTmp[end];
            delete[] allCostTmp;
            delete[] allCost;
            delete[] starts;
            delete[] startsTmp;
            if(final_r>res)
            {
                final_r=res;
            }
        }
    }
    return final_r;
}




subResult efficientAlgorithmFC(path path1, path path2) {
    int lenPath1 = path1.size();
    int lenPath2 = path2.size();
    auto allCost = new double[path2.size()];
    auto allCostTmp = new double[path2.size()];
    auto starts = new int[path2.size()];
    auto startsTmp = new int[path2.size()];
    double empty[lenPath1];
    for (int i = 0; i < lenPath1; ++i) {
        if (i == 0) {
            empty[i] = pointDistance(path2[0], path1[i]);
        } else {
            empty[i] = max(empty[i - 1], pointDistance(path2[0], path1[i]));
        }
    }
    for (int i = 0; i < lenPath1; ++i) {
        for (int j = 0; j < lenPath2; ++j) {
            if (i == 0) {
                starts[j] = j;
                allCost[j] = pointDistance(path1[0], path2[j]);
            } else {
                if (j == 0) {
                    starts[0] = 0;
                    allCost[0] = min(max(allCostTmp[0], pointDistance(path2[0], path1[i])), max(empty[i-1],pointDistance(path1[i], path2[0])));
                    continue;
                }
                if (allCostTmp[j] >= allCostTmp[j - 1]) {
                    starts[j] = startsTmp[j - 1];
                } else {
                    starts[j] = startsTmp[j];
                }
                if (allCostTmp[j] > allCost[j-1] && allCostTmp[j - 1] > allCost[j-1]) {
                    starts[j] = starts[j - 1];
                }
                allCost[j] = max(min(allCostTmp[j], min(allCostTmp[j - 1], allCost[j-1])), pointDistance(path1[i], path2[j]));
            }
        }
        swap(allCost,allCostTmp);
        swap(starts, startsTmp);
    }
    int start;
    int end = 0;
    double res = allCostTmp[0];
    for (int j = 0; j < lenPath2; ++j) {
        if (res > allCostTmp[j]) {
            res = min(res, allCostTmp[j]);
            end = j;
        }
    }
    start = startsTmp[end];
    delete[] allCostTmp;
    delete[] allCost;
    delete[] starts;
    delete[] startsTmp;
    subResult r;
    r.first.first = start;
    r.first.second = end;
    r.second = res;
    return r;
}


double efficientAlgorithm_data_to_data(const string& algorithm, const path& path1, const path& path2) {
    if (matricsType == "dtw") {
        return efficientAlgorithmDTW_data_to_data(path1, path2);
    } 
    else if(matricsType == "FC") {
        return efficientAlgorithmFC_data_to_data(path1, path2);
    } else {
    return efficientAlgorithmWED_data_to_data(path1, path2);
    }
}


subResult efficientAlgorithm(const string& algorithm, const path& path1, const path& path2) {
    if (matricsType == "dtw") {
        return efficientAlgorithmDTW(path1, path2);
    } else if(matricsType == "FC") {
        return efficientAlgorithmFC(path1, path2);
    } else {
        return efficientAlgorithmWED(path1, path2);
    }
}
//}
