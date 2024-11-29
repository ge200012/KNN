//
// Created by 86173 on 2022/2/3.
//

#ifndef TRAJCSIMILAR_CONSTANT_H
#define TRAJCSIMILAR_CONSTANT_H
#include "distance.h"
#include "map"

//自己添加的内容
#include "utility"
#include "vector"
#include "string"

#include<boost/geometry.hpp>
#include<boost/geometry/index/rtree.hpp>
#include<boost/geometry/geometries/linestring.hpp>
#include<boost/geometry/geometries/box.hpp>
namespace bg=boost::geometry;
namespace bgi=boost::geometry::index;





//namespace gnns 
//{
using namespace std;


// 定义轨迹点和轨迹类型
typedef bg::model::point<double, 2, bg::cs::cartesian> rtree_point;
typedef bg::model::linestring<rtree_point> path;
typedef bg::model::box<rtree_point> box_t;

// 创建R-tree的节点类型
typedef std::pair<box_t, unsigned> value;
typedef bgi::rtree< value, bgi::quadratic<16> > rtree_t;


//typedef pair<double,double> point;
//typedef vector<point> path;


typedef pair<pair<int, int>, double> subResult;
typedef pair<double,pair< pair<int,int>, pair<int,int> > > Result;
//自己添加的内容
typedef vector<vector<pair<int, double>>> kNN;



static map<string , pair<vector<double>, vector<double>>> range{
    {"xian", {{34.20, 34.29}, {108.91, 109.00}}},
    {"chengdu", {{30.65, 30.73}, {104.04, 104.13}}},
    {"porto", {{41.11, 41.19}, {-8.67, -8.57}}}
};
static bool generateResult = false;

static double MaxSimilar = 100000000000;

static string filepath = "/home/jiabaojin/project/trajectorySimilarity/data/result/";
static const char *datasource = "/home/jiabaojin/project/trajectorySimilarity/data/data/%s/trajectory_data";

//自己添加的内容
static rtree_point nullPoint = rtree_point(-1000.0f, -1000.0f);
static rtree_point centerPoint = rtree_point(34.44f, 109.95f); //xian
//static rtree_point centerPoint = rtree_point(30.69f, 104.085f); //chengdu
//static rtree_point centerPoint = rtree_point(41.15f, -8.62f); //porto

//static pair<double, double> nullPoint = pair<double, double>{-1000,-1000};
//static pair<double, double> centerPoint = pair<double, double>{34.44, 109.95}; // 西安
//}
#endif //TRAJCSIMILAR_CONSTANT_H
