#ifndef TRAJCSIMILAR_DISTANCE_H
#define TRAJCSIMILAR_DISTANCE_H
#include "iostream"
#include <vector>
#include <cmath>
#include <numeric>
#include "constant.h"
//自己添加的内容
//namespace gnns 
//{

//extern int maxLen;
//extern int minLen;
double distance(rtree_point p1, rtree_point p2);

double pointDistance(rtree_point p1, rtree_point p2);
double levenshteinDistance(rtree_point p1, rtree_point p2);
double edr(rtree_point p1, rtree_point p2, double e=0);
double wedDistance(path path1, path path2, int start, int end);
double dtwDistance(path path1, path path2, int start, int end);

//}
#endif //TRAJCSIMILAR_DISTANCE_H
