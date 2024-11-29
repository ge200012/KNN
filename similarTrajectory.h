//
// Created by 86173 on 2022/1/27.
//

#ifndef TRAJCSIMILAR_SIMILARTRAJECTORY_H
#define TRAJCSIMILAR_SIMILARTRAJECTORY_H

#include "distance.h"

//自己添加的内容
//namespace gnns
//{

using namespace std;
extern string matricsType;
extern string pruningType;
extern string gatherType;
extern string dataType;

//extern int maxLen;
//extern int minLen;
double minSubTrajectory(const string & algorithm, const path& path1, const path& path2);
//自己添加的内容
int edr1(rtree_point p1, rtree_point p2, double e);
double minSubTrajectory1(const string & algorithm, const path& path1, const path& path2);
Result minSubTrajectory2(const string & algorithm, const path& path1, const path& path2);
double minSubTrajectory1query(const path& path1, const path& path2);
double minSubTrajectory2query(const path& path1, const path& path2);
double efficientAlgorithm_data_to_data(const string& algorithm, const path& path1, const path& path2);

subResult efficientAlgorithm(const string& algorithm, const path& path1, const path& path2);

//}
#endif //TRAJCSIMILAR_SIMILARTRAJECTORY_H
