//
// Created by 86173 on 2022/1/27.
//

#include "distance.h"
#include "exactS.h"

//自己添加的内容

//namespace gnns 
//{
string matricsType = "dtw"; // erp, edr, FC, dtw

int minThreshold = 3;


string pruningType = "pointprune";
string gatherType = "gridbase";
string dataType = "xian";


int maxLen = 100;
int minLen = 80;


double AR = 0;
double MR = 0;
double RR = 0;




double pointDistance(rtree_point p1, rtree_point p2) {
    return sqrt((p1.get<0>() - p2.get<0>()) * (p1.get<0>() - p2.get<0>()) + (p1.get<1>() - p2.get<1>()) * (p1.get<1>() - p2.get<1>()));
//    return p1.first == p2.first ? 0:1;
}


double levenshteinDistance(rtree_point p1, rtree_point p2) {
    if (p1.get<0>() == p2.get<0>() && p1.get<1>() == p2.get<1>()) return 0;
    return 1;
}

double lcss(rtree_point p1, rtree_point p2, double e) {
    if (abs(p1.get<0>() - p2.get<0>()) < e && abs(p1.get<1>() - p2.get<1>()) < e) return 0;
    return 1;
}

double erp(rtree_point p1, rtree_point p2) {
    if (p1.get<0>() == nullPoint.get<0>()) {
        p1 = centerPoint;
    }
    if (p2.get<0>() == nullPoint.get<0>()) {
        p2 = centerPoint;
    }
    double distance = sqrt((p1.get<0>() - p2.get<0>()) * (p1.get<0>() - p2.get<0>()) + (p1.get<1>() - p2.get<1>()) * (p1.get<1>() - p2.get<1>()));
    return distance;
}


double edr(rtree_point p1, rtree_point p2, double e) {
    double distance = sqrt((p1.get<0>() - p2.get<0>()) * (p1.get<0>() - p2.get<0>()) + (p1.get<1>() - p2.get<1>()) * (p1.get<1>() - p2.get<1>()));
    if (e < distance) return 1;
    return 0;
}

double distance(rtree_point p1, rtree_point p2) {
    if (matricsType == "levenshteinDistance") return levenshteinDistance(p1, p2);
    if (matricsType == "edr") return edr(p1, p2, 0.0005);
    if (matricsType == "erp") return erp(p1, p2);
    if (matricsType == "lcss") return lcss(p1, p2, 0.0005);
    return 0;
}



double wedDistance(path path1, path path2, int start, int end) {
    auto distanceMatrix = new double[end - start + 1];
    auto distanceMatrixTmp = new double[end - start + 1];
    distanceMatrixTmp[0] = 0;
    //std::cout<<"执行第一个for循环"<<endl;
    for (int j = start + 1; j < end + 1; ++j) {
        distanceMatrixTmp[j - start] = distanceMatrixTmp[j - start - 1] + distance(path2[j - 1], nullPoint);
    }
    //std::cout<<"执行两层for循环"<<endl;
    for (int i = 1; i < path1.size() + 1; ++i) {
        for (int j = start; j < end + 1; ++j) {
            if (j == start) {
                distanceMatrix[0] = distanceMatrixTmp[0] + distance(path1[i-1], nullPoint);
                continue;
            }
            distanceMatrix[j - start] = min(distanceMatrixTmp[j - start - 1] + distance(path1[i-1], path2[j-1]),
                                               min(distanceMatrix[j - start - 1] + distance(nullPoint, path2[j-1]),
                                                   distanceMatrixTmp[j - start] + distance(nullPoint, path1[i-1])));
        }
        swap(distanceMatrixTmp, distanceMatrix);
    }
    //std::cout<<"两层for循环执行结束"<<endl;
    double result = distanceMatrixTmp[end - start];
    delete[] distanceMatrix;
    delete[] distanceMatrixTmp;
    return result;
}

double dtwDistance(path path1, path path2, int start, int end) {
    auto distanceMatrix = new double[end - start + 1];
    auto distanceMatrixTmp = new double[end - start + 1];
    distanceMatrixTmp[0] = 0;
    for (int j = start + 1; j < end + 1; ++j) {
        distanceMatrixTmp[j - start] = distanceMatrixTmp[j - start - 1] + pointDistance(path2[j - 1], path1[0]);
    }
    for (int i = 1; i < path1.size() + 1; ++i) {
        for (int j = start; j < end + 1; ++j) {
            if (j == start) {
                distanceMatrix[0] = distanceMatrixTmp[0] + pointDistance(path1[i-1], path2[start]);
                continue;
            }
            distanceMatrix[j - start] = min(distanceMatrixTmp[j - start - 1] + pointDistance(path1[i-1], path2[j-1]),
                                            min(distanceMatrix[j - start - 1] + pointDistance(path1[i-1], path2[j-1]),
                                                distanceMatrixTmp[j - start] + pointDistance(path1[i-1], path2[j-1])));
        }
        swap(distanceMatrixTmp, distanceMatrix);
    }
    double result = distanceMatrixTmp[end - start];
    delete[] distanceMatrix;
    delete[] distanceMatrixTmp;
    return result;
}
//}

