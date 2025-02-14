#include <algorithm>
#include "pss.h"

//自己添加的内容

//namespace gnns
//{
int pssstart = -1;
int pssend = -1;

subResult pssWED(path p1, path p2) {

    int start = 0;
    int end = p2.size();

    auto distanceMatrix = new double[end - start + 1];
    auto distanceMatrixTmp = new double[end - start + 1];
    distanceMatrixTmp[0] = 0;
    for (int j = start + 1; j < end + 1; ++j) {
        distanceMatrixTmp[j - start] = distanceMatrixTmp[j - start - 1] + distance(p2[j - 1], nullPoint);
    }
    for (int i = 1; i < p1.size() + 1; ++i) {
        for (int j = start; j < end + 1; ++j) {
            if (j == start) {
                distanceMatrix[0] = distanceMatrixTmp[0] + distance(p1[i-1], nullPoint);
                continue;
            }
            distanceMatrix[j - start] = min(distanceMatrixTmp[j - start - 1] + distance(p1[i-1], p2[j-1]),
                                            min(distanceMatrix[j - start - 1] + distance(nullPoint, p2[j-1]),
                                                distanceMatrixTmp[j - start] + distance(nullPoint, p1[i-1])));
        }
        swap(distanceMatrixTmp, distanceMatrix);
    }
    delete[] distanceMatrix;
    distanceMatrix = new double[p1.size() + 1];
    reverse(p1.begin(), p1.end());
    reverse(p2.begin(), p2.end());

    double bestScore = 100000000;
    start = 0;

    auto distanceMatrixReverse = new double[p1.size() + 1];
    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrix[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrix[i] = distanceMatrix[i-1] + distance(p1[i-1], nullPoint);
            }
        }
        distanceMatrixReverse[0] = distanceMatrix[0] + distance(p2[j - 1], nullPoint);
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(distanceMatrix[i-1] + distance(p1[i-1], p2[j-1]),
                                           min(distanceMatrix[i] + distance(nullPoint, p2[j-1]),
                                               distanceMatrixReverse[i-1] + distance(nullPoint, p1[i-1])));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        double theta2 = distanceMatrixTmp[p2.size() - j + 1];
        if (min(theta1,theta2) < bestScore) {
            bestScore = min(theta1,theta2);
            if (theta1 < theta2) {
                pssstart = start;
                pssend = j - 1;
            } else {
                pssstart = j - 1;
                pssend = p2.size() - 1;
            }
            start = j;
        }
        swap(distanceMatrix, distanceMatrixReverse);
    }
    delete[] distanceMatrix;
    delete[] distanceMatrixReverse;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}


subResult pssDTW(path p1, path p2) {

    int start = 0;
    int end = p2.size();
//    reverse(p1.begin(), p1.end());
//    reverse(p2.begin(), p2.end());
    auto distanceMatrix = new double[end - start + 1];
    auto distanceMatrixTmp = new double[end - start + 1];
    distanceMatrixTmp[0] = 0;
    for (int j = start + 1; j < end + 1; ++j) {
        distanceMatrixTmp[j - start] = distanceMatrixTmp[j - start - 1] + pointDistance(p2[j - 1], p1[0]);
    }
    for (int i = 1; i < p1.size() + 1; ++i) {
        for (int j = start; j < end + 1; ++j) {
            if (j == start) {
                distanceMatrix[0] = distanceMatrixTmp[0] + pointDistance(p1[i-1], p2[start]);
                continue;
            }
            distanceMatrix[j - start] = min(distanceMatrixTmp[j - start - 1] + pointDistance(p1[i-1], p2[j-1]),
                                            min(distanceMatrix[j - start - 1] + pointDistance(p1[i-1], p2[j-1]),
                                                distanceMatrixTmp[j - start] + pointDistance(p2[j-1], p1[i-1])));
        }
        swap(distanceMatrixTmp, distanceMatrix);
    }
    reverse(p1.begin(), p1.end());
    reverse(p2.begin(), p2.end());
    double bestScore = 100000000;
    start = 0;

    auto distanceMatrixReverse = new double[p1.size() + 1];
    auto distanceMatrixReverseTmp = new double[p1.size() + 1];

    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrixReverseTmp[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrixReverseTmp[i] = distanceMatrixReverseTmp[i-1] + pointDistance(p1[i-1], p2[start]);
            }
        }
        distanceMatrixReverse[0] = distanceMatrixReverseTmp[0] + pointDistance(p2[j - 1], p1[0]);
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(distanceMatrixReverseTmp[i-1] + pointDistance(p1[i-1], p2[j-1]),
                                           min(distanceMatrixReverseTmp[i] + pointDistance(p1[i-1], p2[j-1]),
                                               distanceMatrixReverse[i-1] + pointDistance(p2[j-1], p1[i-1])));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        double theta2 = distanceMatrixTmp[p2.size() - j + 1];
        if (min(theta1,theta2) < bestScore) {
            bestScore = min(theta1,theta2);
            if (theta1 < theta2) {
                pssstart = start;
                pssend = j - 1;
            } else {
                pssstart = j - 1;
                pssend = p2.size() - 1;
            }
            start = j;
        }
        swap(distanceMatrixReverseTmp, distanceMatrixReverse);

    }
    delete[] distanceMatrixReverse;
    delete[] distanceMatrixReverseTmp;
    delete[] distanceMatrixTmp;
    delete[] distanceMatrix;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}

subResult pssFC(path p1, path p2) {

    int start = 0;
    int end = p2.size();
    //    reverse(p1.begin(), p1.end());
    //    reverse(p2.begin(), p2.end());
    auto distanceMatrix = new double[end - start + 1];
    auto distanceMatrixTmp = new double[end - start + 1];
    distanceMatrixTmp[0] = 0;
    for (int j = start + 1; j < end + 1; ++j) {
        distanceMatrixTmp[j - start] = max(distanceMatrixTmp[j - start - 1], pointDistance(p2[j - 1], p1[0]));
    }
    for (int i = 1; i < p1.size() + 1; ++i) {
        for (int j = start; j < end + 1; ++j) {
            if (j == start) {
                distanceMatrix[0] = max(distanceMatrixTmp[0], pointDistance(p1[i-1], p2[start]));
                continue;
            }
            distanceMatrix[j - start] = min(max(distanceMatrixTmp[j - start - 1], pointDistance(p1[i-1], p2[j-1])),
                                            min(max(distanceMatrix[j - start - 1], pointDistance(p1[i-1], p2[j-1])),
                                                max(distanceMatrixTmp[j - start], pointDistance(p2[j-1], p1[i-1]))));
        }
        swap(distanceMatrixTmp, distanceMatrix);
    }
    reverse(p1.begin(), p1.end());
    reverse(p2.begin(), p2.end());
    double bestScore = 100000000;
    start = 0;

    auto distanceMatrixReverse = new double[p1.size() + 1];
    auto distanceMatrixReverseTmp = new double[p1.size() + 1];

    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrixReverseTmp[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrixReverseTmp[i] = max(distanceMatrixReverseTmp[i-1], pointDistance(p1[i-1], p2[start]));
            }
        }
        distanceMatrixReverse[0] = max(distanceMatrixReverseTmp[0], pointDistance(p2[j - 1], p1[0]));
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(max(distanceMatrixReverseTmp[i-1], pointDistance(p1[i-1], p2[j-1])),
                                           min(max(distanceMatrixReverseTmp[i], pointDistance(p1[i-1], p2[j-1])),
                                               max(distanceMatrixReverse[i-1], pointDistance(p2[j-1], p1[i-1]))));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        double theta2 = distanceMatrixTmp[p2.size() - j + 1];
        if (min(theta1,theta2) < bestScore) {
            bestScore = min(theta1,theta2);
            if (theta1 < theta2) {
                pssstart = start;
                pssend = j - 1;
            } else {
                pssstart = j - 1;
                pssend = p2.size() - 1;
            }
            start = j;
        }
        swap(distanceMatrixReverseTmp, distanceMatrixReverse);

    }
    delete[] distanceMatrixReverse;
    delete[] distanceMatrixReverseTmp;
    delete[] distanceMatrixTmp;
    delete[] distanceMatrix;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}


subResult pss(const path& path1, const path& path2) {
    if (matricsType == "dtw") {
        return pssDTW(path1, path2);
    } else if (matricsType == "FC") {
        return pssFC(path1, path2);
    } else {
        return pssWED(path1, path2);
    }
}


int d = 0;
int max_d = 10;

subResult posWED(path p1, path p2) {
    int start = 0;
    double bestScore = 100000000;
    auto distanceMatrix = new double[p1.size() + 1];
    auto distanceMatrixReverse = new double[p1.size() + 1];
    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrix[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrix[i] = distanceMatrix[i-1] + distance(p1[i-1], nullPoint);
            }
        }
        distanceMatrixReverse[0] = distanceMatrix[0] + distance(p2[j - 1], nullPoint);
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(distanceMatrix[i-1] + distance(p1[i-1], p2[j-1]),
                                           min(distanceMatrix[i] + distance(nullPoint, p2[j-1]),
                                               distanceMatrixReverse[i-1] + distance(nullPoint, p1[i-1])));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        if (theta1 < bestScore) {
            bestScore = theta1;
            pssstart = start;
            pssend = j - 1;
            d = 0;
        } else {
            d ++;
            if (d > max_d && pssstart == start){
                start = pssend + 1;
                j = start;
                d = 0;
            }
        }
        swap(distanceMatrix, distanceMatrixReverse);
    }
    delete[] distanceMatrix;
    delete[] distanceMatrixReverse;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}


subResult posDTW(path p1, path p2) {
    int start = 0;
    double bestScore = 100000000;
    start = 0;

    auto distanceMatrixReverse = new double[p1.size() + 1];
    auto distanceMatrixReverseTmp = new double[p1.size() + 1];
    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrixReverseTmp[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrixReverseTmp[i] = distanceMatrixReverseTmp[i-1] + pointDistance(p1[i-1], p2[start]);
            }
        }
        distanceMatrixReverse[0] = distanceMatrixReverseTmp[0] + pointDistance(p2[j - 1], p1[0]);
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(distanceMatrixReverseTmp[i-1] + pointDistance(p1[i-1], p2[j-1]),
                                           min(distanceMatrixReverseTmp[i] + pointDistance(p1[i-1], p2[j-1]),
                                               distanceMatrixReverse[i-1] + pointDistance(p2[j-1], p1[i-1])));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        if (theta1 < bestScore) {

            bestScore = theta1;
            pssstart = start;
            pssend = j - 1;
            start = j;
        } else {
            d ++;
            if (d > max_d && pssstart == start){
                start = pssend + 1;
                j = start;
                d = 0;
            }
        }
        swap(distanceMatrixReverseTmp, distanceMatrixReverse);
    }
    delete[] distanceMatrixReverseTmp;
    delete[] distanceMatrixReverse;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}


subResult posFC(path p1, path p2) {
    int start = 0;
    double bestScore = 100000000;
    start = 0;

    auto distanceMatrixReverse = new double[p1.size() + 1];
    auto distanceMatrixReverseTmp = new double[p1.size() + 1];
    for (int j = 1; j < p2.size() + 1; ++j) {
        if (j == start + 1) {
            distanceMatrixReverseTmp[0] = 0;
            for (int i = 1; i < p1.size() + 1; ++i) {
                distanceMatrixReverseTmp[i] = max(distanceMatrixReverseTmp[i-1], pointDistance(p1[i-1], p2[start]));
            }
        }
        distanceMatrixReverse[0] = distanceMatrixReverseTmp[0] + pointDistance(p2[j - 1], p1[0]);
        for (int i = 1; i < p1.size() + 1; ++i) {
            distanceMatrixReverse[i] = min(max(distanceMatrixReverseTmp[i-1], pointDistance(p1[i-1], p2[j-1])),
                                           min(max(distanceMatrixReverseTmp[i], pointDistance(p1[i-1], p2[j-1])),
                                               max(distanceMatrixReverse[i-1],pointDistance(p2[j-1], p1[i-1]))));
        }
        double theta1 = distanceMatrixReverse[p1.size()];
        if (theta1 < bestScore) {

            bestScore = theta1;
            pssstart = start;
            pssend = j - 1;
            start = j;
        } else {
            d ++;
            if (d > max_d && pssstart == start){
                start = pssend + 1;
                j = start;
                d = 0;
            }
        }
        swap(distanceMatrixReverseTmp, distanceMatrixReverse);
    }
    delete[] distanceMatrixReverseTmp;
    delete[] distanceMatrixReverse;
    subResult r;
    r.first.first = pssstart;
    r.first.second = pssend;
    r.second = bestScore;
    return r;
}



subResult pos(const path& path1, const path& path2) {
    if (matricsType == "dtw") {
        return posDTW(path1, path2);
    } else if (matricsType == "FC"){

        return posFC(path1, path2);
    }
    else {
        return posWED(path1, path2);
    }
}
//}
