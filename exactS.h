#ifndef TRAJCSIMILAR_EXACTS_H
#define TRAJCSIMILAR_EXACTS_H

#include "distance.h"
//自己添加的内容
//namespace gnns 
//{

subResult exactS(const path& p1, const path& p2);

//自己添加的内容
double exactS_data_to_data(const path& p1, const path& p2);

void calScore(const path& p1, const path& p2, int start, int end);

extern string matricsType;
extern string pruningType;
extern string gatherType;
extern string dataType;
extern int minThreshold;

extern double AR;
extern double RR;
extern double MR;

//}
//extern int maxLen;
//extern int minLen;
#endif //TRAJCSIMILAR_EXACTS_H
