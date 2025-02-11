#ifndef TRAJCSIMILAR_MOSTSIMILAR_H
#define TRAJCSIMILAR_MOSTSIMILAR_H


#include <sys/time.h>
#include "istream"
#include "distance.h"
#include "utils.h"
#include "pss.h"
#include "exactS.h"
#include "similarTrajectory.h"


//自己添加的内容
//namespace gnns
//{

extern string matricsType;
extern string pruningType;
extern string gatherType;
extern string dataType;

extern int maxLen;
extern int minLen;
long nowTime();
//自己添加的内容
double execute_data_to_data(const string & algorithm, const path& p1, const path& p2);

subResult execute(const string & algorithm, const path& p1, const path& p2);
vector<path> mostSimilar(vector<path> paths, int query, const string& algorithm, int limit);
map<int, map<int, pair<vector<double>,subResult>>> multiSimilar(const vector<path>& paths, vector<int> querys, const string& algorithm, int limit);

//}
#endif //TRAJCSIMILAR_MOSTSIMILAR_H
