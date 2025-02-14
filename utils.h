#ifndef TRAJCSIMILAR_UTILS_H
#define TRAJCSIMILAR_UTILS_H
#include <fstream>
#include <algorithm>
#include "distance.h"
#include "constant.h"
#include "iostream"
#include "map"
#include <sstream>
#include <constant.h>
//自己添加的内容

//namespace gnns
//{

using namespace std;
path str2path(const string &s);

extern string matricsType;
extern string pruningType;
extern string dataType;
//extern int maxLen;
//extern int minLen;
//自己添加的内容
vector<path> readfile_rsse(const string& filename);

vector<path> readfile(const string& filename);
void file2paths(const string& filename, const string& target, int bound);
vector<string> split(string str, string pattern);
void paths2file(const string& filename, vector<path>);
void ids2file(const string &filename, const map<int, pair<vector<double>, subResult>>& paths, int queryID);
void dataAnalize(const string& filename);
void dataFilter(const string &filename);
void dataLengthFilter(const string &filename, int minLen, int maxLen, int num);

//}

#endif //TRAJCSIMILAR_UTILS_H
