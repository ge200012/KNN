#ifndef TRAJCSIMILAR_EXPERIMENT_H
#define TRAJCSIMILAR_EXPERIMENT_H

#include <sys/time.h>
#include "mostSimilar.h"

//自己添加的内容
//namespace gnns
//{

extern string matricsType;
extern string pruningType;
extern string gatherType;
extern string dataType;
extern long evaluateTime;
extern long algorithmTime;
extern int datasize;
extern int maxLen;
extern int minLen;

extern double gridSize;
extern double keyNum;
extern double fixRate;
extern double filterNum;


void showResult(const string& algorithm, const path& p1, const path& p2);
void findMostSimilar(const vector<path>& paths, int query, const string& algorithm, const string& targetFile, int limit=1);
void findMultiMostSimilar(const vector<path>& paths, vector<int>& querys, const string& algorithm, string targetFile, int limit);
void varyLength(int l, int r, string dataset);
void varyDataLength(int l, int r);
void varyDataSize(int size);
void varyQuerySize(int size);
void varyGridSize();
void varyRate();
void varyKeyNum();

//}
#endif //TRAJCSIMILAR_EXPERIMENT_H
