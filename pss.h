#ifndef TRAJCSIMILAR_PSS_H
#define TRAJCSIMILAR_PSS_H

#include "distance.h"
#include "istream"

//自己添加的内容
//namespace gnns
//{

using namespace std;
subResult pss(const path& p1, const path& p2);
subResult pos(const path& p1, const path& p2);

extern string matricsType;
extern string pruningType;
extern string gatherType;
extern string dataType;

//}
#endif //TRAJCSIMILAR_PSS_H
