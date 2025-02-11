#ifndef GNNS_EVALUATION_H
#define GNNS_EVALUATION_H

#include "../util/matrix.h"
#include <cassert>
#include <set>

//namespace gnns
//{
    template<typename T>
    float compute_precision(const std::vector<std::vector<int>>& result, const std::vector<std::vector<int>>& ground_truth)
    {
        //assert(result.cols <= ground_truth.cols);
        //assert(result.rows == ground_truth.rows);

        int hit = 0;
        for(int i=0;i<result.size();++i)
        {
            for(int j=0;j<result[0].size();++j)
            {
                for(int k=0;k<result[0].size();++k)
                {
                    if(result[i][j] == ground_truth[i][k]) ++hit;
                }
            }
        }

        return 1.0 * hit / result.size() / result[0].size();
    }
//}


#endif //GNNS_EVALUATION_H
