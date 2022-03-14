//
// Created by NIMS-JUNHONG on 2022/03/04.
//

#ifndef AGM3D_WRITEFILEMULTIPLE_H
#define AGM3D_WRITEFILEMULTIPLE_H


#include "writeFile.h"
#include "StdVector"

namespace AGM {
    template<typename T0, typename T1, typename T2, typename T3>
    class writeFileMultiple {
    public:
        writeFileMultiple();

        writeFileMultiple(std::vector<T0> *pts0, std::vector<T1> *pts1, std::vector<T2> *pts2, std::vector<T3> *pts3);

        void writeResult(const std::string &string);

        virtual ~writeFileMultiple();

    private:
        T0 *pt0{};
        T1 *pt1{};
        T2 *pt2{};
        std::vector<T0> *pts0{};
        std::vector<T1> *pts1{};
        std::vector<T2> *pts2{};
        std::vector<T3> *pts3{};
    };

}


#endif //AGM3D_WRITEFILEMULTIPLE_H
