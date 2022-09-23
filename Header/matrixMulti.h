//
// Created by NIMS-JUNHONG on 2022/03/03.
//

#ifndef AGM3D_MATRIXMULTI_H
#define AGM3D_MATRIXMULTI_H


#include "matrix.h"

namespace AGM {
    template<typename pt>
    class matrixMulti : public matrix<pt> {
    public:
        matrixMulti();

        matrixMulti(std::vector<pt> *pts, std::vector<pt> *pts0, std::vector<pt> *pts1);

        virtual ~matrixMulti();

        void calculateMatrix() override;

        void factorizeMatrix() override;

    private:
        std::vector<pt> *pts0{}, *pts1{};
    };

}

#endif //AGM3D_MATRIXMULTI_H
