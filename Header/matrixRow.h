//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_MATRIXROW_H
#define AGM3D_MATRIXROW_H

#include "plane.h"

namespace AGM {
    struct matrixElement {
        int idx{};
        double value{};
    };

    class matrixRow : public std::vector<matrixElement> {
    public:
        void remove(int i);

        double &operator[](int i);

        matrixRow operator+(const matrixRow &src) const;

        matrixRow operator-(const matrixRow &src) const;

        matrixRow operator*(double d) const;

        matrixRow operator+=(const matrixRow &src);

        matrixRow operator-=(const matrixRow &src);
    };
}


#endif //AGM3D_MATRIXROW_H
