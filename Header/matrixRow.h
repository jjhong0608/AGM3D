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

        auto operator[](int i) -> double &;

        auto operator+(const matrixRow &src) const -> matrixRow;

        auto operator-(const matrixRow &src) const -> matrixRow;

        auto operator*(double d) const -> matrixRow;

        auto operator+=(const matrixRow &src) -> matrixRow;

        auto operator-=(const matrixRow &src) -> matrixRow;
    };
}


#endif //AGM3D_MATRIXROW_H
