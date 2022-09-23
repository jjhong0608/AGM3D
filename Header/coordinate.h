//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_COORDINATE_H
#define AGM3D_COORDINATE_H

#include "GreenfunctionReactionDiffusion.h"

namespace AGM {
    class coordinate : public std::array<double, 3> {
    public:
        coordinate();

        coordinate(double x, double y, double z);

        [[nodiscard]] auto norm() const -> double;

        virtual ~coordinate();

        auto operator+(const coordinate &src) const -> coordinate;

        auto operator-(const coordinate &src) const -> coordinate;

        auto operator*(double d) const -> coordinate;

        auto operator==(const coordinate &src) const -> bool;

        auto operator!=(const coordinate &src) const -> bool;

    };
}

#endif //AGM3D_COORDINATE_H
