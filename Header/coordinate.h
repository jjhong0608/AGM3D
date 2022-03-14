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

        [[nodiscard]] double norm() const;

        virtual ~coordinate();

        coordinate operator+(const coordinate &src) const;

        coordinate operator-(const coordinate &src) const;

        coordinate operator*(double d) const;

        bool operator==(const coordinate &src) const;

        bool operator!=(const coordinate &src) const;

    };
}

#endif //AGM3D_COORDINATE_H
