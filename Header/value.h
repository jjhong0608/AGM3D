//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_VALUE_H
#define AGM3D_VALUE_H

#include "coordinate.h"

namespace AGM {
    class value : public std::array<double, 11> {

    public:
        value();

        virtual ~value();

        auto operator[](const std::string &string) -> double &;

        auto operator[](const std::string &string) const -> const double &;
    };
}

#endif //AGM3D_VALUE_H
