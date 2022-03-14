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

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;
    };
}

#endif //AGM3D_VALUE_H
