//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "coordinate.h"

AGM::coordinate::coordinate() : std::array<double, 3>{} {
}

AGM::coordinate::coordinate(double x, double y, double z) : std::array<double, 3>{x, y, z} {
}

double AGM::coordinate::norm() const {
    return std::sqrt(std::pow(at(0), 2) + std::pow(at(1), 2) + std::pow(at(2), 2));
}

AGM::coordinate AGM::coordinate::operator+(const AGM::coordinate &src) const {
    return {at(0) + src.at(0), at(1) + src.at(1), at(2) + src.at(2)};
}

AGM::coordinate AGM::coordinate::operator-(const AGM::coordinate &src) const {
    return {at(0) - src.at(0), at(1) - src.at(1), at(2) - src.at(2)};
}

AGM::coordinate AGM::coordinate::operator*(double d) const {
    return {at(0) * d, at(1) * d, at(2) * d};
}

bool AGM::coordinate::operator==(const AGM::coordinate &src) const {
    return isclose(at(0), src.at(0)) && isclose(at(1), src.at(1)) && isclose(at(2), src.at(2));
}

bool AGM::coordinate::operator!=(const AGM::coordinate &src) const {
    return !(*this == src);
}

AGM::coordinate::~coordinate() = default;
