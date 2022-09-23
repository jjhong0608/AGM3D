//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "coordinate.h"

AGM::coordinate::coordinate() : std::array<double, 3>{} {
}

AGM::coordinate::coordinate(double x, double y, double z) : std::array<double, 3>{x, y, z} {
}

auto AGM::coordinate::norm() const -> double {
    return std::sqrt(std::pow(at(0), 2) + std::pow(at(1), 2) + std::pow(at(2), 2));
}

auto AGM::coordinate::operator+(const AGM::coordinate &src) const -> AGM::coordinate {
    return {at(0) + src.at(0), at(1) + src.at(1), at(2) + src.at(2)};
}

auto AGM::coordinate::operator-(const AGM::coordinate &src) const -> AGM::coordinate {
    return {at(0) - src.at(0), at(1) - src.at(1), at(2) - src.at(2)};
}

auto AGM::coordinate::operator*(double d) const -> AGM::coordinate {
    return {at(0) * d, at(1) * d, at(2) * d};
}

auto AGM::coordinate::operator==(const AGM::coordinate &src) const -> bool {
    return isclose(at(0), src.at(0)) && isclose(at(1), src.at(1)) && isclose(at(2), src.at(2));
}

auto AGM::coordinate::operator!=(const AGM::coordinate &src) const -> bool {
    return !(*this == src);
}

AGM::coordinate::~coordinate() = default;
