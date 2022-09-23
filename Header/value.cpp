//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "value.h"

AGM::value::value() : std::array<double, 11>{} {

}

AGM::value::~value() = default;

auto AGM::value::operator[](const std::string &string) -> double & {
    return at(valueMap[string]);
}

auto AGM::value::operator[](const std::string &string) const -> const double & {
    return at(valueMap[string]);
}
