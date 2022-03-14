//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "value.h"

AGM::value::value() : std::array<double, 11>{} {

}

AGM::value::~value() = default;

double &AGM::value::operator[](const std::string &string) {
    return at(valueMap[string]);
}

const double &AGM::value::operator[](const std::string &string) const {
    return at(valueMap[string]);
}
