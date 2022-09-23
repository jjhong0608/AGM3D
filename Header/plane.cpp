//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "plane.h"

#include <utility>

AGM::plane::plane() = default;

AGM::plane::plane(std::string mark) : mark(std::move(mark)) {}

AGM::plane::plane(std::string mark, double coordinate) : mark(std::move(mark)), coordinate(coordinate) {}

auto AGM::plane::getMark() const -> const std::string & {
    return mark;
}

void AGM::plane::setMark(const std::string &string) {
    plane::mark = string;
}

auto AGM::plane::getCoordinate() const -> double {
    return coordinate;
}

void AGM::plane::setCoordinate(double d) {
    plane::coordinate = d;
}

void AGM::plane::setPlanePointer() {
    for (auto &item: *this) {
        item->setPlane(this, 0);
    }
}

AGM::plane::~plane() = default;
