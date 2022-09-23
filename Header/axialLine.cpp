//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "axialLine.h"

AGM::axialLine::axialLine() = default;

AGM::axialLine::axialLine(char mark) : mark(mark) {}

auto AGM::axialLine::getMark() const -> char {
    return mark;
}

void AGM::axialLine::setMark(char i) {
    axialLine::mark = i;
}

auto AGM::axialLine::operator[](int i) -> double & {
    return coordinate.at(i);
}

auto AGM::axialLine::getPlane(int i) const -> AGM::plane * {
    return plane.at(i);
}

void AGM::axialLine::setPlane(AGM::plane *pPlane, int i) {
    axialLine::plane.at(i) = pPlane;
}

auto AGM::axialLine::operator-(AGM::axialLine &line) -> double {
    auto subtract = [this, &line](int i0, int i1) -> double {
        if (!iszero(coordinate.at(i0) - line[i0])) {
            return coordinate.at(i0) - line[i0];
        }
        return coordinate.at(i1) - line[i1];
    };
    if (mark == 'x') {
        return subtract(2, 4);
    }
    if (mark == 'y') {
        return subtract(0, 4);
    }
    if (mark == 'z') {
        return subtract(0, 2);
    }
    printError("AGM::axialLine::operator-", "mark (which is %c) is wrong", mark);
    return ZEROVALUE;
}

AGM::axialLine::~axialLine() = default;
