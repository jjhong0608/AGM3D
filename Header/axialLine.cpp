//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "axialLine.h"

AGM::axialLine::axialLine() = default;

AGM::axialLine::axialLine(char mark) : mark(mark) {}

char AGM::axialLine::getMark() const {
    return mark;
}

void AGM::axialLine::setMark(char i) {
    axialLine::mark = i;
}

double &AGM::axialLine::operator[](int i) {
    return coordinate.at(i);
}

AGM::plane *AGM::axialLine::getPlane(int i) const {
    return plane[i];
}

void AGM::axialLine::setPlane(AGM::plane *pPlane, int i) {
    axialLine::plane[i] = pPlane;
}

double AGM::axialLine::operator-(AGM::axialLine &line) {
    auto subtract = [this, &line](int i0, int i1) -> double {
        if (!iszero(coordinate[i0] - line[i0])) {
            return coordinate[i0] - line[i0];
        } else {
            return coordinate[i1] - line[i1];
        }
    };
    if (mark == 'x') return subtract(2, 4);
    else if (mark == 'y') return subtract(0, 4);
    else if (mark == 'z') return subtract(0, 2);
    else printError("AGM::axialLine::operator-", "mark (which is %c) is wrong", mark);
    return ZEROVALUE;
}

AGM::axialLine::~axialLine() = default;
