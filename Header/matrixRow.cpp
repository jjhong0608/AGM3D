//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "matrixRow.h"

void AGM::matrixRow::remove(int i) {
    for (int j = 0; j < size(); ++j) {
        if (at(j).idx == i) {
            erase(begin() + j);
        }
    }
}

double &AGM::matrixRow::operator[](int i) {
    if (empty()) {
        auto args = matrixElement{};
        args.idx = i;
        args.value = ZEROVALUE;
        emplace_back(args);
        return front().value;
    } else {
        for (int j = 0; j < size(); ++j) {
            if (at(j).idx == i) {
                return at(j).value;
            } else if (at(j).idx > i) {
                auto args = matrixElement{};
                args.idx = i;
                args.value = ZEROVALUE;
                emplace(begin() + j, args);
                return at(j).value;
            }
        }
    }
    auto args = matrixElement{};
    args.idx = i;
    args.value = ZEROVALUE;
    emplace_back(args);
    return back().value;
}

AGM::matrixRow AGM::matrixRow::operator+(const AGM::matrixRow &src) const {
    auto row = *this;
    for (const auto &item: src) {
        row[item.idx] += item.value;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator-(const AGM::matrixRow &src) const {
    auto row = *this;
    for (const auto &item: src) {
        row[item.idx] -= item.value;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator*(double d) const {
    auto row = *this;
    for (auto &item: row) {
        item.value *= d;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator+=(const AGM::matrixRow &src) {
    for (const auto &item: src) {
        (*this)[item.idx] += item.value;
    }
    return *this;
}

AGM::matrixRow AGM::matrixRow::operator-=(const AGM::matrixRow &src) {
    for (const auto &item: src) {
        (*this)[item.idx] -= item.value;
    }
    return *this;
}
