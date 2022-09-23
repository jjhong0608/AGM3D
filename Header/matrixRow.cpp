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

auto AGM::matrixRow::operator[](int i) -> double & {
    if (empty()) {
        auto args = matrixElement{};
        args.idx = i;
        args.value = ZEROVALUE;
        emplace_back(args);
        return front().value;
    }
    for (int j = 0; j < size(); ++j) {
        if (at(j).idx == i) {
            return at(j).value;
        }
        if (at(j).idx > i) {
            auto args = matrixElement{};
            args.idx = i;
            args.value = ZEROVALUE;
            emplace(begin() + j, args);
            return at(j).value;
        }
    }

    auto args = matrixElement{};
    args.idx = i;
    args.value = ZEROVALUE;
    emplace_back(args);
    return back().value;
}

auto AGM::matrixRow::operator+(const AGM::matrixRow &src) const -> AGM::matrixRow {
    auto row = *this;
    for (const auto &item: src) {
        row[item.idx] += item.value;
    }
    return row;
}

auto AGM::matrixRow::operator-(const AGM::matrixRow &src) const -> AGM::matrixRow {
    auto row = *this;
    for (const auto &item: src) {
        row[item.idx] -= item.value;
    }
    return row;
}

auto AGM::matrixRow::operator*(double d) const -> AGM::matrixRow {
    auto row = *this;
    for (auto &item: row) {
        item.value *= d;
    }
    return row;
}

auto AGM::matrixRow::operator+=(const AGM::matrixRow &src) -> AGM::matrixRow {
    for (const auto &item: src) {
        (*this)[item.idx] += item.value;
    }
    return *this;
}

auto AGM::matrixRow::operator-=(const AGM::matrixRow &src) -> AGM::matrixRow {
    for (const auto &item: src) {
        (*this)[item.idx] -= item.value;
    }
    return *this;
}
