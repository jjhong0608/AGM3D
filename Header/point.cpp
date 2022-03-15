//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include <functional>
#include "point.h"

int AGM::point::nPts;
std::vector<AGM::axialLine> *AGM::point::xline;
std::vector<AGM::axialLine> *AGM::point::yline;
std::vector<AGM::axialLine> *AGM::point::zline;
std::array<std::vector<AGM::plane>, 2> *AGM::point::xyplane;
std::array<std::vector<AGM::plane>, 2> *AGM::point::yzplane;
std::array<std::vector<AGM::plane>, 2> *AGM::point::xzplane;

AGM::point::point() = default;

AGM::point::point(int idx) : idx(idx) {}

AGM::point::point(const AGM::coordinate &xyz) : xyz(xyz) {}

AGM::point::point(int idx, const AGM::coordinate &xyz) : idx(idx), xyz(xyz) {}

AGM::point::point(const AGM::coordinate &xyz, double mp) : xyz(xyz), mp(mp) {}

AGM::point::point(int idx, const AGM::coordinate &xyz, double mp) : idx(idx), xyz(xyz), mp(mp) {}

int AGM::point::getIdx() const {
    return idx;
}

void AGM::point::setIdx(int i) {
    point::idx = i;
}

const AGM::coordinate &AGM::point::getXyz() const {
    return xyz;
}

void AGM::point::setXyz(const AGM::coordinate &coordinate) {
    point::xyz = coordinate;
}

const AGM::coordinate &AGM::point::getNormal() const {
    return normal;
}

void AGM::point::setNormal(const AGM::coordinate &coordinate) {
    point::normal = coordinate;
}

double AGM::point::getMp() const {
    return mp;
}

void AGM::point::setMp(double d) {
    point::mp = d;
}

char AGM::point::getCondition() const {
    return condition;
}

void AGM::point::setCondition(char i) {
    point::condition = i;
}

const std::array<AGM::point *, 30> &AGM::point::getElement() const {
    return element;
}

void AGM::point::setElement(const std::array<AGM::point *, 30> &array) {
    point::element = array;
}

const std::array<AGM::point *, 30> &AGM::point::getElement1() const {
    return element1;
}

void AGM::point::setElement1(const std::array<AGM::point *, 30> &array) {
    point::element1 = array;
}

const AGM::value &AGM::point::getValue() const {
    return values;
}

void AGM::point::setValue(const AGM::value &value) {
    point::values = value;
}

const std::array<AGM::matrixRow, 3> &AGM::point::getSolMatrixRow() const {
    return solMatrixRow;
}

void AGM::point::setSolMatrixRow(const std::array<AGM::matrixRow, 3> &row) {
    point::solMatrixRow = row;
}

const std::array<AGM::matrixRow, 3> &AGM::point::getDeriMatrixRow() const {
    return deriMatrixRow;
}

void AGM::point::setDeriMatrixRow(const std::array<AGM::matrixRow, 3> &row) {
    point::deriMatrixRow = row;
}

const std::array<double, 3> &AGM::point::getRb() const {
    return rb;
}

void AGM::point::setRb(const std::array<double, 3> &array) {
    point::rb = array;
}

const std::array<AGM::axialLine *, 3> &AGM::point::getAxialLine() const {
    return aline;
}

AGM::axialLine *&AGM::point::getAxialLine(char i) {
    if (i == 'x') return aline[0];
    if (i == 'y') return aline[1];
    if (i == 'z') return aline[2];
    printError("AGM::axialLine *&AGM::point::getAxialLine", "index (which = %c) is wrong", i);
    return aline[0];
}

void AGM::point::setAxialLine(const std::array<AGM::axialLine *, 3> &array) {
    point::aline = array;
}

void AGM::point::setAxialLine(AGM::axialLine *line, char i) {
    if (i == 'x') aline[0] = line;
    if (i == 'y') aline[1] = line;
    if (i == 'z') aline[2] = line;
}

int AGM::point::getNPts() {
    return nPts;
}

void AGM::point::setNPts(int i) {
    point::nPts = i;
}

std::vector<AGM::axialLine> *&AGM::point::getAxialLines(char i) {
    if (i == 'x') return xline;
    if (i == 'y') return yline;
    if (i == 'z') return zline;
    printError("std::vector<AGM::axialLine> *&AGM::point::getAxialLines", "input character (which is %c) is wrong", i);
    return xline;
}

void AGM::point::setAxialLines(std::vector<AGM::axialLine> *line, char i) {
    if (i == 'x') xline = line;
    if (i == 'y') yline = line;
    if (i == 'z') zline = line;
}

std::array<std::vector<AGM::plane>, 2> *&AGM::point::getPlane(std::string &str) {
    if (str == "xy") return xyplane;
    if (str == "yz") return yzplane;
    if (str == "xz") return xzplane;
    printError("std::array<std::vector<AGM::plane>, 2> *&AGM::point::getPlane", "input string (which is %s) is wrong",
               str.c_str());
    return xyplane;
}

void AGM::point::setPlane(std::array<std::vector<AGM::plane>, 2> *pln, const std::string &str) {
    if (str == "xy") xyplane = pln;
    if (str == "yz") yzplane = pln;
    if (str == "xz") xzplane = pln;
}

double &AGM::point::operator[](int i) {
    return xyz[i];
}

const double &AGM::point::operator[](int i) const {
    return xyz[i];
}

AGM::point *&AGM::point::operator[](AGM::STENCIL stencil) {
    return element[stencil];
}

double &AGM::point::operator[](const std::string &string) {
    return values[string];
}

const double &AGM::point::operator[](const std::string &string) const {
    return values[string];
}

double AGM::point::operator-(const AGM::point &src) {
    return (xyz - src.xyz).norm();
}

AGM::point &AGM::point::operator=(const AGM::point &src) {
    if (this != &src) {
        idx = src.idx;
        xyz = src.xyz;
        normal = src.normal;
        mp = src.mp;
        condition = src.condition;
        element = src.element;
        element1 = src.element1;
        values = src.values;
        solMatrixRow = src.solMatrixRow;
        deriMatrixRow = src.deriMatrixRow;
        rb = src.rb;
        dv = src.dv;
        aline = src.aline;
    }
    return *this;
}

void AGM::point::findStencil() {
    switch (condition) {
        case 'D':
        case 'N':
            findStencilBoundary();
            break;
        default:
            break;
    }
}

void AGM::point::findStencilBoundary() {
    auto findLeftPlane = [this](char lineIdx, int planeIdx) -> plane * {
        return this->getAxialLine(lineIdx)->getPlane(planeIdx) - 1;
    };
    auto findRightPlane = [this](char lineIdx, int planeIdx) -> plane * {
        return this->getAxialLine(lineIdx)->getPlane(planeIdx) + 1;
    };
    auto isContainLine = [](axialLine *aln, int lineIdx, double pt) -> bool {
        return isnegativezero((*aln)[lineIdx] - pt) && ispositivezero((*aln)[lineIdx + 1] - pt);
    };
    auto findLeftLine = [this, &isContainLine](plane *pln, char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        std::copy_if(pln->begin(), pln->end(), std::back_inserter(vec),
                     [&isContainLine, &lineIdx, &pt](axialLine *aln) -> bool {
                         return isContainLine(aln, 2 * lineIdx, pt);
                     });
        auto pline = std::find_if(vec.begin(), vec.end(), [this, &lineChar](axialLine *&aln) -> bool {
            return *std::next(&aln) == this->getAxialLine(lineChar);
            // return *(&aln + 1) == this->getAxialLine(lineChar);
        });
        auto vec0 = std::vector<axialLine *>{};
        std::copy_if(pln->begin(), pln->end(), std::back_inserter(vec0), [this, &lineChar](axialLine *aln) -> bool {
            return ispositive(*this->getAxialLine(lineChar) - *aln);
        });
//        auto pline0 = std::find_if(vec0.begin(), vec0.end(), [this, &lineChar](axialLine *&aln) -> bool {
//            return *std::next(&aln) == this->getAxialLine(lineChar);
//            // return *(&aln + 1) == this->getAxialLine(lineChar);
//        });
        auto pline0 = vec0.empty() ? vec0.end() : std::prev(vec0.end());
        if (pline == vec.end() || pline0 == vec0.end() || !iszero(**pline - **pline0)) return nullptr;
        return *pline;
    };
    auto findRightLine = [this, &isContainLine](plane *pln, char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        std::copy_if(pln->begin(), pln->end(), std::back_inserter(vec),
                     [&isContainLine, &lineIdx, &pt](axialLine *aln) -> bool {
                         return isContainLine(aln, 2 * lineIdx, pt);
                     });
        auto pline = std::find_if(vec.begin(), vec.end(), [this, &lineChar](axialLine *&aln) -> bool {
            return *std::prev(&aln) == this->getAxialLine(lineChar);
            // return *(&aln - 1) == this->getAxialLine(lineChar);
        });
        auto vec0 = std::vector<axialLine *>{};
        std::copy_if(pln->begin(), pln->end(), std::back_inserter(vec0), [this, &lineChar](axialLine *aln) -> bool {
            return isnegative(*this->getAxialLine(lineChar) - *aln);
        });
//        auto pline0 = std::find_if(vec0.begin(), vec0.end(), [this, &lineChar](axialLine *&aln) -> bool {
//            return *std::prev(&aln) == this->getAxialLine(lineChar);
//            // return *(&aln - 1) == this->getAxialLine(lineChar);
//        });
        auto pline0 = vec0.empty() ? vec0.end() : vec0.begin();
        if (pline == vec.end() || pline0 == vec0.end() || !iszero(**pline - **pline0)) return nullptr;
        return *pline;
    };
    auto findLeftPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return isnegativezero((*pt)[ptIdx] - xyz[ptIdx]);
        });
        return vec.back();
    };
    auto findRightPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return ispositivezero((*pt)[ptIdx] - xyz[ptIdx]);
        });
        return vec.front();
    };

    auto assignStencil = [&](STENCIL stencil, STENCIL stencil0, STENCIL stencil1, auto func, char lineChar, int lineIdx,
                             int planeIdx, double pt, double sign, double &n) -> void {
        if (element[stencil] == this && ispositive(sign * n)) return;
        auto alin = func(this->getAxialLine(lineChar)->getPlane(planeIdx), lineChar, lineIdx, pt);
        if (alin) {
            auto leftPt = findLeftPt(alin, lineIdx);
            auto rightPt = findRightPt(alin, lineIdx);
            if (leftPt && rightPt && leftPt == rightPt) {
                element[stencil] = leftPt;
                leftPt = nullptr;
                rightPt = nullptr;
                return;
            }
            if (rightPt) {
                element[stencil0] = rightPt;
                element[stencil] = nullptr;
            }
            if (leftPt) {
                element[stencil1] = leftPt;
                element[stencil] = nullptr;
            }

            if (!(rightPt || leftPt)) {
                printError("assignStencil in findStencil", "rightPt & leftPt do not exist");
            }
        } else {
            n = ZEROVALUE;
            double norm{std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2) + std::pow(normal[2], 2))};
            for (int i = 0; i < 3; ++i) {
                normal[i] /= norm;
            }
            element[stencil] = this;
        }
    };

    if (aline[1]) {
        if (!element[R]) assignStencil(R, RF, RB, findRightLine, 'y', 1, 0, xyz[1], UNITVALUE, normal[0]);
        if (!element[L]) assignStencil(L, LF, LB, findLeftLine, 'y', 1, 0, xyz[1], -UNITVALUE, normal[0]);
    } else if (aline[2]) {
        if (!element[R]) assignStencil(R, RU, RD, findRightLine, 'z', 2, 1, xyz[2], UNITVALUE, normal[0]);
        if (!element[L]) assignStencil(L, LU, LD, findLeftLine, 'z', 2, 1, xyz[2], -UNITVALUE, normal[0]);
    }

    if (aline[0]) {
        if (!element[F]) assignStencil(F, FR, FL, findRightLine, 'x', 0, 0, xyz[0], UNITVALUE, normal[1]);
        if (!element[B]) assignStencil(B, BR, BL, findLeftLine, 'x', 0, 0, xyz[0], -UNITVALUE, normal[1]);
    } else if (aline[2]) {
        if (!element[F]) assignStencil(F, FU, FD, findRightLine, 'z', 2, 0, xyz[2], UNITVALUE, normal[1]);
        if (!element[B]) assignStencil(B, BU, BD, findLeftLine, 'z', 2, 0, xyz[2], -UNITVALUE, normal[1]);
    }

    if (aline[0]) {
        if (!element[U]) assignStencil(U, UR, UL, findRightLine, 'x', 0, 1, xyz[0], UNITVALUE, normal[2]);
        if (!element[D]) assignStencil(D, DR, DL, findLeftLine, 'x', 0, 1, xyz[0], -UNITVALUE, normal[2]);
    } else if (aline[1]) {
        if (!element[U]) assignStencil(U, UF, UB, findRightLine, 'y', 1, 1, xyz[1], UNITVALUE, normal[2]);
        if (!element[D]) assignStencil(D, DF, DB, findLeftLine, 'y', 1, 1, xyz[1], -UNITVALUE, normal[2]);
    }
}

void AGM::point::calculateRepresentationFormula() {
    switch (condition) {
        case 'C':
            calculateRepresentationFormulaCross();
            break;
        case 'D':
            calculateRepresentationFormulaDirichlet();
            break;
        case 'N':
            calculateRepresentationFormulaNeumann();
            break;
        default:
            printError("AGM::point::calculateRepresentationFormula()", "boundary condition (which is %c) is wrong",
                       condition);
    }
}

void AGM::point::calculateRepresentationFormulaCross() {
    double xm = element[L]->getXyz()[0];
    double xb = getXyz()[0];
    double xp = element[R]->getXyz()[0];
    double ym = element[B]->getXyz()[1];
    double yb = getXyz()[1];
    double yp = element[F]->getXyz()[1];
    double zm = element[D]->getXyz()[2];
    double zb = getXyz()[2];
    double zp = element[U]->getXyz()[2];

    auto gFuncX{Greenfunction(xm, xb, xp, mp, mp)};
    auto gFuncY{Greenfunction(ym, yb, yp, mp, mp)};
    auto gFuncZ{Greenfunction(zm, zb, zp, mp, mp)};

    std::array<matrixRow, 3> row{};

    row[0][element[L]->getIdx()] = mp * gFuncX.green_function_t(xm);
    row[0][getIdx()] = -UNITVALUE;
    row[0][element[R]->getIdx()] = -mp * gFuncX.green_function_t(xp);

    row[0][element[L]->getIdx() + getNPts()] = gFuncX.green_integral('l');
    row[0][getIdx() + getNPts()] = gFuncX.green_integral('c');
    row[0][element[R]->getIdx() + getNPts()] = gFuncX.green_integral('r');

    row[0][element[L]->getIdx() + 2 * getNPts()] = gFuncX.green_integral('l');
    row[0][getIdx() + 2 * getNPts()] = gFuncX.green_integral('c');
    row[0][element[R]->getIdx() + 2 * getNPts()] = gFuncX.green_integral('r');

    rhsMatrixRow[0][element[L]->getIdx()] = gFuncX.green_integral('l');
    rhsMatrixRow[0][getIdx()] = gFuncX.green_integral('c');
    rhsMatrixRow[0][element[R]->getIdx()] = gFuncX.green_integral('r');

    partMatrixRow[0][element[L]->getIdx()] = gFuncX.green_integral_t('l');
    partMatrixRow[0][getIdx()] = gFuncX.green_integral_t('c');
    partMatrixRow[0][element[R]->getIdx()] = gFuncX.green_integral_t('r');

    row[1][element[B]->getIdx()] = mp * gFuncY.green_function_t(ym);
    row[1][getIdx()] = -UNITVALUE;
    row[1][element[F]->getIdx()] = -mp * gFuncY.green_function_t(yp);

    row[1][element[B]->getIdx() + getNPts()] = -gFuncY.green_integral('l');
    row[1][getIdx() + getNPts()] = -gFuncY.green_integral('c');
    row[1][element[F]->getIdx() + getNPts()] = -gFuncY.green_integral('r');

    rhsMatrixRow[1][element[B]->getIdx() + getNPts()] = gFuncY.green_integral('l');
    rhsMatrixRow[1][getIdx() + getNPts()] = gFuncY.green_integral('c');
    rhsMatrixRow[1][element[F]->getIdx() + getNPts()] = gFuncY.green_integral('r');

    partMatrixRow[1][element[B]->getIdx() + getNPts()] = gFuncY.green_integral_t('l');
    partMatrixRow[1][getIdx() + getNPts()] = gFuncY.green_integral_t('c');
    partMatrixRow[1][element[F]->getIdx() + getNPts()] = gFuncY.green_integral_t('r');

    row[2][element[D]->getIdx()] = mp * gFuncZ.green_function_t(zm);
    row[2][getIdx()] = -UNITVALUE;
    row[2][element[U]->getIdx()] = -mp * gFuncZ.green_function_t(zp);

    row[2][element[D]->getIdx() + 2 * getNPts()] = -gFuncZ.green_integral('l');
    row[2][getIdx() + 2 * getNPts()] = -gFuncZ.green_integral('c');
    row[2][element[U]->getIdx() + 2 * getNPts()] = -gFuncZ.green_integral('r');

    rhsMatrixRow[2][element[D]->getIdx() + 2 * getNPts()] = gFuncZ.green_integral('l');
    rhsMatrixRow[2][getIdx() + 2 * getNPts()] = gFuncZ.green_integral('c');
    rhsMatrixRow[2][element[U]->getIdx() + 2 * getNPts()] = gFuncZ.green_integral('r');

    partMatrixRow[2][element[D]->getIdx() + 2 * getNPts()] = gFuncZ.green_integral_t('l');
    partMatrixRow[2][getIdx() + 2 * getNPts()] = gFuncZ.green_integral_t('c');
    partMatrixRow[2][element[U]->getIdx() + 2 * getNPts()] = gFuncZ.green_integral_t('r');

    solMatrixRow[0] = row[0] + row[1] + row[2];
    solMatrixRow[1] = row[0] - row[1] + row[2];
    solMatrixRow[2] = row[0] + row[1] - row[2];
}

void AGM::point::calculateRepresentationFormulaDirichlet() {
    solMatrixRow[0][getIdx()] = UNITVALUE;
    approximatePhiAndPsiAtBoundary(1);
}

void AGM::point::calculateRepresentationFormulaNeumann() {
    std::array<matrixRow, 3> row{};
    row[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial('x', 0)
                               : calculateRepresentationFormulaNeumannOffAxial('x', 0);
    row[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial('y', 1)
                               : calculateRepresentationFormulaNeumannOffAxial('y', 1);
    row[2] = getAxialLine('z') ? calculateRepresentationFormulaNeumannOnAxial('z', 2)
                               : calculateRepresentationFormulaNeumannOffAxial('z', 2);
    for (int i = 0; i < 3; ++i) {
        if (!row[i].empty() && !iszero(normal[i])) {
            while (row[i].back().idx >= 6 * getNPts()) {
                partMatrixRow[i][row[i].back().idx - 6 * getNPts()] = row[i].back().value * normal[i];
                row[i].pop_back();
            }
            while (row[i].back().idx >= 3 * getNPts()) {
                rhsMatrixRow[i][row[i].back().idx - 3 * getNPts()] = row[i].back().value * normal[i];
                row[i].pop_back();
            }
            solMatrixRow[0] += row[i] * normal[i];
        }
    }
    approximatePhiAndPsiAtBoundary(1);
}

AGM::matrixRow AGM::point::calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) {
    auto Error = []() -> double {
        printError("AGM::point::calculateRepresentationFormulaNeumannOnAxial", "Null pointer");
        return ZEROVALUE;
    };
    point *ptc = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(1) :
                 getAxialLine(axis)->back()->getIdx() == getIdx() ? *std::prev(getAxialLine(axis)->end() - 1) : nullptr;
    point *ptl =
            getAxialLine(axis)->front()->getIdx() == getIdx() ? this : getAxialLine(axis)->back()->getIdx() == getIdx()
                                                                       ? *std::prev(getAxialLine(axis)->end() - 2)
                                                                       : nullptr;
    point *ptr = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(2) :
                 getAxialLine(axis)->back()->getIdx() == getIdx() ? this : nullptr;
    std::string string =
            getAxialLine(axis)->front()->getIdx() == getIdx() ? "ND" : getAxialLine(axis)->back()->getIdx() == getIdx()
                                                                       ? "DN" : "";
    double tm = ptl ? ptl->getXyz()[axisInt] : Error();
    double tb = ptc ? ptc->getXyz()[axisInt] : Error();
    double tp = ptr ? ptr->getXyz()[axisInt] : Error();
    double signPhi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
    double signPsi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
    bool signPhi1 = axis != 'z';
    bool signPsi1 = axis != 'y';
    auto gFunc{Greenfunction(tm, tb, tp, mp, mp)};
    matrixRow row{};
    if (string == "ND") {
        row[ptl->getIdx()] = -mp * gFunc.green_function_ND(tm);
        row[ptc->getIdx()] = -UNITVALUE;
        row[ptr->getIdx()] = UNITVALUE;

        if (signPhi1) {
            row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('l');
            row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('c');
            row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('r');
        }

        if (signPsi1) {
            row[ptl->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_ND('l');
            row[ptc->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_ND('c');
            row[ptr->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_ND('r');
        }

        row[ptl->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_ND('l');
        row[ptc->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_ND('c');
        row[ptr->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_ND('r');

        row[ptl->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_ND('l') + gFunc.green_function_ND(tm);
        row[ptc->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_ND('c');
        row[ptr->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_ND('r');
    } else if (string == "DN") {
        row[ptl->getIdx()] = UNITVALUE;
        row[ptc->getIdx()] = -UNITVALUE;
        row[ptr->getIdx()] = mp * gFunc.green_function_DN(tp);

        if (signPhi1) {
            row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('l');
            row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('c');
            row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('r');
        }

        if (signPsi1) {
            row[ptl->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_DN('l');
            row[ptc->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_DN('c');
            row[ptr->getIdx() + 2 * getNPts()] = signPsi0 * gFunc.green_integral_DN('r');
        }

        row[ptl->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_DN('l');
        row[ptc->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_DN('c');
        row[ptr->getIdx() + (axisInt + 3) * getNPts()] = gFunc.green_integral_DN('r');

        row[ptl->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_DN('l');
        row[ptc->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_DN('c');
        row[ptr->getIdx() + (axisInt + 6) * getNPts()] = gFunc.green_integral_t_DN('r') - gFunc.green_function_DN(tp);
    }
    auto c = -row[getIdx()];
    for (auto &item: row) {
        item.value /= c;
    }
    row.remove(getIdx());
    return row;
}

AGM::matrixRow AGM::point::calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) {
    auto Error = []() -> double {
        printError("AGM::point::calculateRepresentationFormulaNeumannOffAxial", "aline error");
        return ZEROVALUE;
    };
    double tm{}, tb{}, tp{};
    if (axis == 'x') {
        tm = element[L] ? element[L]->getXyz()[0] : aline[1] ? element[LF]->getXyz()[0] : aline[2]
                                                                                          ? element[LU]->getXyz()[0]
                                                                                          : Error();
        tb = getXyz()[0];
        tp = element[R] ? element[R]->getXyz()[0] : aline[1] ? element[RF]->getXyz()[0] : aline[2]
                                                                                          ? element[RU]->getXyz()[0]
                                                                                          : Error();
    } else if (axis == 'y') {
        tm = element[B] ? element[B]->getXyz()[1] : aline[0] ? element[BR]->getXyz()[1] : aline[2]
                                                                                          ? element[BU]->getXyz()[1]
                                                                                          : Error();
        tb = getXyz()[1];
        tp = element[F] ? element[F]->getXyz()[1] : aline[0] ? element[FR]->getXyz()[1] : aline[2]
                                                                                          ? element[FU]->getXyz()[1]
                                                                                          : Error();
    } else if (axis == 'z') {
        tm = element[D] ? element[D]->getXyz()[2] : aline[0] ? element[DR]->getXyz()[2] : aline[1]
                                                                                          ? element[DF]->getXyz()[2]
                                                                                          : Error();
        tb = getXyz()[2];
        tp = element[U] ? element[U]->getXyz()[2] : aline[0] ? element[UR]->getXyz()[2] : aline[1]
                                                                                          ? element[UF]->getXyz()[2]
                                                                                          : Error();
    }
    char realAxis{};
    for (const auto &item: {'z', 'y', 'x'}) {
        if (getAxialLine(item)) realAxis = item;
    }
    double signPhi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
    double signPsi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
    bool signPhi1 = axis != 'z';
    bool signPsi1 = axis != 'y';
    auto gFunc{Greenfunction(tm, tb, tp, mp, mp)};
    auto approximateSol = [this, &realAxis](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
        double m{ptl->getXyz()[i]}, b{xyz[i]}, p{ptr->getXyz()[i]};
        auto func{Greenfunction(m, b, p, d, d)};
        auto mRow{matrixRow()};
        double signPhi0 = realAxis == 'x' ? UNITVALUE : -UNITVALUE;
        double signPsi0 = realAxis == 'x' ? UNITVALUE : -UNITVALUE;
        bool signPhi1 = realAxis != 'z';
        bool signPsi1 = realAxis != 'y';

        mRow[ptl->getIdx()] = d * func.green_function_t(m);
        mRow[ptr->getIdx()] = -d * func.green_function_t(p);

        if (signPhi1) {
            mRow[ptl->getIdx() + getNPts()] = signPhi0 * func.green_integral('L');
            mRow[ptr->getIdx() + getNPts()] = signPhi0 * func.green_integral('R');
        }

        if (signPsi1) {
            mRow[ptl->getIdx() + 2 * getNPts()] = signPsi0 * func.green_integral('L');
            mRow[ptr->getIdx() + 2 * getNPts()] = signPsi0 * func.green_integral('R');
        }

        mRow[ptl->getIdx() + (i + 3) * getNPts()] = func.green_integral('L');
        mRow[ptr->getIdx() + (i + 3) * getNPts()] = func.green_integral('R');

        mRow[ptl->getIdx() + (i + 6) * getNPts()] = func.green_integral_t('L');
        mRow[ptr->getIdx() + (i + 6) * getNPts()] = func.green_integral_t('R');

        return mRow * coefficient;
    };
    auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int i, int plus) -> matrixRow {
        double m{ptl->getXyz()[i]}, b{xyz[i]}, p{ptr->getXyz()[i]};
        auto mRow{matrixRow()};
        mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
        mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

        return mRow * coefficient;
    };
    matrixRow row{};
    auto assignMatrix = [&row, &approximateSol, &linearApproximation, &signPhi0, &signPsi0, &signPhi1, &signPsi1](
            point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, int i, int i0,
            char c) -> void {
        if (pt) {
            row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
            if (signPhi1) row[pt->getIdx() + getNPts()] += signPhi0 * func->green_integral_tau(c);
            if (signPsi1) row[pt->getIdx() + 2 * getNPts()] += signPsi0 * func->green_integral_tau(c);
            row[pt->getIdx() + (i0 + 3) * getNPts()] += func->green_integral_tau(c);
            row[pt->getIdx() + (i0 + 6) * getNPts()] += func->green_integral_ttau(c);
        } else {
            row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
            if (signPhi1) row += linearApproximation(ptr, ptl, signPhi0 * func->green_integral_tau(c), i, 1);
            if (signPsi1) row += linearApproximation(ptr, ptl, signPsi0 * func->green_integral_tau(c), i, 2);
            row += linearApproximation(ptr, ptl, func->green_integral_tau(c), i, i0 + 3);
            row += linearApproximation(ptr, ptl, func->green_integral_ttau(c), i, i0 + 6);
        }
    };
    if (signPhi1) row[getIdx() + getNPts()] += signPhi0 * gFunc.green_integral_tau('c');
    if (signPsi1) row[getIdx() + 2 * getNPts()] += signPsi0 * gFunc.green_integral_tau('c');
    row[getIdx() + (axisInt + 3) * getNPts()] += gFunc.green_integral_tau('c');
    row[getIdx() + (axisInt + 6) * getNPts()] += gFunc.green_integral_ttau('c') + UNITVALUE / mp;

    if (axis == 'x') {
        if (aline[1]) {
            assignMatrix(element[R], element[RF], element[RB], -mp, &gFunc, tp, 1, 0, 'r');
            assignMatrix(element[L], element[LF], element[LB], mp, &gFunc, tm, 1, 0, 'l');
        } else if (aline[2]) {
            assignMatrix(element[R], element[RU], element[RD], -mp, &gFunc, tp, 2, 0, 'r');
            assignMatrix(element[L], element[LU], element[LD], mp, &gFunc, tm, 2, 0, 'l');
        }
    } else if (axis == 'y') {
        if (aline[0]) {
            assignMatrix(element[F], element[FR], element[FL], -mp, &gFunc, tp, 0, 1, 'r');
            assignMatrix(element[B], element[BR], element[BL], mp, &gFunc, tm, 0, 1, 'l');
        } else if (aline[2]) {
            assignMatrix(element[F], element[FU], element[FD], -mp, &gFunc, tp, 2, 1, 'r');
            assignMatrix(element[B], element[BU], element[BD], mp, &gFunc, tm, 2, 1, 'l');
        }
    } else if (axis == 'z') {
        if (aline[0]) {
            assignMatrix(element[U], element[UR], element[UL], -mp, &gFunc, tp, 0, 2, 'r');
            assignMatrix(element[D], element[DR], element[DL], mp, &gFunc, tm, 0, 2, 'l');
        } else if (aline[1]) {
            assignMatrix(element[U], element[UF], element[UB], -mp, &gFunc, tp, 1, 2, 'r');
            assignMatrix(element[D], element[DF], element[DB], mp, &gFunc, tm, 1, 2, 'l');
        }
    }
    return row;
}

void AGM::point::approximatePhiAndPsiAtBoundary(int order) {
    auto findInnerPointOfBoundary = [this]() -> point * {
        for (const auto &item: {'x', 'y', 'z'}) {
            if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
                return getAxialLine(item)->at(2);
            }
            if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
                return *std::prev(getAxialLine(item)->end() - 2);
            }
        }
        printError("AGM::point::approximatePhiAndPsiAtBoundary", "findInnerPointOfBoundary");
        return nullptr;
    };
    auto findStencil = [this]() -> STENCIL {
        if (getAxialLine('x') && getAxialLine('x')->front()->getIdx() == getIdx()) return L;
        if (getAxialLine('x') && getAxialLine('x')->back()->getIdx() == getIdx()) return R;
        if (getAxialLine('y') && getAxialLine('y')->front()->getIdx() == getIdx()) return B;
        if (getAxialLine('y') && getAxialLine('y')->back()->getIdx() == getIdx()) return F;
        if (getAxialLine('z') && getAxialLine('z')->front()->getIdx() == getIdx()) return D;
        if (getAxialLine('z') && getAxialLine('z')->back()->getIdx() == getIdx()) return U;
        printError("AGM::point::approximatePhiAndPsiAtBoundary", "findStencil");
        return R;
    };
    auto pt{*findInnerPointOfBoundary()};
    auto stencil{findStencil()};
    double tm{}, tb{}, tp{};
    if (getAxialLine('x')) {
        tm = pt[L]->getXyz()[0];
        tb = pt.getXyz()[0];
        tp = pt[R]->getXyz()[0];
    } else if (getAxialLine('y')) {
        tm = pt[B]->getXyz()[1];
        tb = pt.getXyz()[1];
        tp = pt[F]->getXyz()[1];
    } else if (getAxialLine('z')) {
        tm = pt[D]->getXyz()[2];
        tb = pt.getXyz()[2];
        tp = pt[U]->getXyz()[2];
    }
    auto secondOrderExtrapolation = [this, &pt](int i, double m, double b, double p, STENCIL stencil1,
                                                STENCIL stencil2) -> void {
        double d{(p - m) * (p - b) * (b - m)}, t0{xyz[i]};
        auto firstTerm = [&t0, &d](double w0, double w1) -> double { return t0 * t0 * (w0 - w1) / d; };
        auto secondTerm = [&t0, &d](double w0, double w1) -> double { return t0 * (w0 * w0 - w1 * w1) / d; };
        auto thirdTerm = [&d](double w0, double w1) -> double { return w0 * w1 * (w0 - w1) / d; };
        auto secondOrder = [&firstTerm, &secondTerm, &thirdTerm](double w0, double w1) -> double {
            return firstTerm(w0, w1) - secondTerm(w0, w1) + thirdTerm(w0, w1);
        };
        solMatrixRow[1][pt[stencil1]->getIdx() + getNPts()] = secondOrder(p, b);
        solMatrixRow[1][getIdx() + getNPts()] = -secondOrder(p, m);
        solMatrixRow[1][pt[stencil2]->getIdx() + getNPts()] = secondOrder(b, m);
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;

        solMatrixRow[2][pt[stencil1]->getIdx() + 2 * getNPts()] = secondOrder(p, b);
        solMatrixRow[2][getIdx() + 2 * getNPts()] = -secondOrder(p, m);
        solMatrixRow[2][pt[stencil2]->getIdx() + 2 * getNPts()] = secondOrder(b, m);
        solMatrixRow[2][getIdx() + 2 * getNPts()] = -UNITVALUE;
    };
    auto firstOrderExtrapolation = [this, &pt](int i, double t1, double t2, STENCIL stencil0) -> void {
        double d{t2 - t1}, t0{xyz[i]};
        auto firstTerm = [&t0, &d]() -> double { return t0 / d; };
        auto secondTerm = [&d](double t) -> double { return t / d; };
        auto firstOrder = [&firstTerm, &secondTerm](double t) -> double { return -firstTerm() + secondTerm(t); };
        solMatrixRow[1][pt[stencil0]->getIdx() + getNPts()] = -firstOrder(t1);
        solMatrixRow[1][pt.getIdx() + getNPts()] = firstOrder(t2);
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;

        solMatrixRow[2][pt[stencil0]->getIdx() + 2 * getNPts()] = -firstOrder(t1);
        solMatrixRow[2][pt.getIdx() + 2 * getNPts()] = firstOrder(t2);
        solMatrixRow[2][getIdx() + 2 * getNPts()] = -UNITVALUE;
    };
    auto zeroOrderExtrapolation = [this, &pt](STENCIL stencil1) -> void {
        solMatrixRow[1][pt[stencil1]->getIdx() + getNPts()] = UNITVALUE;
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;

        solMatrixRow[2][pt[stencil1]->getIdx() + 2 * getNPts()] = UNITVALUE;
        solMatrixRow[2][getIdx() + 2 * getNPts()] = -UNITVALUE;
    };
    if (getAxialLine('x')) {
        if (order == 2) {
            secondOrderExtrapolation(0, tm, tb, tp, L, R);
        } else if (order == 1) {
            if (stencil == R) {
                firstOrderExtrapolation(0, tb, tp, stencil);
            } else if (stencil == L) {
                firstOrderExtrapolation(0, tb, tm, stencil);
            } else {
                printError("AGM::point::approximatePhiAndPsiAtBoundary", "stencil (which is %d) is wrong", stencil);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(stencil);
        } else {
            printError("AGM::point::approximatePhiAndPsiAtBoundary", "order (which is %d) is wrong", order);
        }
    } else if (getAxialLine('y')) {
        if (order == 2) {
            secondOrderExtrapolation(1, tm, tb, tp, B, F);
        } else if (order == 1) {
            if (stencil == F) {
                firstOrderExtrapolation(1, tb, tp, stencil);
            } else if (stencil == B) {
                firstOrderExtrapolation(1, tb, tm, stencil);
            } else {
                printError("AGM::point::approximatePhiAndPsiAtBoundary", "stencil (which is %d) is wrong", stencil);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(stencil);
        } else {
            printError("AGM::point::approximatePhiAndPsiAtBoundary", "order (which is %d) is wrong", order);
        }
    } else if (getAxialLine('z')) {
        if (order == 2) {
            secondOrderExtrapolation(2, tm, tb, tp, D, U);
        } else if (order == 1) {
            if (stencil == U) {
                firstOrderExtrapolation(2, tb, tp, stencil);
            } else if (stencil == D) {
                firstOrderExtrapolation(2, tb, tm, stencil);
            } else {
                printError("AGM::point::approximatePhiAndPsiAtBoundary", "stencil (which is %d) is wrong", stencil);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(stencil);
        } else {
            printError("AGM::point::approximatePhiAndPsiAtBoundary", "order (which is %d) is wrong", order);
        }
    } else {
        printError("AGM::point::approximatePhiAndPsiAtBoundary", "getAxialLine error");
    }
}

void AGM::point::updateRightHandSide(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                     const std::function<double(int)> &h) {
    switch (condition) {
        case 'C':
            updateRightHandSideCross(f, g, h);
            break;
        case 'D':
            updateRightHandSideDirichlet(f, g, h);
            break;
        case 'N':
            updateRightHandSideNeumann(f, g, h);
            break;
        default:
            printError("AGM::point::updateRightHandSide", "condition (which is %d) is wrong", condition);
    }
}

void AGM::point::updateRightHandSideCross(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                          const std::function<double(int)> &h) {
    rb[0] = rb[1] = rb[2] = ZEROVALUE;
    for (const auto &item: rhsMatrixRow[0]) {
        rb[0] -= item.value * f(item.idx);
        rb[1] -= item.value * f(item.idx);
        rb[2] -= item.value * f(item.idx);
    }
    for (const auto &item: rhsMatrixRow[1]) {
        rb[0] -= item.value * g(item.idx - getNPts());
        rb[1] += item.value * g(item.idx - getNPts());
        rb[2] -= item.value * g(item.idx - getNPts());
    }
    for (const auto &item: rhsMatrixRow[2]) {
        rb[0] -= item.value * h(item.idx - 2 * getNPts());
        rb[1] -= item.value * h(item.idx - 2 * getNPts());
        rb[2] += item.value * h(item.idx - 2 * getNPts());
    }
}

void AGM::point::updateRightHandSideDirichlet(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                              const std::function<double(int)> &h) {
    rb[0] = values["bdv"];
    rb[1] = rb[2] = ZEROVALUE;
}

void AGM::point::updateRightHandSideNeumann(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                            const std::function<double(int)> &h) {
    rb[0] = values["bdv"];
    rb[1] = rb[2] = ZEROVALUE;
    for (const auto &item: rhsMatrixRow) {
        for (const auto &item0: item) {
            if (item0.idx < getNPts()) {
                rb[0] -= item0.value * f(item0.idx);
            } else if (item0.idx < 2 * getNPts()) {
                rb[0] -= item0.value * g(item0.idx - getNPts());
            } else if (item0.idx < 3 * getNPts()) {
                rb[0] -= item0.value * h(item0.idx - 2 * getNPts());
            }
        }
    }
}

void AGM::point::updateRightHandSidePart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                         const std::function<double(int)> &h) {
    switch (condition) {
        case 'C':
            updateRightHandSideCrossPart(f, g, h);
            break;
        case 'D':
            updateRightHandSideDirichletPart(f, g, h);
            break;
        case 'N':
            updateRightHandSideNeumannPart(f, g, h);
            break;
        default:
            printError("AGM::point::updateRightHandSidePart", "condition (which is %d) is wrong", condition);
    }
}

void AGM::point::updateRightHandSideCrossPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                              const std::function<double(int)> &h) {
    for (const auto &item: partMatrixRow[0]) {
        rb[0] -= item.value * f(item.idx);
        rb[1] -= item.value * f(item.idx);
        rb[2] -= item.value * f(item.idx);
    }
    for (const auto &item: partMatrixRow[1]) {
        rb[0] -= item.value * g(item.idx - getNPts());
        rb[1] += item.value * g(item.idx - getNPts());
        rb[2] -= item.value * g(item.idx - getNPts());
    }
    for (const auto &item: partMatrixRow[2]) {
        rb[0] -= item.value * h(item.idx - 2 * getNPts());
        rb[1] -= item.value * h(item.idx - 2 * getNPts());
        rb[2] += item.value * h(item.idx - 2 * getNPts());
    }
}

void
AGM::point::updateRightHandSideDirichletPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                             const std::function<double(int)> &h) {

}

void
AGM::point::updateRightHandSideNeumannPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                           const std::function<double(int)> &h) {
    for (const auto &item: partMatrixRow) {
        for (const auto &item0: item) {
            if (item0.idx < getNPts()) {
                rb[0] -= item0.value * f(item0.idx);
            } else if (item0.idx < 2 * getNPts()) {
                rb[0] -= item0.value * g(item0.idx - getNPts());
            } else if (item0.idx < 3 * getNPts()) {
                rb[0] -= item0.value * h(item0.idx - 2 * getNPts());
            }
        }
    }
}

void AGM::point::makeDerivatives() {
    switch (condition) {
        case 'C':
            makeDerivativesCross();
            break;
        case 'D':
        case 'N':
            makeDerivativesBoundary();
            break;
        default:
            printError("AGM::point::makeDerivatives", "condition (which is %c) is wrong", condition);
    }
}

void AGM::point::makeDerivativesCross() {
    double xm = element[L]->getXyz()[0];
    double xb = getXyz()[0];
    double xp = element[R]->getXyz()[0];
    double ym = element[B]->getXyz()[1];
    double yb = getXyz()[1];
    double yp = element[F]->getXyz()[1];
    double zm = element[D]->getXyz()[2];
    double zb = getXyz()[2];
    double zp = element[U]->getXyz()[2];

    auto gFuncX = Greenfunction(xm, xb, xp, mp, mp);
    auto gFuncY = Greenfunction(ym, yb, yp, mp, mp);
    auto gFuncZ = Greenfunction(zm, zb, zp, mp, mp);

    deriMatrixRow[0][element[L]->getIdx()] = mp * gFuncX.green_function_ttau(xm);
    deriMatrixRow[0][element[R]->getIdx()] = -mp * gFuncX.green_function_ttau(xp);

    deriMatrixRow[0][element[L]->getIdx() + getNPts()] = gFuncX.green_integral_tau('l');
    deriMatrixRow[0][getIdx() + getNPts()] = gFuncX.green_integral_tau('c');
    deriMatrixRow[0][element[R]->getIdx() + getNPts()] = gFuncX.green_integral_tau('r');

    deriMatrixRow[0][element[L]->getIdx() + 2 * getNPts()] = gFuncX.green_integral_tau('l');
    deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFuncX.green_integral_tau('c');
    deriMatrixRow[0][element[R]->getIdx() + 2 * getNPts()] = gFuncX.green_integral_tau('r');

    deriMatrixRow[0][element[L]->getIdx() + 3 * getNPts()] = gFuncX.green_integral_tau('l');
    deriMatrixRow[0][getIdx() + 3 * getNPts()] = gFuncX.green_integral_tau('c');
    deriMatrixRow[0][element[R]->getIdx() + 3 * getNPts()] = gFuncX.green_integral_tau('r');

    deriMatrixRow[0][element[L]->getIdx() + 6 * getNPts()] = gFuncX.green_integral_ttau('l');
    deriMatrixRow[0][getIdx() + 6 * getNPts()] = gFuncX.green_integral_ttau('c') + UNITVALUE / mp;
    deriMatrixRow[0][element[R]->getIdx() + 6 * getNPts()] = gFuncX.green_integral_ttau('r');

    deriMatrixRow[1][element[B]->getIdx()] = mp * gFuncY.green_function_ttau(ym);
    deriMatrixRow[1][element[F]->getIdx()] = -mp * gFuncY.green_function_ttau(yp);

    deriMatrixRow[1][element[B]->getIdx() + getNPts()] = -gFuncY.green_integral_tau('l');
    deriMatrixRow[1][getIdx() + getNPts()] = -gFuncY.green_integral_tau('c');
    deriMatrixRow[1][element[F]->getIdx() + getNPts()] = -gFuncY.green_integral_tau('r');

    deriMatrixRow[1][element[B]->getIdx() + 4 * getNPts()] = gFuncY.green_integral_tau('l');
    deriMatrixRow[1][getIdx() + 4 * getNPts()] = gFuncY.green_integral_tau('c');
    deriMatrixRow[1][element[F]->getIdx() + 4 * getNPts()] = gFuncY.green_integral_tau('r');

    deriMatrixRow[1][element[B]->getIdx() + 7 * getNPts()] = gFuncY.green_integral_ttau('l');
    deriMatrixRow[1][getIdx() + 7 * getNPts()] = gFuncY.green_integral_ttau('c') + UNITVALUE / mp;
    deriMatrixRow[1][element[F]->getIdx() + 7 * getNPts()] = gFuncY.green_integral_ttau('r');

    deriMatrixRow[2][element[D]->getIdx()] = mp * gFuncZ.green_function_ttau(zm);
    deriMatrixRow[2][element[U]->getIdx()] = -mp * gFuncZ.green_function_ttau(zp);

    deriMatrixRow[2][element[D]->getIdx() + 2 * getNPts()] = -gFuncZ.green_integral_tau('l');
    deriMatrixRow[2][getIdx() + 2 * getNPts()] = -gFuncZ.green_integral_tau('c');
    deriMatrixRow[2][element[U]->getIdx() + 2 * getNPts()] = -gFuncZ.green_integral_tau('r');

    deriMatrixRow[2][element[D]->getIdx() + 5 * getNPts()] = gFuncZ.green_integral_tau('l');
    deriMatrixRow[2][getIdx() + 5 * getNPts()] = gFuncZ.green_integral_tau('c');
    deriMatrixRow[2][element[U]->getIdx() + 5 * getNPts()] = gFuncZ.green_integral_tau('r');

    deriMatrixRow[2][element[D]->getIdx() + 8 * getNPts()] = gFuncZ.green_integral_ttau('l');
    deriMatrixRow[2][getIdx() + 8 * getNPts()] = gFuncZ.green_integral_ttau('c') + UNITVALUE / mp;
    deriMatrixRow[2][element[U]->getIdx() + 8 * getNPts()] = gFuncZ.green_integral_ttau('r');
}

void AGM::point::makeDerivativesBoundary() {
    deriMatrixRow[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial('x', 0)
                                         : calculateRepresentationFormulaNeumannOffAxial('x', 0);
    deriMatrixRow[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial('y', 1)
                                         : calculateRepresentationFormulaNeumannOffAxial('y', 1);
    deriMatrixRow[2] = getAxialLine('z') ? calculateRepresentationFormulaNeumannOnAxial('z', 2)
                                         : calculateRepresentationFormulaNeumannOffAxial('z', 2);
}

void AGM::point::calculateDerivatives(const std::vector<point> *points, const std::function<double(int)> &f,
                                      const std::function<double(int)> &g, const std::function<double(int)> &h,
                                      const std::function<double(int)> &fp, const std::function<double(int)> &gp,
                                      const std::function<double(int)> &hp) {
    auto assignDerivatives = [&](int i) -> double {
        double d{};
        for (const auto &item: deriMatrixRow[i]) {
            if (item.idx < getNPts()) {
                d += item.value * points->at(item.idx)["sol"];
            } else if (item.idx < 2 * getNPts()) {
                d += item.value * points->at(item.idx - getNPts())["phi"];
            } else if (item.idx < 3 * getNPts()) {
                d += item.value * points->at(item.idx - 2 * getNPts())["psi"];
            } else if (item.idx < 4 * getNPts()) {
                d += item.value * f(item.idx - 3 * getNPts());
            } else if (item.idx < 5 * getNPts()) {
                d += item.value * g(item.idx - 4 * getNPts());
            } else if (item.idx < 6 * getNPts()) {
                d += item.value * h(item.idx - 5 * getNPts());
            } else if (item.idx < 7 * getNPts()) {
                d += item.value * fp(item.idx - 6 * getNPts());
            } else if (item.idx < 8 * getNPts()) {
                d += item.value * gp(item.idx - 7 * getNPts());
            } else if (item.idx < 9 * getNPts()) {
                d += item.value * hp(item.idx - 8 * getNPts());
            } else {
                printError("AGM::point::calculateDerivatives", "item.idx (which is %d) is too large", item.idx);
            }
        }
        return d;
    };
    values["dx"] = assignDerivatives(0);
    values["dy"] = assignDerivatives(1);
    values["dz"] = assignDerivatives(2);
}

void AGM::point::approximateNaNDerivatives(const std::vector<AGM::point> *points) {
    auto findInnerPointOfBoundary = [this]() -> point * {
        for (const auto &item: {'x', 'y', 'z'}) {
            if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
                return getAxialLine(item)->at(1);
            }
            if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
                return *std::prev(getAxialLine(item)->end() - 1);
            }
        }
        printError("AGM::point::approximateNaNDerivatives",
                   "findInnerPointOfBoundary, condition (which is %c) may be not boundary", getCondition());
        return nullptr;
    };
    if (std::isnan(values["dx"])) values["dx"] = findInnerPointOfBoundary()->getValue()["dx"];
    if (std::isnan(values["dy"])) values["dy"] = findInnerPointOfBoundary()->getValue()["dy"];
    if (std::isnan(values["dz"])) values["dz"] = findInnerPointOfBoundary()->getValue()["dz"];
}

void AGM::point::calculateSecondDerivatives(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                            const std::function<double(int)> &h) {
    values["dxx"] = -(f(getIdx()) + values["phi"] + values["psi"]) / mp;
    values["dyy"] = -(g(getIdx()) - values["phi"]) / mp;
    values["dzz"] = -(h(getIdx()) - values["psi"]) / mp;
}

AGM::point::~point() = default;
