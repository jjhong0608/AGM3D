//
// Created by 조준홍 on 2022/02/03.
//

#include "pointHeat.h"

double AGM::pointHeat::time;
double AGM::pointHeat::delta;

double AGM::pointHeat::getTime() {
    return time;
}

void AGM::pointHeat::setTime(double d) {
    pointHeat::time = d;
}

double AGM::pointHeat::getDelta() {
    return delta;
}

void AGM::pointHeat::setDelta(double d) {
    pointHeat::delta = d;
}

void AGM::pointHeat::findStencil(const axialElement *axialElement1, std::vector<AGM::pointHeat> *vector) {
    int n{};
    for (const auto &item: *axialElement1) {
        if (item) {
            element.at(n) = &(vector->at(item->getIdx()));
        }
        ++n;
    }
}

void AGM::pointHeat::calculateRepresentationFormulaCross() {
    double xm = element[L]->getXyz()[0];
    double xb = getXyz()[0];
    double xp = element[R]->getXyz()[0];
    double ym = element[B]->getXyz()[1];
    double yb = getXyz()[1];
    double yp = element[F]->getXyz()[1];
    double zm = element[D]->getXyz()[2];
    double zb = getXyz()[2];
    double zp = element[U]->getXyz()[2];

    auto gFuncX{GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, 2.0 / 3.0 / delta)};
    auto gFuncY{GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, 2.0 / 3.0 / delta)};
    auto gFuncZ{GreenfunctionReactionDiffusion(zm, zb, zp, mp, mp, 2.0 / 3.0 / delta)};

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

AGM::matrixRow AGM::pointHeat::calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) {
    auto Error = []() -> double {
        printError("AGM::pointHeat::calculateRepresentationFormulaNeumannOnAxial", "Null pointer");
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
    auto gFunc{GreenfunctionReactionDiffusion(tm, tb, tp, mp, mp, 2.0 / 3.0 / delta)};
    matrixRow row{};
    if (string == "ND") {
        row[ptl->getIdx()] = -mp * gFunc.green_function_ND(tm);
        row[ptc->getIdx()] = -UNITVALUE;
        row[ptr->getIdx()] = -mp * gFunc.green_function_t_ND(tp);

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
        row[ptl->getIdx()] = mp * gFunc.green_function_t_DN(tm);
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

AGM::matrixRow AGM::pointHeat::calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) {
    auto Error = []() -> double {
        printError("AGM::pointHeat::calculateRepresentationFormulaNeumannOffAxial", "aline error");
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
    auto gFunc{GreenfunctionReactionDiffusion(tm, tb, tp, mp, mp, 2.0 / 3.0 / delta)};
    auto approximateSol = [this, &realAxis](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
        double m{ptl->getXyz()[i]}, b{xyz[i]}, p{ptr->getXyz()[i]};
        auto func{GreenfunctionReactionDiffusion(m, b, p, d, d, 2.0 / 3.0 / delta)};
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
            point *pt, point *ptr, point *ptl, double mp0, GreenfunctionReactionDiffusion *func, double d, int i,
            int i0, char c) -> void {
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

void AGM::pointHeat::makeDerivativesCross() {
    double xm = element[L]->getXyz()[0];
    double xb = getXyz()[0];
    double xp = element[R]->getXyz()[0];
    double ym = element[B]->getXyz()[1];
    double yb = getXyz()[1];
    double yp = element[F]->getXyz()[1];
    double zm = element[D]->getXyz()[2];
    double zb = getXyz()[2];
    double zp = element[U]->getXyz()[2];

    auto gFuncX = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, 2.0 / 3.0 / delta);
    auto gFuncY = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, 2.0 / 3.0 / delta);
    auto gFuncZ = GreenfunctionReactionDiffusion(zm, zb, zp, mp, mp, 2.0 / 3.0 / delta);

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

void AGM::pointHeat::calculateDerivatives(const std::vector<pointHeat> *points, const std::function<double(int)> &f,
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
                printError("AGM::pointHeat::calculateDerivatives", "item.idx (which is %d) is too large", item.idx);
            }
        }
        return d;
    };
    values["dx"] = assignDerivatives(0);
    values["dy"] = assignDerivatives(1);
    values["dz"] = assignDerivatives(2);
}

void AGM::pointHeat::approximateNaNDerivatives(const std::vector<AGM::pointHeat> *points) {
    auto findInnerPointOfBoundary = [this]() -> point * {
        for (const auto &item: {'x', 'y', 'z'}) {
            if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
                return getAxialLine(item)->at(1);
            }
            if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
                return *std::prev(getAxialLine(item)->end() - 1);
            }
        }
        printError("AGM::pointHeat::approximateNaNDerivatives", "findInnerPointOfBoundary");
        return nullptr;
    };
    if (std::isnan(values["dx"])) values["dx"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dx"];
    if (std::isnan(values["dy"])) values["dy"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dy"];
    if (std::isnan(values["dz"])) values["dz"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dz"];
}
