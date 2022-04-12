//
// Created by NIMS-JUNHONG on 2022/02/16.
//

#include "function.h"
#include "StdVector"

AGM::ellipticFunction::ellipticFunction() = default;

double AGM::ellipticFunction::u(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return std::sin(x) + std::cos(y) + std::sin(z);
}

double AGM::ellipticFunction::phi(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return HALFVALUE * (std::sin(x) - std::cos(y) + std::sin(z));
}

double AGM::ellipticFunction::psi(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return HALFVALUE * (std::sin(x) + std::cos(y) - std::sin(z));
}

double AGM::ellipticFunction::f(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return u(pt);
}

double AGM::ellipticFunction::ux(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return std::cos(x);
}

double AGM::ellipticFunction::uy(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return -std::sin(y);
}

double AGM::ellipticFunction::uz(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return std::cos(z);
}

void AGM::ellipticFunction::assignBoundaryValue(AGM::point &pt) {
    if (pt.getCondition() == 'D') {
        pt["bdv"] = u(pt);
    } else if (pt.getCondition() == 'N') {
        pt["bdv"] = ux(pt) * pt.getNormal()[0] + uy(pt) * pt.getNormal()[1] + uz(pt) * pt.getNormal()[2];
    }
    pt["rhs"] = f(pt);
}

AGM::ellipticFunction::~ellipticFunction() = default;

AGM::heatFunction::heatFunction() = default;

double AGM::heatFunction::initialTime() {
    return UNITVALUE;
}

double AGM::heatFunction::terminalTime() {
    return 1.25;
}

double AGM::heatFunction::deltaTime() {
    return 0.01;
}

double AGM::heatFunction::u(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    double r{std::sqrt(x * x + y * y + z * z)};
    return std::exp(-std::pow(r, 2) / (4 * t)) / std::pow(4 * M_PI * t, 1.5);
}

double AGM::heatFunction::phi(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return std::pow(M_PI, -1.5) * std::pow(t, -3.5) *
           (-0.010416666666666666 * std::pow(x, 2) + 0.020833333333333336 * std::pow(y, 2) -
            0.010416666666666666 * std::pow(z, 2)) *
           std::exp(-1.0 / 4.0 * (std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)) / t);
}

double AGM::heatFunction::psi(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return std::pow(M_PI, -1.5) * std::pow(t, -3.5) *
           (-0.010416666666666666 * std::pow(x, 2) - 0.010416666666666666 * std::pow(y, 2) +
            0.020833333333333336 * std::pow(z, 2)) *
           std::exp(-1.0 / 4.0 * (std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)) / t);
}

double AGM::heatFunction::f(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return ZEROVALUE;
}

double AGM::heatFunction::ux(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return -0.0625 * std::pow(M_PI, -1.5) * std::pow(t, -2.5) * x *
           exp(-1.0 / 4.0 * (std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)) / t);

}

double AGM::heatFunction::uy(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return -0.0625 * std::pow(M_PI, -1.5) * std::pow(t, -2.5) * y *
           exp(-1.0 / 4.0 * (std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)) / t);
}

double AGM::heatFunction::uz(double t, const AGM::pointHeat &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return -0.0625 * std::pow(M_PI, -1.5) * std::pow(t, -2.5) * z *
           exp(-1.0 / 4.0 * (std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)) / t);
}

void AGM::heatFunction::assignPreviousValue(AGM::value &value, AGM::pointHeat &pt) {
    value["sol"] = u(pointHeat::getTime(), pt);
    value["phi"] = phi(pointHeat::getTime(), pt);
    value["psi"] = psi(pointHeat::getTime(), pt);
    value["rhs"] = f(pointHeat::getTime(), pt);
    value["dx"] = ux(pointHeat::getTime(), pt);
    value["dy"] = uy(pointHeat::getTime(), pt);
    value["dz"] = uz(pointHeat::getTime(), pt);
}

void AGM::heatFunction::assignBoundaryValue(AGM::pointHeat &pt) {
    if (pt.getCondition() == 'D') {
        pt["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), pt);
    } else if (pt.getCondition() == 'N') {
        pt["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[0] +
                    uy(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[1] +
                    uz(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[2];
    }
    pt["rhs"] = f(pointHeat::getTime() + pointHeat::getDelta(), pt);
}

AGM::heatFunction::~heatFunction() = default;

AGM::NavierStokesFunction::NavierStokesFunction() = default;

double AGM::NavierStokesFunction::initialTime() {
    return ZEROVALUE;
}

double AGM::NavierStokesFunction::terminalTime() {
    return 3e1;
}

double AGM::NavierStokesFunction::deltaTime() {
    return 1e-2;
}

double AGM::NavierStokesFunction::u(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    if (isclose(z, UNITVALUE)) {
        return UNITVALUE;
    } else {
        return ZEROVALUE;
    }
}

double AGM::NavierStokesFunction::v(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return ZEROVALUE;
}

double AGM::NavierStokesFunction::w(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return ZEROVALUE;
}

double AGM::NavierStokesFunction::p(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::phiU(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::psiU(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::phiV(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::psiV(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::phiW(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::psiW(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::ux(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::uy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::uz(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::vx(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::vy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::vz(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::wx(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::wy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::wz(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::px(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::py(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::pz(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::f1(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::f2(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

double AGM::NavierStokesFunction::f3(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return 0;
}

void AGM::NavierStokesFunction::loadPreviousValue(const std::string &filename, std::vector<value> *pu,
                                                  std::vector<value> *pv, std::vector<value> *pw,
                                                  std::vector<value> *pp) {
    int idx{}, bc{};
    double x{}, y{}, z{};
    std::ifstream f(filename);
    for (int i = 0; i < point::getNPts(); ++i) {
        f >> idx >> x >> y >> z;
        if (idx >= pu->size()) {
            printError("AGM::NavierStokesFunction::loadPreviousValue",
                       "idx (which is %d) is greater(or equal) than size of the point (which is %d)", idx, pu->size());
        }
        f >> pu->at(idx)["sol"];
        f >> pv->at(idx)["sol"];
        f >> pw->at(idx)["sol"];
        f >> pp->at(idx)["sol"];
        f >> pu->at(idx)["dx"];
        f >> pu->at(idx)["dy"];
        f >> pu->at(idx)["dz"];
        f >> pv->at(idx)["dx"];
        f >> pv->at(idx)["dy"];
        f >> pv->at(idx)["dz"];
        f >> pw->at(idx)["dx"];
        f >> pw->at(idx)["dy"];
        f >> pw->at(idx)["dz"];
        f >> pu->at(idx)["phi"];
        f >> pu->at(idx)["psi"];
        f >> pv->at(idx)["phi"];
        f >> pv->at(idx)["psi"];
        f >> pw->at(idx)["phi"];
        f >> pw->at(idx)["psi"];
        f >> bc;
    }
    f.close();
}

void AGM::NavierStokesFunction::assignPreviousValue(AGM::value &pu, AGM::value &pv, AGM::value &pw, AGM::value &pp,
                                                    AGM::point &uvel, AGM::point &vvel, AGM::point &wvel,
                                                    AGM::point &pres) {
    pu["sol"] = u(pointHeat::getTime(), uvel);
    pv["sol"] = v(pointHeat::getTime(), vvel);
    pw["sol"] = w(pointHeat::getTime(), wvel);
    pu["phi"] = phiU(pointHeat::getTime(), uvel);
    pv["phi"] = phiV(pointHeat::getTime(), vvel);
    pw["phi"] = phiW(pointHeat::getTime(), wvel);
    pu["psi"] = psiU(pointHeat::getTime(), uvel);
    pv["psi"] = psiV(pointHeat::getTime(), vvel);
    pw["psi"] = psiW(pointHeat::getTime(), wvel);
    pu["rhs"] = f1(pointHeat::getTime(), uvel);
    pv["rhs"] = f2(pointHeat::getTime(), vvel);
    pw["rhs"] = f3(pointHeat::getTime(), wvel);
    pu["dx"] = ux(pointHeat::getTime(), uvel);
    pv["dx"] = vx(pointHeat::getTime(), vvel);
    pw["dx"] = wx(pointHeat::getTime(), wvel);
    pp["dx"] = px(pointHeat::getTime(), pres);
    pu["dy"] = uy(pointHeat::getTime(), uvel);
    pv["dy"] = vy(pointHeat::getTime(), vvel);
    pw["dy"] = wy(pointHeat::getTime(), wvel);
    pp["dy"] = py(pointHeat::getTime(), pres);
    pu["dz"] = uz(pointHeat::getTime(), uvel);
    pv["dz"] = vz(pointHeat::getTime(), vvel);
    pw["dz"] = wz(pointHeat::getTime(), wvel);
    pp["dz"] = pz(pointHeat::getTime(), pres);
}

void AGM::NavierStokesFunction::assignBoundaryValue(AGM::point &uvel, AGM::point &vvel, AGM::point &wvel) {
    if (uvel.getCondition() == 'D') {
        uvel["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    } else if (uvel.getCondition() == 'N') {
        uvel["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] +
                      uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1] +
                      uz(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[2];
    }
    uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    if (vvel.getCondition() == 'D') {
        vvel["bdv"] = v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else if (vvel.getCondition() == 'N') {
        vvel["bdv"] = vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] +
                      vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1] +
                      vz(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[2];
    }
    vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    if (wvel.getCondition() == 'D') {
        wvel["bdv"] = w(pointHeat::getTime() + pointHeat::getDelta(), wvel);
    } else if (wvel.getCondition() == 'N') {
        wvel["bdv"] = wx(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[0] +
                      wy(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[1] +
                      wz(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[2];
    }
    wvel["rhs"] = f3(pointHeat::getTime() + pointHeat::getDelta(), wvel);
}

AGM::NavierStokesFunction::~NavierStokesFunction() = default;
