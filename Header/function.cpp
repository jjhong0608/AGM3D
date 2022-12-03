//
// Created by NIMS-JUNHONG on 2022/02/16.
//

#include "function.h"
#include "StdVector"

AGM::ellipticFunction::ellipticFunction() = default;

auto AGM::ellipticFunction::u(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return std::sin(x) + std::cos(y) + std::sin(z);
}

auto AGM::ellipticFunction::phi(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return HALFVALUE * (std::sin(x) - std::cos(y) + std::sin(z));
}

auto AGM::ellipticFunction::psi(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return HALFVALUE * (std::sin(x) + std::cos(y) - std::sin(z));
}

auto AGM::ellipticFunction::f(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return u(pt);
}

auto AGM::ellipticFunction::ux(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return std::cos(x);
}

auto AGM::ellipticFunction::uy(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return -std::sin(y);
}

auto AGM::ellipticFunction::uz(const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
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

auto AGM::heatFunction::initialTime() -> double {
    return UNITVALUE;
}

auto AGM::heatFunction::terminalTime() -> double {
    const auto time{1.25};
    return time;
}

auto AGM::heatFunction::deltaTime() -> double {
    const auto time{1e-2};
    return time;
}

auto AGM::heatFunction::u(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    double r{std::sqrt(x * x + y * y + z * z)};
    const auto oneHalf{1.5e0};
    return std::exp(-std::pow(r, 2) / (4 * t)) / std::pow(4 * M_PI * t, oneHalf);
}

auto AGM::heatFunction::phi(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::heatFunction::psi(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::heatFunction::f(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]}, y{pt[1]}, z{pt[2]};
    return ZEROVALUE;
}

auto AGM::heatFunction::ux(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;

}

auto AGM::heatFunction::uy(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::heatFunction::uz(double t, const AGM::pointHeat &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
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

auto AGM::NavierStokesFunction::initialTime() -> double {
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::terminalTime() -> double {
    const auto time{5e1};
    return time;
}

auto AGM::NavierStokesFunction::deltaTime() -> double {
    const auto time{1e-3};
    return time;
}

auto AGM::NavierStokesFunction::writeTime() -> double {
    const auto time{1e-1};
    return time;
}

auto AGM::NavierStokesFunction::u(double t, const AGM::point &pt) -> double {
    const auto x{pt[0]};
    const auto y{pt[1]};
    const auto z{pt[2]};
    const auto r{std::sqrt(std::pow(y, 2) + std::pow(z, 2))};
    const auto radius{2e0};
    const auto maxVelocity{2e1 / (2 * radius)};
    if (isclose(pt.getNormal()[0], -UNITVALUE)) {
        return maxVelocity * (radius + r) * (radius - r);
    }
    return ZEROVALUE;

    /* cavity flow */
    if (isclose(z, UNITVALUE)) {
        return UNITVALUE;
    }
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::v(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::w(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::p(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::phiU(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::psiU(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::phiV(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::psiV(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::phiW(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::psiW(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::ux(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::uy(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::uz(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::vx(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::vy(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::vz(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::wx(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::wy(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::wz(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::px(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::py(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::pz(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::f1(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::f2(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::f3(double t, const AGM::point &pt) -> double {
    double x{pt[0]};
    double y{pt[1]};
    double z{pt[2]};
    return ZEROVALUE;
}

void AGM::NavierStokesFunction::loadPreviousValue(const std::string &filename, std::vector<value> *pu,
                                                  std::vector<value> *pv, std::vector<value> *pw,
                                                  std::vector<value> *pp) {
    int idx{};
    int bc{};
    double x{};
    double y{};
    double z{};
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
        uvel["bdv"] = 2 * u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    } else if (uvel.getCondition() == 'N') {
        uvel["bdv"] = 2 * ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] +
                      2 * uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1] +
                      2 * uz(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[2];
    }
    uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    if (vvel.getCondition() == 'D') {
        vvel["bdv"] = 2 * v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else if (vvel.getCondition() == 'N') {
        vvel["bdv"] = 2 * vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] +
                      2 * vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1] +
                      2 * vz(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[2];
    }
    vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    if (wvel.getCondition() == 'D') {
        wvel["bdv"] = 2 * w(pointHeat::getTime() + pointHeat::getDelta(), wvel);
    } else if (wvel.getCondition() == 'N') {
        wvel["bdv"] = 2 * wx(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[0] +
                      2 * wy(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[1] +
                      2 * wz(pointHeat::getTime() + pointHeat::getDelta(), wvel) * wvel.getNormal()[2];
    }
    wvel["rhs"] = f3(pointHeat::getTime() + pointHeat::getDelta(), wvel);
}

AGM::NavierStokesFunction::~NavierStokesFunction() = default;
