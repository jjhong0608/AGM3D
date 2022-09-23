//
// Created by NIMS-JUNHONG on 2022/02/16.
//

#ifndef AGM3D_FUNCTION_H
#define AGM3D_FUNCTION_H

#include "pointHeat.h"
#include "StdVector"

namespace AGM {
    class ellipticFunction {
    public:
        ellipticFunction();

        virtual ~ellipticFunction();

        auto u(const point &pt) -> double;

        auto phi(const point &pt) -> double;

        auto psi(const point &pt) -> double;

        auto f(const point &pt) -> double;

        auto ux(const point &pt) -> double;

        auto uy(const point &pt) -> double;

        auto uz(const point &pt) -> double;

        void assignBoundaryValue(point &pt);
    };

    class heatFunction {
    public:
        heatFunction();

        virtual ~heatFunction();

        static auto initialTime() -> double;

        static auto terminalTime() -> double;

        static auto deltaTime() -> double;

        auto u(double t, const pointHeat &pt) -> double;

        auto phi(double t, const pointHeat &pt) -> double;

        auto psi(double t, const pointHeat &pt) -> double;

        auto f(double t, const pointHeat &pt) -> double;

        auto ux(double t, const pointHeat &pt) -> double;

        auto uy(double t, const pointHeat &pt) -> double;

        auto uz(double t, const pointHeat &pt) -> double;

        void assignPreviousValue(value &value, pointHeat &pt);

        void assignBoundaryValue(pointHeat &pt);

    };

    class NavierStokesFunction {
    public:
        NavierStokesFunction();

        virtual ~NavierStokesFunction();

        static auto initialTime() -> double;

        static auto terminalTime() -> double;

        static auto deltaTime() -> double;

        static auto writeTime() -> double;

        auto u(double t, const point &pt) -> double;

        auto v(double t, const point &pt) -> double;

        auto w(double t, const point &pt) -> double;

        auto p(double t, const point &pt) -> double;

        auto phiU(double t, const point &pt) -> double;

        auto psiU(double t, const point &pt) -> double;

        auto phiV(double t, const point &pt) -> double;

        auto psiV(double t, const point &pt) -> double;

        auto phiW(double t, const point &pt) -> double;

        auto psiW(double t, const point &pt) -> double;

        auto ux(double t, const point &pt) -> double;

        auto uy(double t, const point &pt) -> double;

        auto uz(double t, const point &pt) -> double;

        auto vx(double t, const point &pt) -> double;

        auto vy(double t, const point &pt) -> double;

        auto vz(double t, const point &pt) -> double;

        auto wx(double t, const point &pt) -> double;

        auto wy(double t, const point &pt) -> double;

        auto wz(double t, const point &pt) -> double;

        auto px(double t, const point &pt) -> double;

        auto py(double t, const point &pt) -> double;

        auto pz(double t, const point &pt) -> double;

        auto f1(double t, const point &pt) -> double;

        auto f2(double t, const point &pt) -> double;

        auto f3(double t, const point &pt) -> double;

        void loadPreviousValue(const std::string &filename, std::vector<value> *pu,
                               std::vector<value> *pv, std::vector<value> *pw,
                               std::vector<value> *pp);

        void assignPreviousValue(value &pu, value &pv, value &pw, value &pp, point &uvel, point &vvel, point &wvel,
                                 point &pres);

        void assignBoundaryValue(point &uvel, point &vvel, point &wvel);

    };
}

#endif //AGM3D_FUNCTION_H
