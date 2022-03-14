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

        double u(const point &pt);

        double phi(const point &pt);

        double psi(const point &pt);

        double f(const point &pt);

        double ux(const point &pt);

        double uy(const point &pt);

        double uz(const point &pt);

        void assignBoundaryValue(point &pt);
    };

    class heatFunction {
    public:
        heatFunction();

        virtual ~heatFunction();

        static double initialTime();

        static double terminalTime();

        static double deltaTime();

        double u(double t, const pointHeat &pt);

        double phi(double t, const pointHeat &pt);

        double psi(double t, const pointHeat &pt);

        double f(double t, const pointHeat &pt);

        double ux(double t, const pointHeat &pt);

        double uy(double t, const pointHeat &pt);

        double uz(double t, const pointHeat &pt);

        void assignPreviousValue(value &value, pointHeat &pt);

        void assignBoundaryValue(pointHeat &pt);

    };

    class NavierStokesFunction {
    public:
        NavierStokesFunction();

        virtual ~NavierStokesFunction();

        static double initialTime();

        static double terminalTime();

        static double deltaTime();

        double u(double t, const point &pt);

        double v(double t, const point &pt);

        double w(double t, const point &pt);

        double p(double t, const point &pt);

        double phiU(double t, const point &pt);

        double psiU(double t, const point &pt);

        double phiV(double t, const point &pt);

        double psiV(double t, const point &pt);

        double phiW(double t, const point &pt);

        double psiW(double t, const point &pt);

        double ux(double t, const point &pt);

        double uy(double t, const point &pt);

        double uz(double t, const point &pt);

        double vx(double t, const point &pt);

        double vy(double t, const point &pt);

        double vz(double t, const point &pt);

        double wx(double t, const point &pt);

        double wy(double t, const point &pt);

        double wz(double t, const point &pt);

        double px(double t, const point &pt);

        double py(double t, const point &pt);

        double pz(double t, const point &pt);

        double f1(double t, const point &pt);

        double f2(double t, const point &pt);

        double f3(double t, const point &pt);

        void loadPreviousValue(const std::string &filename, std::vector<value> *pu,
                               std::vector<value> *pv, std::vector<value> *pw,
                               std::vector<value> *pp);

        void assignPreviousValue(value &pu, value &pv, value &pw, value &pp, point &uvel, point &vvel, point &wvel,
                                 point &pres);

        void assignBoundaryValue(point &uvel, point &vvel, point &wvel);

    };
}

#endif //AGM3D_FUNCTION_H
