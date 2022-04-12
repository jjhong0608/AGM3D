//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_POINT_H
#define AGM3D_POINT_H

#include <functional>
#include "matrixRow.h"
#include "StdVector"

namespace AGM {
    class point;

    class axialLine;

    using axialElement = std::array<point *, 30>;

    struct interfacePoint {
        std::string pos{};
        point *point{};
    };

    using interfaceElement = std::vector<interfacePoint>;

    class point {
    protected:
        int idx{};
        coordinate xyz{}, normal{};
        double mp{};
        char condition{};
        axialElement element{};
        interfaceElement ielement{};
        value values{};
        std::array<matrixRow, 3> solMatrixRow{}, deriMatrixRow{}, rhsMatrixRow{}, partMatrixRow{};
        std::array<double, 3> rb{}, dv{};
        std::array<axialLine *, 3> aline{};
        static int nPts;
        static std::vector<axialLine> *xline, *yline, *zline;
        static std::array<std::vector<plane>, 2> *xyplane, *yzplane, *xzplane;

    public:
        point();

        explicit point(int idx);

        explicit point(const coordinate &xyz);

        point(int idx, const coordinate &xyz);

        point(const coordinate &xyz, double mp);

        point(int idx, const coordinate &xyz, double mp);

        virtual ~point();

        int getIdx() const;

        void setIdx(int i);

        const coordinate &getXyz() const;

        void setXyz(const coordinate &coordinate);

        const coordinate &getNormal() const;

        void setNormal(const coordinate &coordinate);

        double getMp() const;

        void setMp(double d);

        char getCondition() const;

        void setCondition(char i);

        const std::array<point *, 30> &getElement() const;

        void setElement(const std::array<point *, 30> &array);

        const std::vector<interfacePoint> &getIelement() const;

        void setIelement(const std::vector<interfacePoint> &vector);

        const value &getValue() const;

        void setValue(const value &value);

        const std::array<matrixRow, 3> &getSolMatrixRow() const;

        void setSolMatrixRow(const std::array<matrixRow, 3> &row);

        const std::array<matrixRow, 3> &getDeriMatrixRow() const;

        void setDeriMatrixRow(const std::array<matrixRow, 3> &row);

        const std::array<double, 3> &getRb() const;

        void setRb(const std::array<double, 3> &array);

        const std::array<AGM::axialLine *, 3> &getAxialLine() const;

        axialLine *&getAxialLine(char i);

        void setAxialLine(const std::array<AGM::axialLine *, 3> &array);

        void setAxialLine(AGM::axialLine *line, char i);

        static int getNPts();

        static void setNPts(int i);

        static std::vector<axialLine> *&getAxialLines(char i);

        static void setAxialLines(std::vector<axialLine> *line, char i);

        static std::array<std::vector<plane>, 2> *&getPlane(std::string &str);

        static void setPlane(std::array<std::vector<plane>, 2> *pln, const std::string &str);

        double &operator[](int i);

        const double &operator[](int i) const;

        point *&operator[](STENCIL stencil);

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;

        double operator-(const point &src);

        point &operator=(const point &src);

        void findStencil();

        void findStencilBoundary();

        void calculateRepresentationFormula();

        virtual void calculateRepresentationFormulaCross();

        void calculateRepresentationFormulaDirichlet();

        void calculateRepresentationFormulaNeumann();

        virtual matrixRow calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt);

        virtual matrixRow calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt);

        void approximatePhiAndPsiAtBoundary(int order);

        void updateRightHandSide(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                 const std::function<double(int)> &h);

        void updateRightHandSideCross(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                      const std::function<double(int)> &h);

        void updateRightHandSideDirichlet(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                          const std::function<double(int)> &h);

        void updateRightHandSideNeumann(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                        const std::function<double(int)> &h);

        void updateRightHandSidePart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                     const std::function<double(int)> &h);

        void updateRightHandSideCrossPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                          const std::function<double(int)> &h);

        void updateRightHandSideDirichletPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                              const std::function<double(int)> &h);

        void updateRightHandSideNeumannPart(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                            const std::function<double(int)> &h);

        void makeDerivatives();

        virtual void makeDerivativesCross();

        void makeDerivativesBoundary();

        virtual void calculateDerivatives(const std::vector<point> *points, const std::function<double(int)> &f,
                                          const std::function<double(int)> &g, const std::function<double(int)> &h,
                                          const std::function<double(int)> &fp, const std::function<double(int)> &gp,
                                          const std::function<double(int)> &hp);

        virtual void approximateNaNDerivatives(const std::vector<point> *points);

        void calculateSecondDerivatives(const std::function<double(int)> &f, const std::function<double(int)> &g,
                                        const std::function<double(int)> &h);

    };

}


#endif //AGM3D_POINT_H
