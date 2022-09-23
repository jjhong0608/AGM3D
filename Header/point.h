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

    static const int elementNumber = 30;
    using axialElement = std::array<point *, elementNumber>;

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

        auto getIdx() const -> int;

        void setIdx(int i);

        auto getXyz() const -> const coordinate &;

        void setXyz(const coordinate &coordinate);

        auto getNormal() const -> const coordinate &;

        void setNormal(const coordinate &coordinate);

        auto getMp() const -> double;

        void setMp(double d);

        auto getCondition() const -> char;

        void setCondition(char i);

        auto getElement() const -> const std::array<point *, elementNumber> &;

        void setElement(const std::array<point *, elementNumber> &array);

        auto getIelement() const -> const std::vector<interfacePoint> &;

        void setIelement(const std::vector<interfacePoint> &vector);

        auto getValue() const -> const value &;

        void setValue(const value &value);

        auto getSolMatrixRow() const -> const std::array<matrixRow, 3> &;

        void setSolMatrixRow(const std::array<matrixRow, 3> &row);

        auto getDeriMatrixRow() const -> const std::array<matrixRow, 3> &;

        void setDeriMatrixRow(const std::array<matrixRow, 3> &row);

        auto getRb() const -> const std::array<double, 3> &;

        void setRb(const std::array<double, 3> &array);

        auto getAxialLine() const -> const std::array<AGM::axialLine *, 3> &;

        auto getAxialLine(char i) -> axialLine *&;

        void setAxialLine(const std::array<AGM::axialLine *, 3> &array);

        void setAxialLine(AGM::axialLine *line, char i);

        static auto getNPts() -> int;

        static void setNPts(int i);

        static auto getAxialLines(char i) -> std::vector<axialLine> *&;

        static void setAxialLines(std::vector<axialLine> *line, char i);

        static auto getPlane(std::string &str) -> std::array<std::vector<plane>, 2> *&;

        static void setPlane(std::array<std::vector<plane>, 2> *pln, const std::string &str);

        auto operator[](int i) -> double &;

        auto operator[](int i) const -> const double &;

        auto operator[](STENCIL stencil) -> point *&;

        auto operator[](const std::string &string) -> double &;

        auto operator[](const std::string &string) const -> const double &;

        auto operator-(const point &src) -> double;

        auto operator=(const point &src) -> point &;

        void findStencil();

        void findStencilBoundary();

        void calculateRepresentationFormula();

        virtual void calculateRepresentationFormulaCross();

        void calculateRepresentationFormulaDirichlet();

        void calculateRepresentationFormulaNeumann();

        virtual auto calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> matrixRow;

        virtual auto calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> matrixRow;

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
