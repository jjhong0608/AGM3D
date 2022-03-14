//
// Created by 조준홍 on 2022/02/03.
//

#ifndef AGM3D_POINTHEAT_H
#define AGM3D_POINTHEAT_H

#include "point.h"

namespace AGM {
    class pointHeat : public point {
    protected:
        static double time, delta;

    public:
        static double getTime();

        static void setTime(double d);

        static double getDelta();

        static void setDelta(double d);

        void findStencil(const axialElement *axialElement1, std::vector<pointHeat> *vector);

        void calculateRepresentationFormulaCross() override;

        matrixRow calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) override;

        matrixRow calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) override;

        void makeDerivativesCross() override;

        void calculateDerivatives(const std::vector<pointHeat> *points, const std::function<double(int)> &f,
                                  const std::function<double(int)> &g, const std::function<double(int)> &h,
                                  const std::function<double(int)> &fp, const std::function<double(int)> &gp,
                                  const std::function<double(int)> &hp);

        void approximateNaNDerivatives(const std::vector<pointHeat> *points);
    };

}


#endif //AGM3D_POINTHEAT_H
