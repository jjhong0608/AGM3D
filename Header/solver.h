//
// Created by NIMS-JUNHONG on 2022/01/06.
//

#ifndef AGM3D_SOLVER_H
#define AGM3D_SOLVER_H

#include "writeFileMultiple.h"

namespace AGM {
    class solver {
    private:
        std::vector<point> *pts{};

    public:
        explicit solver(std::vector<point> *pts);

        virtual ~solver();

        std::vector<point> *getPts() const;

        void setPts(std::vector<point> *vector);

        void ellipticSolver();

        void heatSolver();

        void NavierStokesSolver();
    };
}


#endif //AGM3D_SOLVER_H
