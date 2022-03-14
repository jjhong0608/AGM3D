//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_AXIALLINE_H
#define AGM3D_AXIALLINE_H

#include "value.h"

namespace AGM {
    class point;

    class plane;

    class axialLine : public std::vector<point *> {
    private:
        char mark{};
        std::array<double, 6> coordinate{};
        std::array<plane *, 2> plane{};

    public:
        axialLine();

        explicit axialLine(char mark);

        virtual ~axialLine();

        [[nodiscard]] char getMark() const;

        void setMark(char i);

        AGM::plane *getPlane(int i) const;

        void setPlane(AGM::plane *pPlane, int i);

        double &operator[](int i);

        double operator-(axialLine &line);
    };
}


#endif //AGM3D_AXIALLINE_H
