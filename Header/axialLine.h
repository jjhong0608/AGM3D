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
        static const int coordinateNumber{6};
        std::array<double, coordinateNumber> coordinate{};
        std::array<plane *, 2> plane{};

    public:
        axialLine();

        explicit axialLine(char mark);

        virtual ~axialLine();

        [[nodiscard]] auto getMark() const -> char;

        void setMark(char i);

        auto getPlane(int i) const -> AGM::plane *;

        void setPlane(AGM::plane *pPlane, int i);

        auto operator[](int i) -> double &;

        auto operator-(axialLine &line) -> double;
    };
}


#endif //AGM3D_AXIALLINE_H
