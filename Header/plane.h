//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_PLANE_H
#define AGM3D_PLANE_H

#include "axialLine.h"

namespace AGM {
    class plane : public std::vector<axialLine *> {
    private:
        std::string mark{};
        double coordinate{};

    public:
        plane();

        explicit plane(std::string mark);

        plane(std::string mark, double coordinate);

        virtual ~plane();

        [[nodiscard]] auto getMark() const -> const std::string &;

        void setMark(const std::string &string);

        [[nodiscard]] auto getCoordinate() const -> double;

        void setCoordinate(double d);

        void setPlanePointer();
    };
}


#endif //AGM3D_PLANE_H
