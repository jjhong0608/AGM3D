//
// Created by NIMS-JUNHONG on 2022/01/07.
//

#ifndef AGM3D_WRITEFILE_H
#define AGM3D_WRITEFILE_H

#include "readFile.h"

namespace AGM {
    template<typename T>
    class writeFile {
    protected:
        T *pt{};
        std::vector<T> *pts{};

    public:
        writeFile();

        explicit writeFile(std::vector<T> *pts);

        T *getPt() const;

        void setPt(T *t);

        std::vector<T> *getPts() const;

        void setPts(std::vector<T> *vector);

        auto calculateErrorAtPoint(const std::string &string);

        double calculateError(const std::string &string);

        void writeResult(const std::string &string);
    };

}


#endif //AGM3D_WRITEFILE_H
