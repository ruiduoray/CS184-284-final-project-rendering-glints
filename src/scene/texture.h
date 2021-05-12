#ifndef TEXTURE_H
#define TEXTURE_H

#include "CGL/vector3D.h"

namespace CGL {

class texture {
    public:
        virtual Vector3D value(float u, float v, const Vector3D& p) const = 0;
};

};

#endif //TEXTURE_H