#ifndef IMAGE_TEXTURE_H
#define IMAGE_TEXTURE_H
#include "texture.h"

namespace CGL {

class image_texture : public texture {
    public:
        image_texture() {}
        image_texture(unsigned char *pixels, int A, int B) : data(pixels), nx(A), ny(B) {}

        virtual Vector3D value(float u, float v, const Vector3D& p) const;

        unsigned char *data;
        int nx, ny;
};

Vector3D image_texture::value(float u, float v, const Vector3D& p) const {
    int i = u * nx;
    int j = (1 - v) * ny - 0.001;

    if (i < 0) i = 0;
    if (j < 0) j = 0;

    if (i > nx - 1) i = nx - 1;
    if (j > ny - 1) j = ny - 1;

    float r = int(data[3 * i + 3 * nx * j]) / 255.0;
    float g = int(data[3 * i + 3 * nx * j + 1]) / 255.0;
    float b = int(data[3 * i + 3 * ny * j + 2]) / 255.0;

    return Vector3D(r, g, b);
}

}

#endif