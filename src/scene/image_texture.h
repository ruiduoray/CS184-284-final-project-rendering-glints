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
    // int j = (1 - v) * ny - 0.001;
    int j = (1 - v) * ny;

    if (i < 0) return Vector3D(.5f, .5f, .5f);
    if (j < 0) return Vector3D(.5f, .5f, .5f);

    if (i > nx - 1) return Vector3D(.5f, .5f, .5f);
    if (j > ny - 1) return Vector3D(.5f, .5f, .5f);

    float r = int(data[3 * i + 3 * nx * j]) / 255.0;
    float g = int(data[3 * i + 3 * nx * j + 1]) / 255.0;
    float b = int(data[3 * i + 3 * nx * j + 2]) / 255.0;

    return Vector3D(r, g, b);
}

}

#endif