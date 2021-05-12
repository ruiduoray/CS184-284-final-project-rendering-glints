#ifndef CONSTANT_TEXTURE_H
#define CONSTANT_TEXTURE_H
#include "texture.h"
#include "CGL/vector3D.h"

namespace CGL {

class constant_texture : public texture {
public:
	constant_texture() {}
	constant_texture(Vector3D c) : color(c) {}

	virtual Vector3D value(float u, float v, const Vector3D& p) const {
		return color;
	}

	Vector3D color;
};

};

#endif // CONSTANT_TEXTURE_H