#ifndef CGL_INTERSECT_H
#define CGL_INTERSECT_H

#include <vector>

#include "CGL/vector3D.h"
#include "CGL/misc.h"

#include "pathtracer/bsdf.h"

namespace CGL {
    
namespace SceneObjects {

class Primitive;

/**
 * A record of an intersection point which includes the time of intersection
 * and other information needed for shading
 */
struct Intersection {

  Intersection() : t (INF_D), primitive(NULL), bsdf(NULL) { }

  double t;    ///< time of intersection

  const Primitive* primitive;  ///< the primitive intersected

  Vector3D n;  ///< normal at point of intersection

  BSDF* bsdf; ///< BSDF of the surface at point of intersection

  Vector3D hit_p;

  float u, v;
  
  public:
    void get_uv(const Vector3D& p, float& u, float& v) {
      float phi = atan2(p.z, p.x);
      float theta = asin(p.y);
      u = 1 - (phi + M_PI) / (2 * M_PI);
      v = (theta + M_PI / 2) / M_PI;
    }

  // More to follow.
};

} // namespace SceneObjects
} // namespace CGL

#endif // CGL_INTERSECT_H
