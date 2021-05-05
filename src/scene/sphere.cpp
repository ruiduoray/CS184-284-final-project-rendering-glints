#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.

  double a = dot(r.d, r.d);
  Vector3D diff_o = r.o - o;
  double b = 2 * dot(diff_o, r.d);
  double c = dot(diff_o, diff_o) - r2;
  double root_test = pow(b, 2) - 4 * a * c;

  if (root_test < 0) return false;
  else if (root_test > 0) {
    double t_lo = (-1 * b - sqrt(root_test)) / (2 * a);
    double t_hi = (-1 * b + sqrt(root_test)) / (2 * a);
    t1 = t_lo, t2 = t_hi;
    return true;
  }
  else {
    double t_lo = -1 * b / (2 * a);
    t1 = t_lo, t2 = t_lo;
    return true;
  }

  return false;
}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.

  double t1, t2;
  bool is_i = test(r, t1, t2);

  if (is_i) {
    if (t1 >= r.min_t && t1 <= r.max_t) {
      r.max_t = t1;
      return true;
    } else if (t2 >= r.min_t && t2 <= r.max_t) {
      r.max_t = t2;
      return true;
    } else return false;
  }
  else return false;
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  
  double t1, t2;
  bool is_i = test(r, t1, t2);

  if (is_i) {
    // Update t
    if (t1 >= r.min_t && t1 <= r.max_t) {
      i->t = t1;
      r.max_t = t1;
    } else if (t2 >= r.min_t && t2 <= r.max_t) {
      i->t = t2;
      r.max_t = t2;
    } else return false;
    Vector3D i_pt = r.o + i->t * r.d;
    i->n = (i_pt - o);
    i->n.normalize();
    i->primitive = this;
    i->bsdf = get_bsdf();
    i->hit_p = i_pt;
    
    return true;
  }

  return false;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
