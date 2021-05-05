#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

Vector3D Triangle::moller_trumbore(const Ray &r) const {
  // Calculate t using Moller Trumbore Algorithm
  Vector3D e1 = p2 - p1;
  Vector3D s = r.o - p1;
  Vector3D e2 = p3 - p1;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);
  double t = dot(s2, e2) / dot(s1, e1);
  double b1 = dot(s1, s) / dot(s1, e1);
  double b2 = dot(s2, r.d) / dot(s1, e1);

  return Vector3D(t, b1, b2);
}

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.
  
  // Calculate t, b1, b2 using Moller Trumbore Algorithm
  Vector3D mt = moller_trumbore(r);

  // Intersect only if t falls the range min_t ~ max_t
  if (mt.x >= r.min_t && mt.x <= r.max_t && mt.y >= 0 && mt.y <= 1 && mt.z >= 0 && mt.z <= 1 && mt.y + mt.z <= 1) {
    r.max_t = mt.x;
    return true;
  } else {
    return false;
  }
}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  
  Vector3D mt = moller_trumbore(r);

  // Updating the intersection point parameters
  if (mt.x >= r.min_t && mt.x <= r.max_t && mt.y >= 0 && mt.y <= 1 && mt.z >= 0 && mt.z <= 1 && mt.y + mt.z <= 1) {
    // Updating t
    r.max_t = mt.x;
    isect->t = mt.x;
    // Calculate the intersection coordinates
    Vector3D isect_pt = r.o + mt.x * r.d;
    // Calculate barycentric coordinate of the intersection point
    double area_a = 0.5 * cross(p3 - isect_pt, p2 - isect_pt).norm();
    double area_b = 0.5 * cross(p1 - isect_pt, p3 - isect_pt).norm();
    double area_c = 0.5 * cross(p1 - isect_pt, p2 - isect_pt).norm();
    double alp = area_a / (area_a + area_b + area_c);
    double bet = area_b / (area_a + area_b + area_c);
    double gam = 1 - alp - bet;
    // Use barycentric coordinates to update the normal vector of the intersection point
    isect->n = alp * n1 + bet * n2 + gam * n3;
    isect->n.normalize();
    // Update the primitive and bsdf
    isect->primitive = this;
    isect->bsdf = this->get_bsdf();
    isect->hit_p = isect_pt;

    return true;
  } else {
    return false;
  }

}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
