#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>
#include <limits>
#include <vector>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bbox;
  int p_count = 0;

  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
    p_count ++;
  }

  BVHNode *node = new BVHNode(bbox);
  if (p_count <= max_leaf_size) {
    // Return leaf node if less or equal to max_leaf_size primitives
    node->start = start;
    node->end = end;
  } else {
    // calculate average of centroids
    Vector3D avg_centroid = Vector3D(0, 0, 0);
    for (auto p = start; p != end; p++) {
      Vector3D centroid = (*p)->get_bbox().centroid();
      avg_centroid += centroid;
    }
    avg_centroid /= p_count;
    // Using average of centroids along the longest axis as the split point
    Vector3D extent = bbox.extent;
    double largest_extent = -1 * std::numeric_limits<double>::infinity();
    int split_axis;
    for(int i = 0; i < 3; ++i) {
      if (extent[i] > largest_extent) {
        largest_extent = extent[i];
        split_axis = i;
      }
    }
    double split_point = avg_centroid[split_axis];
    // Split the primitives using the split point and axis
    std::vector<Primitive *> left;
    std::vector<Primitive *> right;
    for (auto p = start; p != end; p++) {
      Vector3D centroid = (*p)->get_bbox().centroid();
      if (centroid[split_axis] <= split_point) left.push_back(*p);
      else right.push_back(*p);
    }

    // if the split results in one child having no primitives, simply return as leaf node,
    // disregarding the max_leaf_size restriction to prevent infinite recursion
    if (left.size() == 0 || right.size() == 0) {
      node->start = start;
      node->end = end;
      return node;
    }
    
    int offset = 0;
    for (auto p = left.begin(); p != left.end(); p++) {
      *(start + offset) = *p;
      offset += 1;
    }
    for (auto p = right.begin(); p != right.end(); p++) {
      *(start + offset) = *p;
      offset += 1;
    }
    node->l = construct_bvh(start, std::next(start, left.size()), max_leaf_size);
    node->r = construct_bvh(std::next(start, left.size()), end, max_leaf_size);

  }

  return node;

}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  // Handle the base case where node is a leaf node
  if (node->l == NULL && node->r == NULL) {
    for (auto p = node->start; p != node->end; p++) {
      total_isects++;
      if ((*p)->has_intersection(ray)) {
        return true;
      }
    }
    return false;
  }

  // Handle the inner nodes recursively
  BBox bbox = node->bb;
  double tmin, tmax;
  total_isects++;
  if (bbox.intersect(ray, tmin, tmax)) {
    if (tmax < ray.min_t || tmin > ray.max_t) return false;
    bool hitl = has_intersection(ray, node->l);
    bool hitr = has_intersection(ray, node->r);
    return hitl || hitr;
  } else {
    return false;
  }

  // for (auto p : primitives) {
  //   total_isects++;
  //   if (p->has_intersection(ray))
  //     return true;
  // }
  // return false;


}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.

  // Handle the base case where node is a leaf node
  if (node->l == NULL && node->r == NULL) {
    bool hit = false;
    for (auto p = node->start; p != node->end; p++) {
      total_isects++;
      hit = (*p)->intersect(ray, i) || hit;
    }
    return hit;
  }

  // Handle the inner nodes recursively
  BBox bbox = node->bb;
  double tmin, tmax;
  total_isects++;
  if (bbox.intersect(ray, tmin, tmax)) {
    if (tmax < ray.min_t || tmin > ray.max_t) return false;
    bool hitl = intersect(ray, i, node->l);
    bool hitr = intersect(ray, i, node->r);
    return hitl || hitr;
  } else {
    return false;
  }

  // bool hit = false;
  // for (auto p : primitives) {
  //   total_isects++;
  //   hit = p->intersect(ray, i) || hit;
  // }
  // return hit;


}

} // namespace SceneObjects
} // namespace CGL
