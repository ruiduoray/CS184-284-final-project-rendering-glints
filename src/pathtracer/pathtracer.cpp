#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"

#include <cmath>

using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D acc = Vector3D(0);
  for (int i = 0; i < num_samples; ++i) {
    // First, draw a random ramdom object-space vector on the hemisphere
    Vector3D sample = hemisphereSampler->get_sample();
    Vector3D sample_world = o2w * sample;
    Vector3D f = isect.bsdf->f(w_out, sample, isect.hit_p); 
    Intersection is;
    Ray r1 = Ray(hit_p + EPS_F * sample_world, sample_world);
    r1.min_t = EPS_F;
    Vector3D li;
    if(!bvh->intersect(r1, &is)) continue;
    li = is.bsdf->get_emission();
    double cosine = cos_theta(sample);
    acc += (f * li * cosine);
  }
  Vector3D L_out = acc * 2 * M_PI / num_samples;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 

  return L_out;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;
  Vector3D wi;
  double distToLight, pdf;
  for (auto p : scene->lights) {
    int num_samples;
    if (p->is_delta_light()) num_samples = 1;
    else num_samples = ns_area_light;
    Vector3D acc;
    for (int i = 0; i < ns_area_light; i++) {
      Vector3D l = p->sample_L(hit_p, &wi, &distToLight, &pdf);
      Vector3D wi_obj = w2o * wi;
      Vector3D f = isect.bsdf->f(w_out, wi_obj, isect.hit_p);
      if (wi_obj.z < 0) continue;
      Ray r1 = Ray(hit_p + EPS_F * wi, wi);
      r1.min_t = EPS_F;
      r1.max_t = distToLight - EPS_F;
      Intersection is;
      Vector3D li;
      if(bvh->intersect(r1, &is)) li = is.bsdf->get_emission();
      else li = l;
      double cosine = cos_theta(wi_obj);
      acc += (f * li * cosine) / pdf;
    }
    L_out += (acc / ns_area_light);
  }

  return L_out;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light
  
  return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
  
  if (direct_hemisphere_sample) return estimate_direct_lighting_hemisphere(r, isect);
  else return estimate_direct_lighting_importance(r, isect);


}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out = one_bounce_radiance(r, isect);

  Vector3D wi; 
  double pdf;
  Vector3D f = isect.bsdf->sample_f(w_out, &wi, &pdf);
  if (pdf == 0) return L_out;
  Vector3D wi_world = o2w * wi;
  Intersection is;
  Ray r1 = Ray(hit_p + EPS_F * wi_world, wi_world);
  r1.d.normalize();
  r1.min_t = EPS_F;
  if(!bvh->intersect(r1, &is)) return L_out;
  r1.depth = r.depth - 1;

  if (r.depth == max_ray_depth) L_out += (at_least_one_bounce_radiance(r1, is) * f * cos_theta(wi) / pdf);
  
  if (r.depth > 0 && coin_flip(0.65)) L_out += (at_least_one_bounce_radiance(r1, is) * f * cos_theta(wi) / pdf / 0.65);


  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  // if (!bvh->intersect(r, &isect))
  //   return envLight ? envLight->sample_dir(r) : L_out;
  // return L_out;

  // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  if (!bvh->intersect(r, &isect)) return L_out;

  if (r.depth == 0) return zero_bounce_radiance(r, isect);
  else if (r.depth == 1) return zero_bounce_radiance(r, isect) + one_bounce_radiance(r, isect);
  else return zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect);

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"
  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D acc = Vector3D(0, 0, 0);

  // Keep track of s1 and s2 for adaptive sampling
  double s1 = 0, s2 = 0;

  // Generate ns_aa number of random rays from pixel (x, y)
  int i = 0;
  for (; i < num_samples; ++i) {
    Vector2D sample = gridSampler->get_sample();
    Vector2D coord = Vector2D(x + sample.x, y + sample.y);
    Ray r = camera->generate_ray(coord.x / sampleBuffer.w, coord.y / sampleBuffer.h);
    r.depth = max_ray_depth;
    Vector3D radiance = est_radiance_global_illumination(r);
    acc += radiance;
    // Perform adaptive sampling
    s1 += radiance.illum();
    s2 += pow(radiance.illum(), 2);
    int n = i + 1;
    // Check convergence every samplesPerBatch samples
    if (n % (int)samplesPerBatch == 0) {
      double mean = s1 / n;
      double var = (1 / ((double)n - 1)) * (s2 - pow(s1, 2) / (double)n);
      double I = 1.96 * sqrt(var) / sqrt(n);
      // Check convergence
      if (I <= maxTolerance * mean) {
        break;
      }
    }
  }
  acc /= (i + 1);

  sampleBuffer.update_pixel(acc, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = i + 1;  
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
