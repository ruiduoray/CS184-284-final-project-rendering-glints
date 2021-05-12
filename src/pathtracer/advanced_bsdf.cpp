#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi, Vector3D hit_p) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1;
  return reflectance / abs_cos_theta(*wi);
  // return Vector3D();
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
    // TODO Project 3-2: Part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    double tan_t = tan(acos(h.z));
    double alpha_2 = alpha * alpha;
    return exp(-tan_t * tan_t / alpha_2 - log(PI * alpha_2 * pow(h.z, 4)));
    //return 1.0;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.


  // compute the help variable of eta ^ 2 + k ^ 2
  Vector3D eta_k_srq_sum = eta * eta + k * k;
  Vector3D two_eta_cos_t = 2 * eta * wi.z;
  //render_debugger_node();

  double cos_t_2 = wi.z * wi.z;

  Vector3D R_s = (eta_k_srq_sum - two_eta_cos_t + cos_t_2) / (eta_k_srq_sum + two_eta_cos_t + cos_t_2);
  Vector3D R_p = (eta_k_srq_sum * cos_t_2 - two_eta_cos_t + 1) / (eta_k_srq_sum * cos_t_2 + two_eta_cos_t + 1);

  return (R_s + R_p) / 2;
  return Vector3D();
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi, Vector3D hit_p) {
  // TODO Project 3-2: Part 2
  // Implement microfacet model here.
    if (wo.z > 0 && wi.z > 0) {
        Vector3D h = wo + wi;
        h.normalize();
        return F(wi) * G(wo, wi) * D(h) / (4 * wo.z * wi.z);
    } else
        return Vector3D();

    // return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  // Here is the default setting, cosine cosineHemisphereSampler
/*
  *wi = cosineHemisphereSampler.get_sample(pdf);
  return MicrofacetBSDF::f(wo, *wi);
  */


  // Here is the importance sampling

  double theta_h = atan(sqrt(-alpha * alpha * log(1 - random_uniform())));
  double phi_h = 2 * PI * random_uniform();
  Vector3D h = {sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};
  *wi = 2 * (dot(h, wo)) * h - wo;

  if (wi->z > 0) {
      double p_theta = 2 * sin(theta_h) / (alpha * alpha * pow(cos(theta_h), 3))
              * exp(-tan(theta_h) * tan(theta_h) / alpha / alpha);
      double p_phi = 0.5 / PI;
      *pdf = (float) (p_theta * p_phi / (4 * sin(theta_h) * dot(*wi, h)));

      return f(wo, *wi, Vector3D());
  }
  // return the empty vector
  else {
      *pdf = 0;
      return Vector3D();
  }



}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi, Vector3D hit_p) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
    // TODO Project 3-2: Part 1
    // Implement RefractionBSDF

    float eta = (wo.z >= 0) ? 1 / ior : ior;
    bool flag = refract(wo, wi, ior);
    *pdf = 1;
    if (flag == false){
        return Vector3D();
    }
    else{
        return (transmittance / abs_cos_theta(*wi) / pow(eta, 2));
    }
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi, Vector3D hit_p) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
    // TODO Project 3-2: Part 1
    // Compute Fresnel coefficient and either reflect or refract based on it.

    // compute Fresnel coefficient and use it as the probability of reflection
    // - Fundamentals of Computer Graphics page 305
    if (refract(wo, wi, ior)) {
        float eta = (wo.z >= 0) ? 1 / ior : ior;
        auto r0 = (float) pow((1 - eta) / (1 + eta), 2);
        auto r = r0 + (1 - r0) * (float) pow(1 - abs(wo.z), 5);
        if (coin_flip(r)) { // reflection
            reflect(wo, wi);
            *pdf = r;
            return r * reflectance / abs_cos_theta(*wi);
        } else { // refraction
            *pdf = 1 - r;
            return (1 - r) * transmittance / abs_cos_theta(*wi) / eta / eta;
        }
    } else { // total internal refraction, treat as mirror
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
    }
    //return Vector3D();
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Project 3-2: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = {-wo.x, -wo.y, wo.z};


}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {
    // TODO Project 3-2: Part 1
    // Use Snell's Law to refract wo surface and store result ray in wi.
    // Return false if refraction does not occur due to total internal reflection
    // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
    // ray entering the surface through vacuum.
    float eta = (wo.z >= 0) ? 1 / ior : ior;
    auto discriminant = (float) (1 - eta * eta * (1 - wo.z * wo.z));
    if (discriminant < 0) // total internal reflection
        return false;
    else {
        *wi = {-eta * wo.x, -eta * wo.y, -wo.z / abs(wo.z) * sqrt(discriminant)};
        return true;
    }

    //return true;

}

} // namespace CGL
