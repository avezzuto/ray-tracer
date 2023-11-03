#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    const float_t increment = 0.05, startP = 0.2;
 std::vector<float> axisReducedCenters;
 size_t psize = primitives.size();
 // Calculate centers and reduce to the axis coordinate only
 for (int i = 0; i < psize; i++) {
     glm::vec3 center = computePrimitiveCentroid(primitives[i]);
     axisReducedCenters.push_back(center[axis]);
 }
 std::sort(axisReducedCenters.begin(), axisReducedCenters.end());

 // geometric computations
 float_t SA;
 switch (axis) {
 case 0:
     SA = (aabb.upper[1] - aabb.lower[1]) * (aabb.upper[2] - aabb.lower[2]);
     break;
 case 1:
     SA = (aabb.upper[0] - aabb.lower[0]) * (aabb.upper[2] - aabb.lower[2]);
     break;
 default:
     SA = (aabb.upper[0] - aabb.lower[0]) * (aabb.upper[1] - aabb.lower[1]);
     break;
 }
 float_t axisSpan = abs(aabb.upper[axis] - aabb.lower[axis]);
 float_t totalVolume = axisSpan * SA;
 int bucketsIndex = 0;
 size_t buckets[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 for (int i = 0; i < axisReducedCenters.size(); i++) {
     float_t ac = axisReducedCenters[i] - aabb.lower[axis];
     if (ac > axisSpan * ((bucketsIndex + 1) * increment) + startP) {
         if (bucketsIndex == 10) {
             buckets[11] += psize - 1 - i;
             break;
         }
         bucketsIndex++;
     } else {
         buckets[bucketsIndex]++;
     }
 }

 size_t numPrimitivesA = 0;
 size_t numPrimitivesB = psize;
 float_t volumeA = 0;
 float_t volumeB = totalVolume;
 float_t splitMax = INFINITY;
 size_t finalI = 0;
 for (int i = 0; i < 12; i++) {

     // Update the number of primitives in bins A and B
     numPrimitivesA += buckets[i];
     numPrimitivesB -= buckets[i];

     // Update the surface area of AABBs A and B
     volumeA += increment * totalVolume;
     volumeB -= increment * totalVolume;

     // Calculate the cost for bins A and B
     float costA = numPrimitivesA * volumeA / totalVolume;
     float costB = numPrimitivesB * volumeB / totalVolume;

     // Calculate the total cost for the pair
     if (costA + costB < splitMax) {
         splitMax = costA + costB;
         finalI = increment * axisSpan;
     }
 }
 if (finalI == 0 || finalI == psize - 1)
     return (psize + 1) / 2;
 return finalI;
}