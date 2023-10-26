#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// DONE: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // Calculate position using the sample which is linearly interpolated between the endpoints
    position = (light.endpoint1 - light.endpoint0) * sample + light.endpoint0;

    // Calculate color using the sample which is linearly interpolated between the colors
    color = (light.color1 - light.color0)*sample + light.color0;
}

// DONE: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // Linear interpolation of position based on the parallelogram using the length of the edges and starting at v0
    position = light.v0 + light.edge01 * sample.x + light.edge02 * sample.y;

    // Linear interpolation of color based on the parallelogram using the color 'length' of the edges and starting at color0
    glm::vec3 xOne = light.color0 + (light.color2 - light.color0) * sample.x;
    glm::vec3 xTwo = light.color1 + (light.color3 - light.color1) * sample.x;
    glm::vec3 y = xOne + (xTwo - xOne) * sample.y;
    color = y;
}

// DONE: Standard feature
// Given a sampled position on some light, an+d the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer

        // Intersection of given ray with plane
        glm::vec3 intersection = { ray.origin + ray.t * ray.direction };

        // Normalised light direction from light to intersection with plane
        glm::vec3 lightDirection = glm::normalize(lightPosition - intersection);

        // Light is on other side of the normal, light does not see the point
        if (glm::dot(hitInfo.normal, lightDirection) <= 0) {
            return false;
        }

        float epsilon = 1e-6;
        glm::vec3 lightOrigin = { intersection + epsilon * lightDirection };

        // Time when the ray will hit the light
        int t = glm::length(lightPosition - intersection);

        // Create the ray towards the light starting from the plane, in the direction of the light
        Ray lightRay = Ray(lightOrigin, lightDirection, t);
        HitInfo hit = hitInfo;

        // Check if there is anything between lightRay and the light
        return !state.bvh.intersect(state, lightRay, hit);
    }
}

// DONE: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    float epsilon = 1e-6;

    // Intersection of given ray with plane
    glm::vec3 intersection = { ray.origin + ray.t * ray.direction };
    
    // Normalised light direction from light to intersection with plane
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersection);

    // Time when the ray will hit the light
    int t = glm::length(lightPosition - intersection);

    // Ray of light from the poistion of the light to the intersection point
    Ray lightRay = Ray(intersection + epsilon * lightDirection, lightDirection, t);

    // Variable to store the information of the current hit in the recursion
    HitInfo currentHit = hitInfo;

    // Variable for the color which will be returned
    glm::vec3 currentColor = lightColor;

    // Go along the ray and intersect, update the color based on the info of the object that has been hit
    while (state.bvh.intersect(state, lightRay, currentHit)) {
        // Update current color with given formula
        currentColor = currentColor * currentHit.material.kd * (1 - currentHit.material.transparency);

        // Calculate new origin of intersection with current object
        glm::vec3 origin = { lightRay.origin + lightRay.t * lightRay.direction };

        // Move the origin along the ray to slightly beyond the point where we have intersected an object
        lightRay.origin = { lightRay.origin + (lightRay.t + epsilon) * lightRay.direction };

        // Update t to be correct for the new origin
        lightRay.t = glm::length(lightPosition - origin);
    }
    
    return currentColor;
}

// DONE: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 visible = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);
    // Check if light is visible to ray, if it is then compute the shading using this light
    if (visible != glm::vec3{0,0,0}) {
        glm::vec3 intersection = ray.origin + ray.t * ray.direction;
        glm::vec3 lightDirection = glm::normalize(light.position - intersection);
        glm::vec3 cameraDirection = -ray.direction;
        return computeShading(state, cameraDirection, lightDirection, visible, hitInfo);
    }
    return {0, 0, 0};
}

// DONE: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // DONE: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model

    // Position at the start of the segment
    glm::vec3 position = light.endpoint0;

    // Color at the start of the segment
    glm::vec3 color = light.color0;

    // Variable where all the contributions of each point of the segment are added
    glm::vec3 contributions = { 0, 0, 0 };
    for (int i = 1; i <= numSamples; i++) {
        // Take a sample along the segment
        sampleSegmentLight(state.sampler.next_1d(), light, position, color);
        PointLight point = PointLight(position, color);

        // Compute the contribution of this point light and add it to the contributions
        contributions += computeContributionPointLight(state, point, ray, hitInfo);
    }
    return contributions;
}

// DONE: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // DONE: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
    // 
    // Position at the edge of the parallelogram
    glm::vec3 position = light.v0;

    // Color at the start of the parallelogram
    glm::vec3 color = light.color0;

    // Variable where all the contributions of each point of the parallelogram are added
    glm::vec3 contributions = { 0, 0, 0 };
    for (int i = 1; i <= numSamples; i++) {
        // Take a sample in the parallelogram
        sampleParallelogramLight(state.sampler.next_2d(), light, position, color);
        PointLight point = PointLight(position, color);

        // Compute the contribution of this point light and add it to the contributions
        contributions += computeContributionPointLight(state, point, ray, hitInfo);
    }
    return contributions;
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forwards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}