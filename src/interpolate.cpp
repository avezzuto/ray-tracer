#include "interpolate.h"
#include <glm/geometric.hpp>

// DONE Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    //https://math.stackexchange.com/questions/2455538/finding-barycentric-coordinates-of-a-point-p-in-a-triangle

    glm::vec3 v1v0 = v1 - v0;
    glm::vec3 v2v0 = v2 - v0;
    glm::vec3 pv0 = p - v0;

    glm::vec3 c1 = cross(v1v0, pv0);
    glm::vec3 c2 = cross(v1v0, v2v0);

    float gamma = dot(c1, c2) / dot(c2, c2);
    float beta = dot(cross(pv0, v2v0), c2) / dot(c2, c2);
    float alpha = 1.0f - beta - gamma;

    return glm::vec3(alpha, beta, gamma);
}

// DONE Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    glm::vec3 result = (bc.x * n0 + bc.y * n1 + bc.z * n2);
    return result;
}

// DONE Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    glm::vec2 result = bc.x * t0 + bc.y * t1 + bc.z * t2;
    return result;
}
