#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    float x = (texCoord.x * image.width) - 0.5f;
    float y = (texCoord.y * image.height) - 0.5f;

    int i = static_cast<int>(std::round(x));
    int j = static_cast<int>(std::round(y));

    i = std::clamp(i, 0, image.width - 1);
    j = std::clamp(j, 0, image.height - 1);

    int pixelIndex = i + j * image.width;

    return image.pixels[pixelIndex];

    return image.pixels[0];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    float x = texCoord.x * (image.width - 1); // Convert to 0-based index
    float y = texCoord.y * (image.height - 1); // Convert to 0-based index

    int i0 = static_cast<int>(x);
    int j0 = static_cast<int>(y);
    int i1 = std::min(i0 + 1, image.width - 1);
    int j1 = std::min(j0 + 1, image.height - 1);

    float dx = x - i0;
    float dy = y - j0;

    glm::vec3 texel00 = image.pixels[i0 + j0 * image.width];
    glm::vec3 texel01 = image.pixels[i0 + j1 * image.width];
    glm::vec3 texel10 = image.pixels[i1 + j0 * image.width];
    glm::vec3 texel11 = image.pixels[i1 + j1 * image.width];

    glm::vec3 interpolatedTexel = (1.0f - dy) * ((1.0f - dx) * texel00 + dx * texel10) + dy * ((1.0f - dx) * texel01 + dx * texel11);

    return interpolatedTexel;
}