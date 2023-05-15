#ifdef saturate
// saturate is being scoped to the module so alias
#define saturate_0(x) saturate(x)
#else
#define saturate(x) clamp(x, 0.0, 1.0)
#endif

// #ifdef maxDistance
// // maxDistance is being scoped to the module so alias
// #define maxDistance_0(x) maxDistance(x)
// #else
// #define maxDistance(x) 8.0
// #endif

// source: https://www.shadertoy.com/view/ml2GWc

// Parameters
#define NUM_STEPS 32
// Dithering to smooth volume rays
#define DITHERING

// Lighting
#define VOLUME_DENSITY 0.5
#define VOLUME_ABSORBTION 0.5

// light function
// return the direction and the length of the light vector
vec4 getLight(vec3 lightPos, vec3 p) {
    vec3 lig = lightPos - p; // light vector
    float l = length(lig); // length of the light vector
    lig = normalize(lig); // normalize it
    return vec4(lig,l);
}

// volume rendering
// depth is the depth buffer (distance to the scene)
vec4 renderVolume(vec3 ro, vec3 rd, in vec2 uv, float depth, in float time) {
  float tmax = min(10., depth); // max distance

  vec4 volumeColorSum = vec4(vec3(0), 1); // color and opacity

  float stepSize = tmax / float(NUM_STEPS); // step size
  float t = 0.; // distance travelled

#ifdef DITHERING
  t += stepSize * noise(uv);
#endif

  for (int i = 0; i < NUM_STEPS; i++) { // raymarching loop
    vec3 p = ro + rd * t; // current point
    float h = VOLUME_DENSITY * volumeNoise(8.*p); // density of the fog

    // lighting
    vec4 lighting = getLight(sunPos, p); // light direction & length of the light vector

    float sha = shadow(p, lighting.xyz, 0., lighting.w, time); // shadow of the fractal (god rays)

    // coloring
    vec3 col = sunColor(p) * sha / (lighting.w*lighting.w); // inverse square law

    volumeColorSum.rgb += h * stepSize * volumeColorSum.a * col; // add the color to the final result
    volumeColorSum.a *= exp(-h * stepSize * VOLUME_ABSORBTION); // beer's law

    if (volumeColorSum.a < .01) break; // optimization
    t += stepSize; // march
  }

  // output
  return volumeColorSum;
}

vec4 godRays ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv, in vec4 color, in float generalT) {
  vec4 volume = renderVolume(rayOrigin, rayDirection, uv, t.x, generalT);

  // // TODO Figure out how to factor in alpha for background
  color.rgb = color.rgb * volume.a + volume.rgb;

  // color = vec4(volume.rgb, 1.);

  // color = mix(color, volume, volume.a);
  // color += volume;

  return color;
}

#pragma glslify: export(godRays)
