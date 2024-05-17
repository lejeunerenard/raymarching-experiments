#define TWO_PI 6.2831853072

#ifndef DIVIDE_ITERS
#define DIVIDE_ITERS 3.
#endif

#pragma glslify: rotMat2 = require(../rotation-matrix2)

vec3 subdivide (inout vec2 q, in float seed, in float t) {
  // Square sizing
  float size = 1.15;
  vec2 dMin = vec2(-size, -size);
  vec2 dMax = vec2( size,  size);

  // // Custom Sizing
  // vec2 dMin = vec2(-0.4, -0.8);
  // vec2 dMax = vec2( 0.4,  0.8);

  float id = 0.;

  vec2 dim = dMax - dMin;

  // Params
#ifndef MIN_SIZE
  float MIN_SIZE = 0.001;
#endif
#ifndef MIN_ITERS
  float MIN_ITERS = 1.;
#endif

  // Include rotation
  // TODO Figure out if this is a meaningful technique
  float accumRot = 0.;

  vec2 currentCenter = (dMin + dMax) * 0.5;
  for (float i = 0.; i < DIVIDE_ITERS; i++) {
    vec2 generationVariability = currentCenter;
    // Noise point to divide box into 4 pieces
    vec2 divHash = 0.5 + 0.3 * vec2(
      noise(vec2(i + 0. * id, seed) + generationVariability),
      noise(vec2(i + 0. * id + 2.782, seed) + generationVariability)
    );

    // Move the dividing line
    divHash += 0.15 * cos(TWO_PI * t + i + dot(generationVariability, vec2(1)));

// #define GOOD_START 1
#ifdef GOOD_START
    if (i == 0.) {
      divHash = vec2(0.49, 0.501);
    }
#endif

    vec2 divide = divHash * dim + dMin;

    // Clamp division Line
    divide = clamp(divide, dMin + MIN_SIZE, dMax - MIN_SIZE);

    // Find the minimum dimension size
    vec2 minAxis = min(abs(dMin - divide), abs(dMax - divide));
    float minSize = vmin(minAxis);

    // check if too small
    bool smallEnough = minSize < MIN_SIZE;
    if (smallEnough && i + 1. > MIN_ITERS) { break; }

    // Rotate space on iterating
    accumRot += 0.065 * TWO_PI;

    // update the box domain
    vec2 overLine = step(q, divide);
    dMax = mix(  dMax, divide, overLine);
    dMin = mix(divide,   dMin, overLine);

    // id
    vec2 diff2 = mix(-divide, divide, overLine);
    id = length(diff2 + 10.) * 100.;

    // recalculate the dimension
    dim = dMax - dMin;
    currentCenter = (dMin + dMax) * 0.5;
  }

  vec2 center = (dMin + dMax) * 0.5;
  q -= center;

  // q *= rotMat2(accumRot);

  return vec3(dim, id);
}

vec3 subdivide (inout vec2 q, in float seed) {
  return subdivide(q, seed, 0.);
}
#pragma glslify: export(subdivide)
