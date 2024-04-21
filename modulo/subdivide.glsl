#define TWO_PI 6.2831853072

#ifndef DIVIDE_ITERS
#define DIVIDE_ITERS 5.
#endif

vec3 subdivide (inout vec2 q, in float seed, in float t) {
  vec2 dMin = vec2(-1.5, -1.5);
  vec2 dMax = vec2( 1.5,  1.5);

  float id = 0.;

  vec2 dim = dMax - dMin;

  // Params
#ifndef MIN_SIZE
  float MIN_SIZE = 0.01;
#endif
#ifndef MIN_ITERS
  float MIN_ITERS = 1.;
#endif

  for (float i = 0.; i < DIVIDE_ITERS; i++) {
    // Noise point to divide box into 4 pieces
    vec2 divHash = vec2(
      noise(vec2(i + id, seed)),
      noise(vec2(i + id + 2.782, seed))
    );

#define GOOD_START 1
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

    // update the box domain
    dMax = mix(  dMax, divide, step(q, divide));
    dMin = mix(divide,   dMin, step(q, divide));

    // id
    vec2 diff2 = mix(-divide, divide, step(q, divide));
    id = length(diff2 + 10.) * 100.;

    // recalculate the dimension
    dim = dMax - dMin;
  }

  vec2 center = (dMin + dMax) * 0.5;
  q -= center;

  return vec3(dim, id);
}

vec3 subdivide (inout vec2 q, in float seed) {
  return subdivide(q, seed, 0.);
}
#pragma glslify: export(subdivide)
