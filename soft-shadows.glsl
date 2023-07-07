#define RAY_STEPS maxSteps
// TECHNIQUE:
// 0 : Default proportion of closest distance / distance from point
// 1 : More pessimistic closest point based on previous & current sphere
// 2 : Proportion but supports negative distances
#define TECHNIQUE 2

// Sources
// - https://iquilezles.org/articles/rmshadows/
// - https://www.shadertoy.com/view/lsKcDD
// - https://www.shadertoy.com/view/Xds3zN

float softshadow( in vec3 ro, in vec3 rd, in float mint, in float maxt, in float k, in float generalT ) {
  float res = 1.0;
#if TECHNIQUE == 1
  float ph = 1e10; // big so that y = 0 on the first iteration
#endif
  float t = mint;

  for( int i=0; i<RAY_STEPS; i++ ) {
    float h = map(ro + rd*t, generalT).x;

    // Anti-banding by using intersection of previous sphere & current sphere to
    // find the minimum possible distance from object to ray.
#if TECHNIQUE == 1
      // Skip first iteration if you are getting artifact on the first iteration,
      // or unroll the first iteration out of the loop.
      // Note: Not sure why this is necessary as the ph being big should force
      // the same effect.
      // float y = (i == 0) ? 0. : h * h / (2. * ph);

      float y = h * h / (2. * ph); // Standard version

      float d = sqrt(max(h * h - y * y, 0.));
      res = min( res, d / (k * max(0., t - y)) );
      ph = h;
#else
      res = min( res, h / (k * t) );
#endif

#if TECHNIQUE == 2
      // Prevent back stepping
      t += clamp(h , 0.005, 0.5);

      if(res < -1. || t >= maxt) break;
#else
      t += h;

      if(res < 0.0001 || t >= maxt) break;
#endif
  }

#if TECHNIQUE == 2
  res = max(res, -1.);
  return 0.25 * (1. + res) * (1. + res) * (2. - res);
#else
  res = clamp( res, 0.0, 1.0 );

  // Found this in iq's example implementations. Why it's there, I have no
  // clue. I don't notice a difference.
  // res = res * res * (3. - 2. * res);
  return res;
#endif
}

float softshadow( in vec3 ro, in vec3 rd, in float mint, in float maxt ) {
  return softshadow( ro, rd, mint, maxt, 0.125, 0. );
}

#pragma glslify: export(softshadow)
