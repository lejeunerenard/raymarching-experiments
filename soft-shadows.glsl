#define NEW_TECHNIQUE 1
// source: https://www.shadertoy.com/view/lsKcDD

// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax, in float generalT ) {
  const float k = 16.0;

  float res = 1.0;
  float t = mint;
  float ph = 1e10;

  for( int i=0; i<32; i++ ) {
    float h = map(ro + rd*t, generalT).x;
    if( h<0.0001) return 0.;

#ifdef NEW_TECHNIQUE
    // Skip first iteration if you are getting artifact on the first iteration,
    // or unroll the first iteration out of the loop
    // float y = (i == 0) ? 0. : h * h / (2. * ph);
    float y = h * h / (2. * ph); // Standard version

    float d = sqrt(h * h - y * y);
    res = min( res, k * d / max(0., t - y) );
#else
    res = min( res, k * h / t );
#endif
    ph = h;

    t += h;
  }

  res = clamp( res, 0.0, 1.0 );

  // Found this in iq's example implementations. Why it's there, I have no
  // clue. I don't notice a difference.
  res = res * res * (3. - 2. * res);
  return res;
}

float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax ) {
  return softshadow( ro, rd, mint, tmax, 0. );
}

#pragma glslify: export(softshadow)
