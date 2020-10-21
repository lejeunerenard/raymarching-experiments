#define NEW_TECHNIQUE 1
// source: https://www.shadertoy.com/view/lsKcDD

// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax ) {
  const float k = 8.0;

  float res = 1.0;
  float t = mint;
  float ph = 1e10;

  for( int i=0; i<32; i++ ) {
    float h = map(ro + rd*t).x;

#ifdef NEW_TECHNIQUE
    float y = (i == 0) ? 0. : h * h / (2. * ph);
    float d = sqrt(h * h - y * y);
    res = min( res, k * d / max(0., t - y) );
#else
    res = min( res, k * h / t );
#endif
    ph = h;

    t += h;
    if( h<0.0001 || t>tmax ) break;
  }

  return clamp( res, 0.0, 1.0 );
}

#pragma glslify: export(softshadow)
