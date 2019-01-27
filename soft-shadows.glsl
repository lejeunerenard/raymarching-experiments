// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax ) {
  float res = 1.0;
    float t = mint;
    for( int i=0; i<64; i++ ) {
      vec3 h = map(ro + rd*t);
      res = min( res, 10.0 * h.x / t );
      t += h.x;
      if( h.x<0.0001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

#pragma glslify: export(softshadow)
