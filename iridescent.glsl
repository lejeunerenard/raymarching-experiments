#pragma glslify: hsv = require(glsl-hsv2rgb)

// @cabbibo
// source: https://www.shadertoy.com/view/4llGDB

#define STEPS 10
vec3 iridescant( vec3 ro , vec3 rd, float lumShift ){

  float lum = 1.;


  vec3 col = vec3(0.);
  for( int i = 0; i < STEPS; i++ ){
    vec3 p = ro + rd * .1  * float( i );
    float period = lumPeriod(p, i);

    lum += abs(sin( p.x * period * 20. ) + sin( p.y * period * 20. ));

    col += hsv( vec3(lum / 10. + lumShift, 1. , 1.) ) / lum;
  }

  return col/float(STEPS);
}
vec3 iridescant( vec3 ro , vec3 rd ){
  return iridescant(ro, rd, 0.);
}

#pragma glslify: export(iridescant)
