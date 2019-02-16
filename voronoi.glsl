#pragma glslify: rotMat2 = require(./rotation-matrix2)
#ifndef PI
#define PI 3.1415926536
#endif
#ifndef TWO_PI
#define TWO_PI 6.2831853072
#endif

vec2 hash( vec2 p  ) {
  p = vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3)));
  return fract(sin(p)*18.5453);
}

vec3 hash( vec3 x  )
{
  x = vec3( dot(x,vec3(127.1,311.7, 74.7)),
        dot(x,vec3(269.5,183.3,246.1)),
        dot(x,vec3(113.5,271.9,124.6)) );

  return fract(sin(x)*43758.5453123);
}

vec2 voronoi(in vec2 x, in float time) {
  vec2 p = floor(x);
  vec2 f = fract(x);

  vec3 res = vec3(8.0);
  for (int j=-1; j<=1; j++) {
    for (int i=-1; i<=1; i++) {
      vec2 b = vec2(i, j);
      vec2 offset = hash(p + b + time);
      offset = 0.5 + 0.5 * cos( time + TWO_PI * offset );
      vec2 r = vec2(b) - f + offset;
      float d = dot( r, r );
      // float d = dot( abs(r), vec2(1) );
      // float d = min( abs(r.x), abs(r.y) );
      // vec2 absR = abs(r);
      // float d = max( absR.x, absR.y );

      // Mask output
      // float v = mask((p + b) * 0.090909);
      // if (d < res.x && v < 0.00125) {

      if (d < res.x) {
        res = vec3( d, offset );
      }
    }
  }
  // Add mask to distance calculation
  // float v = mask(x * 0.090909);
  // res.x = min(res.x, -24. * v);
  // res.x = mix(1., res.x, smoothstep(0.00125, 0., v));

  return vec2(sqrt( res.x ), dot(res.yz, vec2(1.0)));
}
vec2 voronoi(in vec2 x) {
  return voronoi(x, 0.0);
}

float metric (in vec3 r) {
  return dot(r, r);
}

vec3 voronoi(in vec3 x, in float time) {
  vec3 p = floor(x);
  vec3 f = fract(x);

  float id = 0.0;
  vec2 res = vec2(100.0);
  for (int k=-1; k<=1; k++)
  for (int j=-1; j<=1; j++)
  for (int i=-1; i<=1; i++)
  {
    vec3 b = vec3(i, j, k);
    vec3 offset = hash(p + b);
    vec3 r = b - f + offset;
    float d = metric(r);

    if (d < res.x) {
      id = dot( p + b, vec3(1.0, 57.0, 113.0) );
      res = vec2(d, res.x);
    } else if ( d < res.y ) {
      res.y = d;
    }
  }

  return vec3(sqrt(res), abs(id));
}

#pragma glslify: export(voronoi)
