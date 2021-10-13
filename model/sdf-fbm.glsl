// #ifndef REPS
// #define REPS 11
// #endif

// source: https://www.iquilezles.org/www/articles/fbmsdf/fbmsdf.htm
float hash( in vec3 x ) {
  return sin(1.5*x.x)*sin(1.5*x.y)*sin(1.5*x.z);
}

float sph ( in ivec3 i, in vec3 f, in ivec3 c ) {
  float rad = 0.65 * hash(vec3(i + c));
  return length(f - vec3(c)) - rad;
}

float sdBase ( in vec3 p ) {
  ivec3 i = ivec3(floor(p));
   vec3 f =       fract(p);

   return min(min(min(sph(i,f,ivec3(0,0,0)),
                      sph(i,f,ivec3(0,0,1))),
                  min(sph(i,f,ivec3(0,1,0)),
                      sph(i,f,ivec3(0,1,1)))),
              min(min(sph(i,f,ivec3(1,0,0)),
                      sph(i,f,ivec3(1,0,1))),
                  min(sph(i,f,ivec3(1,1,0)),
                      sph(i,f,ivec3(1,1,1)))));
}

float sdFBM (in vec3 p, float d, float th){
  float s = 1.;

  for (int i = 0; i < REPS; i++) {
    float n = s * sdBase(p);

#define FBM_SUBTRACT 1
#ifdef FBM_SUBTRACT
    d = smax(d, -n, smoothness * s);
#else
    n = smax(n, d - clipFactor * s, smoothness * s);
    d = smin(n, d                 , smoothness * s);
#endif

    p = mat3( 0.00, 1.60, 1.20,
             -1.60, 0.72,-0.96,
             -1.20,-0.96, 1.28) * p;

    s = 0.5 * s;
    if (s < th) break;
  }

  return d;
}

float sdFBM (in vec3 p, float d) {
  return sdFBM(p, d, 1e-3); // ? not sure about this default number
}

#pragma glslify: export(sdFBM)
