#define _PHI_ (1.+sqrt(5.))/2.

#define _IVNORM_ (0.5/_PHI_)
#define _PHI1_ (_PHI_*_IVNORM_)
#define _1PHI_ (_IVNORM_)
#define _PHI2_ (_PHI_*_PHI_*_IVNORM_)

#define _IKVNORM_ 1. / sqrt(pow(_PHI_*(1.+_PHI_),2.) + pow(pow(_PHI_, 2.) - 1., 2.) + pow(1.+_PHI_, 2.))
#define _C1_ (_PHI_*(1.+_PHI_)*_IKVNORM_)
#define _C2_ ((_PHI_*_PHI_-1.)*_IKVNORM_)
#define _1C_ ((1.+_PHI_)*_IKVNORM_)

#pragma glslify: foldNd = require(../foldNd)

vec3 dodecahedronFold (in vec3 p, inout float minD) {
  for(int i=0; i < Iterations; i++) {
    p=abs(p);

    vec3 axis = vec3(_PHI2_, _1PHI_, -_PHI1_);
    foldNd(p, axis);

    axis = vec3(-_PHI1_, _PHI2_, _1PHI_);
    foldNd(p, axis);
    p += 0.10 * cos( 5.0 * p.yzx);
    p += 0.05 * cos(11.0 * p.yzx);

    axis = vec3(_1PHI_, -_PHI1_, _PHI2_);
    foldNd(p, axis);

    axis = vec3(-_C1_, _C2_, _1C_);
    foldNd(p, axis);

    axis = vec3(_1C_, -_C1_, _C2_);
    foldNd(p, axis);

    // Stretch
    p = (vec4(p, 1.) * kifsM).xyz;

    minD = min(length(p), minD);
  }

  return p;
}

#pragma glslify: export(dodecahedronFold)
