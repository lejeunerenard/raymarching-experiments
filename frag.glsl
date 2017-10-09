#extension GL_OES_standard_derivatives : enable
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)
#define saturate(x) clamp(x, 0.0, 1.0)

// #define debugMapCalls
// #define debugMapMaxed
// #define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform bool BLOOM;
uniform vec3 cOffset;
uniform vec3 cameraRo;
uniform vec4 offsetC;
uniform mat4 cameraMatrix;
uniform mat4 orientation;
uniform mat4 projectionMatrix;
uniform sampler2D tMatCap;
uniform sampler2D audioTexture;

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 384
#define maxDistance 50.0
#define fogMaxDistance maxDistance

#define slowTime time * .2

vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: cnoise4 = require(glsl-noise/classic/4d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
//#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
//#pragma glslify: pnoise3 = require(glsl-noise/periodic/3d)
#pragma glslify: vmax = require(./hg_sdf/vmax)

#define combine(v1, v2, t, p) mix(v1, v2, t/p)

#pragma glslify: rotationMatrix = require(./rotation-matrix3)

float ncnoise2(in vec2 x) {
  return smoothstep(-1.00, 1.00, cnoise2(x));
}
// float ncnoise3(in vec3 x) {
//   return smoothstep(-1.00, 1.00, cnoise3(x));
// }

// 3D noise function (IQ)
float noise(vec3 p) {
  vec3 ip=floor(p);
    p-=ip;
    vec3 s=vec3(7,157,113);
    vec4 h=vec4(0.,s.yz,s.y+s.z)+dot(ip,s);
    p=p*p*(3.-2.*p);
    h=mix(fract(sin(h)*43758.5),fract(sin(h+s.x)*43758.5),p.x);
    h.xy=mix(h.xz,h.yw,p.y);
    return mix(h.x,h.y,p.z);
}
// source: https://www.shadertoy.com/view/lsl3RH
float noise( in vec2 x ) {
  return sin(1.5*x.x)*sin(1.5*x.y);
}
float sinoise3( in vec3 x ) {
  return sin(1.5 * x.x) * sin(1.51 * x.y) * sin(1.52 * x.z * x.x);
}

float iqFBM (vec2 p) {
  float f = 0.0;

  f += 0.500000*cnoise2( p ); p = p*2.02;
  f += 0.250000*cnoise2( p ); p = p*2.03;
  f += 0.125000*cnoise2( p ); p = p*2.01;
  f += 0.062500*cnoise2( p ); p = p*2.025;

  return f * 1.066667;
}

float vfbm4 (vec2 p) {
  float f = 0.0;
  float a = PI * 0.173;
  mat2 m = mat2(
    cos(a), sin(a),
    -sin(a), cos(a));

  f += 0.500000 * noise( p ); p *= m * 2.02;
  f += 0.250000 * noise( p ); p *= m * 2.03;
  f += 0.125000 * noise( p ); p *= m * 2.01;
  f += 0.062500 * noise( p ); p *= m * 2.025;

  return f * 0.9375;
}

float vfbm4 (vec3 p) {
  float f = 0.0;
  const float a = 0.523;
  mat3 m = rotationMatrix(vec3(1, 0, 0), a);

  f += 0.500000 * noise( p ); p *= m * 2.02;
  f += 0.250000 * noise( p ); p *= m * 2.03;
  f += 0.125000 * noise( p ); p *= m * 2.01;
  f += 0.062500 * noise( p ); p *= m * 2.025;

  return f * 0.9375;
}

float vfbm6 (vec2 p) {
  float f = 0.0;
  float a = 1.123;
  mat2 m = mat2(
    cos(a), sin(a),
    -sin(a), cos(a));

  f += 0.500000 * (0.5 + 0.5 * noise( p )); p *= m * 2.02;
  f += 0.250000 * (0.5 + 0.5 * noise( p )); p *= m * 2.03;
  f += 0.125000 * (0.5 + 0.5 * noise( p )); p *= m * 2.01;
  f += 0.062500 * (0.5 + 0.5 * noise( p )); p *= m * 2.025;
  f += 0.031250 * (0.5 + 0.5 * noise( p )); p *= m * 2.011;
  f += 0.015625 * (0.5 + 0.5 * noise( p )); p *= m * 2.0232;

  return f * 0.9375;
}

float vfbm6 (vec3 p) {
  float f = 0.0;
  const float a = 1.123;
  mat3 m = rotationMatrix(vec3(1, 0, 0), a);

  f += 0.500000 * (0.5 + 0.5 * noise( p )); p *= m * 2.02;
  f += 0.250000 * (0.5 + 0.5 * noise( p )); p *= m * 2.03;
  f += 0.125000 * (0.5 + 0.5 * noise( p )); p *= m * 2.01;
  f += 0.062500 * (0.5 + 0.5 * noise( p )); p *= m * 2.025;
  f += 0.031250 * (0.5 + 0.5 * noise( p )); p *= m * 2.011;
  f += 0.015625 * (0.5 + 0.5 * noise( p )); p *= m * 2.0232;

  return f * 0.9375;
}

float iqFBM (vec3 p) {
  float f = 0.0;

  f += 0.500000*noise( p ); p = p*2.02;
  f += 0.250000*noise( p ); p = p*2.03;
  f += 0.125000*noise( p ); p = p*2.01;
  // f += 0.062500*noise( p ); p = p*2.025;

  return f * 1.066667;
}

float fbmWarp (vec2 p, out vec2 q, out vec2 s, out vec2 r) {
  const float scale = 4.0;

  q = vec2(
        iqFBM(p + vec2(0.0, 0.0)),
        iqFBM(p + vec2(3.2, 34.5)));

  s = vec2(
        iqFBM(p + scale * q + vec2(23.9, 234.0)),
        iqFBM(p + scale * q + vec2(3.2, 852.0)));

  r = vec2(
        iqFBM(p + scale * s + vec2(23.9, 234.0)),
        iqFBM(p + scale * s + vec2(3.2, 852.0)));

  return iqFBM(p + scale * r);
}
float fbmWarp (vec2 p, out vec2 q) {
  vec2 s = vec2(0);
  vec2 r = vec2(0);
  return fbmWarp(p, q, s, r);
}
float fbmWarp (vec2 p) {
  vec2 q = vec2(0);
  vec2 s = vec2(0);
  vec2 r = vec2(0);
  return fbmWarp(p, q, s, r);
}

float fbmWarp (vec3 p, out vec3 q, out vec3 s, vec3 r) {
  const float scale = 4.0;

  q = vec3(
        iqFBM(p + vec3(0.0, 0.0, 0.0)),
        iqFBM(p + vec3(3.2, 34.5, .234)),
        iqFBM(p + vec3(7.0, 2.9, -2.42)));

  s = vec3(
        iqFBM(p + scale * q + vec3(23.9, 234.0, -193.0)),
        iqFBM(p + scale * q + vec3(3.2, 852.0, 23.42)),
        iqFBM(p + scale * q + vec3(7.0, -232.0, -2.42)));

  return iqFBM(p + scale * s);
}
float fbmWarp (vec3 p, out vec3 q) {
  vec3 s = vec3(0);
  vec3 r = vec3(0);
  return fbmWarp(p, q, r, s);
}
float vfbmWarp (vec3 p, out vec3 q, out vec3 s, vec3 r) {
  const float scale = 4.0;

  q = vec3(
        vfbm4(p + vec3(0.0, 0.0, 0.0)),
        vfbm4(p + vec3(3.2, 34.5, .234)),
        vfbm4(p + vec3(7.0, 2.9, -2.42)));

  // s = vec3(
  //       vfbm4(p + scale * q + vec3(23.9, 234.0, -193.0)),
  //       vfbm4(p + scale * q + vec3(3.2, 852.0, 23.42)),
  //       vfbm4(p + scale * q + vec3(7.0, -232.0, -2.42)));

  return vfbm6(p + scale * q);
}
float vfbmWarp (vec3 p, out vec3 q) {
  vec3 s = vec3(0);
  vec3 r = vec3(0);
  return vfbmWarp(p, q, r, s);
}
float vfbmWarp (vec3 p) {
  vec3 q = vec3(0);
  vec3 s = vec3(0);
  vec3 r = vec3(0);
  return vfbmWarp(p, q, r, s);
}

#pragma glslify: import(./background)

// Orbit Trap
float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
}

float fOpIntersectionRound(float a, float b, float r) {
  vec2 u = max(vec2(r + a,r + b), vec2(0));
  return min(-r, max (a, b)) + length(u);
}

float fOpDifferenceRound (float a, float b, float r) {
  return fOpIntersectionRound(a, -b, r);
}

// IQ's capsule
float sdCapsule( vec3 p, vec3 a, vec3 b, float r ) {
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}
float sdCone( vec3 p, vec2 c ) {
    // c must be normalized
    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}

// Cone with correct distances to tip and base circle. Y is up, 0 is in the middle of the base.
float fCone(vec3 p, float radius, float height) {
  vec2 q = vec2(length(p.xz), p.y);
  vec2 tip = q - vec2(0.0, height);
  vec2 mantleDir = normalize(vec2(height, radius));
  float mantle = dot(tip, mantleDir);
  float d = max(mantle, -q.y);
  float projected = dot(tip, vec2(mantleDir.y, -mantleDir.x));

  // distance to tip
  if ((q.y > height) && (projected < 0.0)) {
    d = max(d, length(tip));
  }

  // distance to base ring
  if ((q.x > radius) && (projected > length(vec2(height, radius)))) {
    d = max(d, length(q - vec2(radius, 0.0)));
  }
  return d;
}

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}
float udRoundBox( vec3 p, vec3 b, float r ) {
  return length(max(abs(p)-b,0.0))-r;
}

float length16 (in vec2 p) {
  return pow(pow(p.x, 16.0) + pow(p.y, 16.0), 0.125);
}
float length8 (in vec2 p) {
  return pow(p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x
+ p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y, 0.125);
}
float length8 (in vec3 p) {
  return pow(p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x
+ p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y
+ p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z, 0.125);
}

float sdTorus88( vec3 p, vec2 t ) {
  vec2 q = vec2(length16(p.xz)-t.x,p.y);
  return length16(q)-t.y;
}

float sdHexPrism( vec3 p, vec2 h ) {
    vec3 q = abs(p);
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}
float sdTorus( vec3 p, vec2 t ) {
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdPlane( vec3 p, vec4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

// Endless "corner"
float fCorner (vec2 p) {
  return length(max(p, vec2(0))) + vmax(min(p, vec2(0)));
}

#define Iterations 12
// #pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
// #pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM, Iterations=Iterations)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
// #pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

mat4 octaM = mat4(
scale, 0., 0., 0.,
0., scale, 0., 0.,
0., 0., scale, 0.,
1., 1., 1., 1.) * mat4(
1., 0., 0.1, .0,
0., 1., 0., 0.,
0., 0.2, 1., 0.,
0., 0., 0., 1.);
#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=2, kifsM=kifsM, trapCalc=trapCalc)
// 
// #pragma glslify: fold = require(./folds)
#pragma glslify: foldNd = require(./foldNd)
#pragma glslify: twist = require(./twist)

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
  vec2 u = max(vec2(r - a,r - b), vec2(0));
  return max(r, min (a, b)) - length(u);
}

vec2 dMin (vec2 d1, vec2 d2) {
  return (d1.x < d2.x) ? d1 : d2;
}
vec3 dMin (vec3 d1, vec3 d2) {
  return (d1.x < d2.x) ? d1 : d2;
}

// Smooth versions
vec2 dSMin (vec2 d1, vec2 d2, in float r) {
  float d = fOpUnionRound(d1.x, d2.x, r);
  return vec2(d, (d1.x < d2.x) ? d1.y : d2.y);
}
vec3 dSMin (vec3 d1, vec3 d2, in float r) {
  float d = fOpUnionRound(d1.x, d2.x, r);
  return vec3(d, (d1.x < d2.x) ? d1.yz : d2.yz);
}

vec3 dMax (vec3 d1, vec3 d2) {
  return (d1.x > d2.x) ? d1 : d2;
}

float gRAngle = TWO_PI * 0.025 * time;
float gRc = cos(gRAngle);
float gRs = sin(gRAngle);
mat3 globalRot = mat3(
  gRc, 0.0, -gRs,
  0.0, 1.0,  0.0,
  gRs, 0.0,  gRc);

#pragma glslify: rotMat2 = require(./rotation-matrix2)

// IQ
float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

// p as usual, e exponent (p in the paper), r radius or something like that
// #pragma glslify: octahedral = require(./model/octahedral)
// #pragma glslify: dodecahedral = require(./model/dodecahedral)
// #pragma glslify: icosahedral = require(./model/icosahedral)

bool isMaterial( float m, float goal ) {
  return m < goal + 1. && m > goal - .1;
}
float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

// #pragma glslify: pModInterval1 = require(./hg_sdf/p-mod-interval1)
// #pragma glslify: pMod1 = require(./hg_sdf/p-mod1.glsl)
// #pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
// #pragma glslify: pMod3 = require(./hg_sdf/p-mod3.glsl)
// #pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
// #pragma glslify: quad = require(glsl-easings/quintic-in-out)
// #pragma glslify: cub = require(glsl-easings/cubic-in-out)
// #pragma glslify: circ = require(glsl-easings/circular-in-out)
// #pragma glslify: quart = require(glsl-easings/quadratic-in-out)
// #pragma glslify: ease = require(glsl-easings/elastic-in-out)
// #pragma glslify: voronoi = require(./voronoi)
#pragma glslify: band = require(./band-filter)
// #pragma glslify: tetrahedron = require(./model/tetrahedron)
float nsin (in float t) {
  return 0.5 + 0.5 * sin(TWO_PI * t);
}
vec3 nsin (in vec3 t) {
  return 0.5 + 0.5 * sin(TWO_PI * t);
}


// Logistic function
float sigmoid ( in float x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

vec3 theColor (vec2 uv) {
  return mix(#10FFFF, #03E8A8, uv.y);
}

#define jTrap 9
float julia (in vec3 p, in vec4 c) {
    vec4 z = vec4(p, 0.);

    float dz2 = 1.;
    float mz2 = dot(z,z);
    for (int i = 0; i < jTrap; i++) {
        // ∆z[n+1]^2 = 4z[n]^2∆z[n]^2
        dz2 *= 4. * mz2;

        // [n+1] = z[n]^2 + c
        z = vec4(z.x*z.x - dot(z.yzw, z.yzw),
               2.*z.x*z.yzw) + c;

        mz2 = dot(z,z);
        if (mz2 > 4.) break;
    }

  return 0.25 * sqrt(mz2 / dz2) * log(mz2);
}

// Return value is (distance, material, orbit trap)
vec3 mPos = vec3(0);
vec3 map (in vec3 p) {
  vec3 d = vec3(maxDistance, 0, 0);
  float minD = maxDistance;

  const float period = 4.0;
  const float transitionTime = 0.25;
  float modTime = mod(slowTime, period);
  float mixT = saturate((modTime - transitionTime) / (period - transitionTime));

  p *= globalRot;
  vec3 q = p;

  vec3 qRot = q;

  vec3 cosT = PI * 0.5 * vec3(-slowTime, 0.0, slowTime);
  qRot += 0.2500 * cos( 3.0 * qRot.yzx + cosT);
  qRot += 0.1250 * cos( 7.0 * qRot.yzx + cosT);
  qRot += 0.0625 * cos(11.0 * qRot.yzx + cosT);
  float a = atan(qRot.x, qRot.y);

  qRot *= rotationMatrix(normalize(vec3(1, 1, 1)), PI + a);

  // q = qRot;
  // q = mix(q, qRot, nsin(q.y + 1.0));
  q = mix(q, qRot, q.y + 1.0);

  float l = length(q.xy);
  a = atan(q.x, q.y);
  a -= PI * slowTime;

  mPos = q.xyz;

  float rMod = sin(TWO_PI * 27.0 * l + 13.0 * a) + noise(vec2(9.0 * l, sin(4.0 * a)));
  vec3 c = vec3(sdTorus(q.xzy, vec2(0.5, 0.2 + 0.0125 * rMod)), 0.0, rMod);
  c.x *= 0.075;
  d = dMin(d, c);

  // Crop
  float s = length(p) - 0.8;
  d.x = max(d.x, s);

  return d;
}

vec4 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  for (int i = 0; i < maxSteps; i++) {
    vec3 d = map(rayOrigin + rayDirection * t);
    if (d.x < epsilon) return vec4(t + d.x, d.y, float(i), d.z);
    t += d.x;
    maxI = float(i);
    trap = d.z;
    if (t > maxDistance) break;
  }
  return vec4(-1., 0., maxI, trap);
}

#pragma glslify: getNormal = require(./get-normal, map=map)
vec3 getNormal2 (in vec3 p, in float eps) {
  vec2 e = vec2(1.,0.) * 0.01 * eps;
  return normalize(vec3(
    map(p + e.xyy).x - map(p - e.xyy).x,
    map(p + e.yxy).x - map(p - e.yxy).x,
    map(p + e.yyx).x - map(p - e.yyx).x));
}

// Material Functions
float diffuse (in vec3 nor, in vec3 lightPos) {
  return saturate(dot(lightPos, nor));
}

#pragma glslify: softshadow = require(./soft-shadows, map=map)
#pragma glslify: calcAO = require(./ao, map=map)

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#043210, color, 1. - l * .0625);
  // Dark
  color = mix(#ef78FF, color, clamp(exp(l) * .325, 0., 1.));
}

#pragma glslify: checker = require(glsl-checker)
#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: rgb2hsv = require(./rgb2hsv.glsl)
#pragma glslify: debugColor = require(./debug-color-clip)

const float n1 = 1.0;
const float n2 = 1.4;
const float amount = 0.0625;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float spread = 1.0; // saturate(dot(rd, gRd));
  // float n = smoothstep(0.75, 1.0, sin(25.0 * rd.x + 0.01 * noise(433.0 * rd)));
  float n = cnoise3(1.25 * rd + 0.0125 * noise(433.0 * rd));
  float v = n;

  color = vec3(v * spread);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd, in float ior) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=amount)

float dispersionMarch (in vec3 rayDirection) {
  vec3 rayOrigin = gPos + -gNor * 0.001;
  rayDirection = normalize(rayDirection);

  float t = 0.0;

  for (int i = 0; i < 50; i++) {
    float d = map(rayOrigin + rayDirection * t).x;
    if (d >= 0.0) break;
    d = min(d, -0.0001);

    t += abs(d);
  }
  return t;
}

vec3 secondRefraction (in vec3 rd, in float ior) {
  float d = 0.0;

  #if 1
  d = dispersionMarch(rd);
  #else
  const int samples = 12;
  for (int i = 0; i < samples; i++) {
    vec3 lightDir = rd;
    lightDir += 0.1 * vec3(
    noise(gPos),
    noise(1.3 * gPos + 340.0),
    noise(-1.9 * gPos + 640.0));

    d += dispersionMarch(lightDir);
  }
  d /= float(samples);
  #endif

  vec3 reflectionPoint = gPos + rd * d;
  vec3 reflectionPointNor = getNormal2(reflectionPoint, 0.1);
  dNor = reflectionPointNor;

  float sss = saturate(1.0 - pow(length(reflectionPoint - gPos) * 0.48, 10.0));
  vec3 disp = scene(refract(rd, reflectionPointNor, ior), ior);
  disp *= min(1.0, 1.0 / (d * d * 5.0));

  return disp;
}

#pragma glslify: dispersionStep1 = require(./glsl-dispersion, scene=secondRefraction, amount=amount)

#pragma glslify: gradient = require(./gradient)

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  vec3 color = vec3(0);

  float dNor1 = dot(nor, vec3(1)) + 0.1 * noise(4293.0 * mPos);
  color += #FF0000 * nsin(vec3(
    1.1 * dNor1 + 0.03423,
    1.0 * dNor1 + 0.20343,
    0.9 * dNor1 + 0.23034));
  color += #00FFFF * nsin(vec3(
    1.15 * dNor1 + 0.93,
    1.00 * dNor1 + 0.732,
    1.05 * dNor1 + 0.2348));

  trap *= 0.2;

  color += #FFFF00 * nsin(vec3(
    1.05 * trap + 3.23490,
    1.00 * trap + 2.82734,
    0.95 * trap + 7.234034));
  color += #0000FF * nsin(vec3(
    1.05 * trap + 2.9376034,
    1.00 * trap + 1.9376034,
    0.95 * trap + 6.9376034));

  return color;
}

#pragma glslify: reflection = require(./reflection, getNormal=getNormal2, diffuseColor=baseColor, map=map, maxDistance=maxDistance, epsilon=epsilon, maxSteps=maxSteps)

const vec3 glowColor = pow(#ED4F2C, vec3(2.2));

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)
#pragma glslify: matCap = require(./matCap, texture=tMatCap)

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    gPos = pos;
    if (t.x>0.) {
      vec3 color = vec3(0.0);

      vec3 nor = getNormal2(pos, 0.28 * t.x);
      nor += 0.05 * vec3(
          cnoise3(690.0 * mPos),
          cnoise3(870.0 * mPos + 234.634),
          cnoise3(910.0 * mPos + 23.4634));
      nor = normalize(nor);
      gNor = nor;

      vec3 ref = reflect(rayDirection, nor);
      ref = normalize(ref);

      gRd = rayDirection;

      // Basic Diffusion
      vec3 diffuseColor = baseColor(pos, nor, rayDirection, t.y, t.w);

      // Declare lights
      struct light {
        vec3 position;
        vec3 color;
        float intensity;
      };
      const int NUM_OF_LIGHTS = 3;
      const float repNUM_OF_LIGHTS = 0.333333;
      light lights[NUM_OF_LIGHTS];
      lights[0] = light(normalize(vec3(0.4, -0.5, 1.)), #FFFFFF, 1.0);
      lights[1] = light(normalize(vec3(-1., .75, 0.5)), #FFFFFF, 1.0);
      lights[2] = light(normalize(vec3(-0.75, 1.0, 1.0)), #FFFFFF, 1.0);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      const float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 1.0;
      float specCo = 1.0;
      float disperCo = 0.5;

      float specAll = 0.0;
      for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
        vec3 lightPos = lights[i].position;
        float dif = diffuse(nor, lightPos);
        float spec = pow(clamp( dot(ref, (lightPos)), 0., 1. ), 16.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        dif *= max(0.5, saturate(softshadow(pos, lightPos, 0.15, 1.75)));
        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        lin += specCo * spec * (1. - fre);
        specAll += specCo * spec * (1. - fre);

        // Ambient
        lin += 0.2 * amb * diffuseColor;

        color +=
          saturate((occ * dif * lights[i].intensity) * lights[i].color * diffuseColor)
          + saturate(lights[i].intensity * mix(lights[i].color, vec3(1), 0.1) * lin * mix(diffuseColor, #ffffff, 0.4));
      }

      color *= 1.0 / float(NUM_OF_LIGHTS);
      // color += 1.0 * vec3(pow(specAll, 8.0));

      // color += 0.01 * isMaterialSmooth(t.y, 1.0) * matCap(reflect(rayDirection, nor));

      vec3 reflectColor = vec3(0);
      vec3 reflectionRd = reflect(rayDirection, nor);
      reflectColor += 1.0 * reflection(pos, reflectionRd) * isMaterialSmooth(t.y, 1.0);
      color += reflectColor;

      // color += 0.015625 * dispersionStep1(nor, rayDirection, n2, n1);
      color += 0.25 * dispersion(nor, rayDirection, n2, n1);
      // color += 0.005 + 0.005 * sin(TWO_PI * (dot(nor, -rayDirection) + vec3(0, 0.33, 0.67)));

      // Fog
      color = mix(background, color, (fogMaxDistance - t.x) / fogMaxDistance);
      color = mix(background, color, exp(-t.x * 0.005));

      // Inner Glow
      // color += 0.5 * innerGlow(5.0 * t.w);

      // Post process
      // vec3 colorBefore = color;
      // colorMap(color);
      // color = mix(color, colorBefore, 0.5);

      // Debugging
      #ifdef debugMapCalls
      color = vec3(t.z / float(maxSteps));
      #endif

      #ifdef debugMapMaxed
      if (t.z / float(maxSteps) > 0.9) {
        color = vec3(1., 0., 1.);
      }
      #endif

      color = saturate(color);
      return vec4(color, 1.);
    } else {
      vec4 color = vec4(0.);
      if (!BLOOM) {
        color.a = 1.0;
      }

      // Radial Gradient
      // color = mix(vec4(theColor(uv), 1.0), vec4(background, 1), pow(length(uv) * 1.7, 8.0));

      // Glow
      vec3 glowColor = vec3(1);
      float i = 1. - 0.9 * pow(clamp(t.z / float(maxSteps), 0., 1.), 9.0);
      color = mix(vec4(glowColor, 1.0), color, i);

      return color;
    }
}

vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  vec4 t = march(ro, rd);
  return shade(ro, rd, t, uv);
}

void main() {
    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;
    background = getBackground(uv);

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec2 R = resolution * 2.;

    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 rd = getRayDirection(vec2(
                  float(x) / R.y + uv.x,
                  float(y) / R.y + uv.y),
                  projectionMatrix);
            rd = (vec4(rd, 1.) * cameraMatrix).xyz;
            rd = normalize(rd);
            color += saturate(sample(ro, rd, uv));
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(uv, projectionMatrix);
    rd = (vec4(rd, 1.) * cameraMatrix).xyz;
    rd = normalize(rd);
    gl_FragColor = saturate(sample(ro, rd, uv));
    #endif

    // gamma
    gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(0.454545));

    // Gradient effect
    // float brightness = length(gl_FragColor.rgb);
    // vec2 angle = normalize(vec2(1.0, 1.0));
    // gl_FragColor.rgb *= mix(
    //   vec3(1),
    //   mix(
    //     #330000,
    //     #00aaaa,
    //     saturate(-0.25 + dot(angle, uv.xy)))
    //   , 0.22);
    // gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(1.0 - 0.3 * brightness));
    // gl_FragColor.rgb *= 1.1;

    // Go to white as it gets brighter
    // float brightness = length(gl_FragColor.rgb);
    // gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(1.0 - 0.4 * brightness));

    // 'Film' Noise
    gl_FragColor.rgb += 0.015 * (cnoise2((500. + 60.1 * time) * uv + sin(uv + time)) + cnoise2((500. + 300.0 * time) * uv + 253.5));
}
