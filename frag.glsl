#extension GL_OES_standard_derivatives : enable
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)
#define saturate(x) clamp(x, 0.0, 1.0)

// #define debugMapCalls
// #define debugMapMaxed
// #define SS 2
// #define ORTHO 1

// @TODO Why is dispersion shitty on lighter backgrounds? I can see it blowing
// out, but it seems more than it is just screened or overlayed by the
// background instead of correctly fused into it.
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
uniform sampler2D textTex;

uniform float angle1C;
uniform float angle2C;
uniform float angle3C;

// KIFS
uniform mat4 kifsM;
uniform mat4 kifsM2;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 2048
#define maxDistance 100.0
#define fogMaxDistance 70.0

#define slowTime time * 0.2
// v3
// #define slowTime time * 0.06666667

vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);
const float totalT = 16.0;
float modT = mod(time, totalT);
float norT = modT / totalT;
float cosT = TWO_PI / totalT * modT;
const float edge = 0.01;
const float thickness = 0.75;

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: cnoise4 = require(glsl-noise/classic/4d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
//#pragma glslify: pnoise3 = require(glsl-noise/periodic/3d)
#pragma glslify: vmax = require(./hg_sdf/vmax)

#define combine(v1, v2, t, p) mix(v1, v2, t/p)

#pragma glslify: rotationMatrix = require(./rotation-matrix3)

float ncnoise2(in vec2 x) {
  return smoothstep(-1.00, 1.00, cnoise2(x));
}
float ncnoise3(in vec3 x) {
  return smoothstep(-1.00, 1.00, cnoise3(x));
}

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
  const float a = 0.823;
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
float vfbmWarp (vec3 p, out vec3 q, out vec3 s) {
  const float scale = 4.0;

  q = vec3(
        vfbm4(p + vec3(0.0, 0.0, 0.0)),
        vfbm4(p + vec3(3.2, 34.5, .234)),
        vfbm4(p + vec3(7.0, 2.9, -2.42)));

  s = vec3(
        vfbm4(p + scale * q + vec3(23.9, 234.0, -193.0)),
        vfbm4(p + scale * q + vec3(3.2, 852.0, 23.42)),
        vfbm4(p + scale * q + vec3(7.0, -232.0, -2.42)));

  return vfbm6(p + scale * s);
}
float vfbmWarp (vec3 p, out vec3 q) {
  vec3 s = vec3(0);
  vec3 r = vec3(0);
  return vfbmWarp(p, q, r);
}
float vfbmWarp (vec3 p) {
  vec3 q = vec3(0);
  vec3 s = vec3(0);
  vec3 r = vec3(0);
  return vfbmWarp(p, q, r);
}
float vfbmWarp (vec2 p, out vec2 q, out vec2 s, vec2 r) {
  const float scale = 4.0;
  const float angle = 0.01 * PI;
  const float si = sin(angle);
  const float c = cos(angle);
  const mat2 rot = mat2(c, si, -si, c);

  q = vec2(
        vfbm4(p + vec2(0.0, 0.0)),
        vfbm4(p + vec2(3.2, 34.5)));
  q *= rot;

  s = vec2(
        vfbm4(p + scale * q + vec2(23.9, 234.0)),
        vfbm4(p + scale * q + vec2(7.0, -232.0)));
  s *= rot;

  // r = vec2(
  //       vfbm4(p + scale * s + vec2(23.9, 234.0)),
  //       vfbm4(p + scale * s + vec2(7.0, -232.0)));
  // r *= rot;

  return vfbm6(p + scale * s);
}
float vfbmWarp (vec2 p) {
  vec2 q = vec2(0);
  vec2 s = vec2(0);
  vec2 r = vec2(0);

  return vfbmWarp(p, q, s, r);
}

// vec3 hsv(vec3 c);
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
  vec2 q = vec2(length8(p.xz)-t.x,p.y);
  return length8(q)-t.y;
}

float sdEllipsoid( in vec3 p, in vec3 r ) {
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
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

#define Iterations 5
// #pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
// #pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM, Iterations=Iterations)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
// #pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

// #pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=3, kifsM=kifsM, trapCalc=trapCalc)
// 
#pragma glslify: fold = require(./folds)
#pragma glslify: foldNd = require(./foldNd)
#pragma glslify: twist = require(./twist)

void opCheapBend (inout vec3 p, float a) {
    float c = cos(a*p.y);
    float s = sin(a*p.y);
    mat2 m = mat2(c,-s,s,c);
    p = vec3(m*p.xy,p.z);
}
void opCheapBend (inout vec3 p) {
  opCheapBend(p, 20.0);
}

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
vec3 dSMax (vec3 d1, vec3 d2, in float r) {
  float h = saturate(0.5 + 0.5 * (d1.x - d2.x) / r);
  float d = mix(d2.x, d1.x, h) + h * ( 1.0 - h ) * r;
  return vec3(d, (d1.x < d2.x) ? d1.yz : d2.yz);
}

vec3 dMax (vec3 d1, vec3 d2) {
  return (d1.x > d2.x) ? d1 : d2;
}

mat3 globalRot;
mat3 globalLRot;

#pragma glslify: rotMat2 = require(./rotation-matrix2)

// IQ
float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}
float sdCappedCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// p as usual, e exponent (p in the paper), r radius or something like that
#pragma glslify: octahedral = require(./model/octahedral)
#pragma glslify: dodecahedral = require(./model/dodecahedral)
#pragma glslify: icosahedral = require(./model/icosahedral)

#pragma glslify: sdTriPrism = require(./model/tri-prism)

bool isMaterial( float m, float goal ) {
  return m < goal + 1. && m > goal - .1;
}
float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

#pragma glslify: pModInterval1 = require(./hg_sdf/p-mod-interval1)
#pragma glslify: pMod1 = require(./hg_sdf/p-mod1.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: pMod3 = require(./hg_sdf/p-mod3.glsl)
#pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: quad = require(glsl-easings/quintic-in-out)
// #pragma glslify: cub = require(glsl-easings/cubic-in-out)
#pragma glslify: bounceIn = require(glsl-easings/bounce-in)
#pragma glslify: bounceOut = require(glsl-easings/bounce-out)
#pragma glslify: cubic = require(glsl-easings/cubic-in-out)
#pragma glslify: cubicOut = require(glsl-easings/cubic-out)
#pragma glslify: cubicIn = require(glsl-easings/cubic-in)
#pragma glslify: circ = require(glsl-easings/circular-in-out)
#pragma glslify: circIn = require(glsl-easings/circular-in)
#pragma glslify: circOut = require(glsl-easings/circular-out)
#pragma glslify: expo = require(glsl-easings/exponential-in-out)
#pragma glslify: expoOut = require(glsl-easings/exponential-out)
#pragma glslify: elastic = require(glsl-easings/elastic-in-out)
#pragma glslify: sine = require(glsl-easings/sine-in-out)
#pragma glslify: quart = require(glsl-easings/quadratic-in-out)
#pragma glslify: quint = require(glsl-easings/quintic-in-out)
#pragma glslify: quintIn = require(glsl-easings/quintic-in)
#pragma glslify: quintOut = require(glsl-easings/quintic-out)
// #pragma glslify: elasticInOut = require(glsl-easings/elastic-in-out)
#pragma glslify: elasticOut = require(glsl-easings/elastic-out)
// #pragma glslify: elasticIn = require(glsl-easings/elastic-in)
#pragma glslify: voronoi = require(./voronoi)
// #pragma glslify: band = require(./band-filter)

#pragma glslify: tetrahedron = require(./model/tetrahedron)
#pragma glslify: cellular = require(./cellular-tile)

// Starts at 0.5 goes towards 1.0
float nsin (in float t) {
  return 0.5 + 0.5 * sin(TWO_PI * t);
}
vec2 nsin (in vec2 t) {
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

#define jTrap 14
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

vec3 travel (in float t) {
  const float scale = 1.0;
  t *= scale;
  return vec3(mod(t, scale) - scale * 0.5, 0, 0);
}
// Source: https://www.shadertoy.com/view/MdcXzn
const float X_REPEAT_DIST = 0.90;
const float Z_REPEAT_DIST = 1.05;
vec3 DF_repeatHex(vec3 p)
{
    //Repetition
    float xRepeatDist = X_REPEAT_DIST;
    float zRepeatDist = Z_REPEAT_DIST*0.5;
    float latticeX = (fract(p.x/xRepeatDist+0.5)-0.5)*xRepeatDist;
    float latticeY = (fract(p.z/zRepeatDist+0.5)-0.5)*zRepeatDist;
    vec2 anchorPosXZ = p.xz-vec2(latticeX,latticeY);
    p.x = latticeX; //Cyclic coords.
    p.z = latticeY;

    return p;
}

vec3 mPos = vec3(0);
mat3 mRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);

vec3 crystal (vec3 q, vec3 offset) {
  vec3 d = vec3(maxDistance, 0, 0);

  vec3 base = vec3(sdHexPrism(q, vec2(0.1, 0.5)), 0, 0);
  d = dMin(d, base);

  // Cap
  vec3 dQ = q.xzy * vec3(1.5, 0.5, 1.6);
  vec3 cap = vec3(dodecahedral(dQ, 62., 0.3));
  cap -= 0.05 * cellular(1.2 * q + offset);
  d = dMin(d, cap);

  // Inclusions
  d.x = max(d.x, -0.2 * cellular(q));

  return d;
}

vec3 crystal(vec3 q) {
  return crystal(q, vec3(0));
}

vec3 rowOfBoxes (in vec3 q, in float size, in float r) {
  vec3 d = vec3(maxDistance, 0, 0);

  vec3 b;
  b = vec3(sdBox(q + 2. * (size + r) * vec3( 0, 0,  0), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3( 0, 0,  1), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3( 1, 0,  1), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3( 1, 0,  0), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3( 1, 0, -1), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3( 0, 0, -1), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3(-1, 0, -1), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3(-1, 0,  0), vec3(size)), 0, 0);
  d = dMin(d, b);

  b = vec3(sdBox(q - 2. * (size + r) * vec3(-1, 0,  1), vec3(size)), 0, 0);
  d = dMin(d, b);

  return d;
}

// Return value is (distance, material, orbit trap)
vec3 map (in vec3 p, in float dT) {
  vec3 d = vec3(maxDistance, 0, 0);

  vec3 q = p;

  q += 0.000500 * snoise3(73. * q.yzx);

  mPos = q;

  /*
  const float size = 0.1;
  // Floor
  vec3 f = vec3(sdPlane(q, vec4(0, 1, 0, 0)), 0, 0);
  d = dMin(d, f);

  // Boxes
  vec2 c = pMod2(q.xz, vec2(size));

  const float vOffset = 0.05;
  q.y += 0.75 * vOffset * sin(length(c) - cosT) + vOffset;

  const float lateralScale = size * (0.5 - 0.3);
  vec3 b = vec3(sdBox(q, vec3(lateralScale, 0.5 * vOffset, lateralScale)), 0, 0);
  d = dSMin(d, b, 0.05 * scale);

  d.x *= 0.60;
  */

  // Box
  const float size = 1.2;
  vec3 b = vec3(sdBox(q + vec3(0, 0, size), vec3(size)), 0, 0);
  d = dMin(d, b);

  // --- Cutout ---
  float cutoutT = 1. - expo(2. * norT - step(0.5, norT) * 4. * (norT - 0.5));

  q += 0.2 * snoise3(3. * q.yzx - slowTime);
  float n = length(q) - 0.8 * smoothstep(0., 0.9, sin(cosT));
  n *= -1.;
  d.x = max(d.x, n);
  /*
  // Center Sphere
  const float sRadius = 0.4;
  q.z -= sRadius * cutoutT;
  vec3 c = vec3(length(q) - sRadius, 0, 0);
  c.x *= -1.; // Invert
  d = dMax(d, c);

  // Bar
  const float bHeight = sRadius * 2.25;
  const float bWidth = 0.05;
  q.xy *= rotMat2(PI * 0.25);
  vec3 bar = vec3(sdBox(q - vec3(0, 0, 0.65 * bWidth), vec3(bWidth, bHeight, bWidth)), 0, 0);

  // ... Buffer Circle
  const float bSRadius = sRadius * 1.1;;
  c = vec3(length(q) - bSRadius, 0, 0);
  c.x *= -1.; // Invert
  bar = dMax(bar, c);

  bar.x *= -1.; // Invert
  d = dMax(d, bar);
  */

  d.x *= 0.2;

  return d;
}
vec3 map (in vec3 p) {
  return map(p, 0.);
}

vec4 march (in vec3 rayOrigin, in vec3 rayDirection, in float deltaT) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  const float deltaTDelta = 0.000;
// #define OVERSTEP 1
#ifdef OVERSTEP
  // Source: https://www.shadertoy.com/view/MdcXzn
  const int halfMax = (maxSteps / 2);
  for( int i = 0; i < halfMax; i++ ) {
       vec3 d = map(rayOrigin + rayDirection * t, deltaT);

       if( d.x < epsilon * 100.0 ) break;
       t += d.x;
       if (t > maxDistance) break;
      deltaT += deltaTDelta;
   }

   t -= Z_REPEAT_DIST * 0.5;

   for( int i = 0; i < halfMax; i++ ) {
       vec3 d = map(rayOrigin + rayDirection * t);

       if( d.x<epsilon ) return vec4(t + d.x, d.y, float(i), d.z);

       t += min(d.x, Z_REPEAT_DIST * 0.2);
      deltaT += deltaTDelta;
   }
#else
  for (int i = 0; i < maxSteps; i++) {
    vec3 d = map(rayOrigin + rayDirection * t, deltaT);
    if (d.x < epsilon) return vec4(t + d.x, d.y, float(i), d.z);
    t += d.x;
    maxI = float(i);
    trap = d.z;
    if (t > maxDistance) break;
    deltaT += deltaTDelta;
  }
#endif
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

// #pragma glslify: checker = require(glsl-checker)
#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: rgb2hsv = require(./rgb2hsv.glsl)
// #pragma glslify: debugColor = require(./debug-color-clip)
#pragma glslify: hsb2rgb = require(./color-map/hsb2rgb)

const float n1 = 1.0;
const float n2 = 1.5;
const float amount = 0.08;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float spread = 1.0; // saturate(1.0 - 2.0 * dot(-rd, gNor));
  // float n = smoothstep(0.75, 1.0, sin(250.0 * rd.x + 0.01 * noise(433.0 * rd)));

  float startPoint = 0.1;

  vec3 spaceScaling = vec3(0.734, 1.14, 0.2);
  float n = ncnoise3(spaceScaling * rd + startPoint);
  n = smoothstep(0.0, 0.80, n);

  // vec3 spaceScaling = vec3(0.5);
  // float n = vfbmWarp(spaceScaling * rd + startPoint);
  // n = smoothstep(0.65, 0.85, n);

  // float n = smoothstep(0.9, 1.0, sin(TWO_PI * (dot(vec2(8), rd.xz) + 2.0 * cnoise3(1.5 * rd)) + time));

  // float n = cnoise3(0.5 * rd);
  // n = smoothstep(-0.1, 0.9, n);

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

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=amount, time=time)

#ifndef dispersionMap
#define dispersionMap map
#endif

float dispersionMarch (in vec3 rayDirection) {
  vec3 rayOrigin = gPos + -gNor * 0.01;
  rayDirection = normalize(rayDirection);

  float t = 0.0;

  for (int i = 0; i < 20; i++) {
    float d = dispersionMap(rayOrigin + rayDirection * t).x;
    if (d >= 0.0) break;
    d = min(d, -0.0001);

    t += abs(d);
  }
  return t;
}

vec3 secondRefraction (in vec3 rd, in float ior) {
  float d = 0.0;

  float angle = ncnoise3(rd);

  rd *= rotationMatrix(vec3(1, 4, 2), smoothstep(0.5, 0.6, angle));

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

  vec3 reflectionPoint = gPos - gNor * 0.1 + rd * d;
  vec3 reflectionPointNor = getNormal2(reflectionPoint, 0.001);
  dNor = reflectionPointNor;
  reflectionPointNor = normalize(reflectionPointNor);

  vec3 disp = min(1.5, 1.0 / d) * dispersion(reflectionPointNor, rd, n2, n1);

  return disp;
}

#pragma glslify: dispersionStep1 = require(./glsl-dispersion, scene=secondRefraction, amount=amount, time=time)

// #pragma glslify: rainbow = require(./color-map/rainbow)

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  vec3 color = vec3(0.80);
  color += 0.175 * snoise3(1345. * vec3(1, 1, 0.01) * pos);

  return color;
}

#pragma glslify: reflection = require(./reflection, getNormal=getNormal2, diffuseColor=baseColor, map=map, maxDistance=maxDistance, epsilon=epsilon, maxSteps=512, getBackground=getBackground)

const vec3 glowColor = pow(#ED4F2C, vec3(2.2));

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)

vec4 shade ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    gPos = pos;

    // Declare lights
    struct light {
      vec3 position;
      vec3 color;
      float intensity;
    };
    const int NUM_OF_LIGHTS = 3;
    const float repNUM_OF_LIGHTS = 0.33333;
    light lights[NUM_OF_LIGHTS];

    // vec2 lightPosRef = vec2(0.95, 0.2);
    // mat2 lightPosRefInc = rotMat2(TWO_PI * repNUM_OF_LIGHTS);
    // lightPosRef *= rotMat2(TWO_PI * mod(time * 0.1, 1.));

    // for (int i = 0; i < NUM_OF_LIGHTS; i++) {
    //   vec3 lightColor = hsb2rgb(vec3(float(i) * 1.1 * repNUM_OF_LIGHTS, 1., 1));
    //   float greenish = 0.0; // dot(normalize(lightColor), #00FF00);
    //   lights[i] = light(vec3(lightPosRef, 0.25), lightColor, mix(1.8, 0.8, greenish));
    //   lightPosRef *= lightPosRefInc;
    // }

    const float lightPosScale = 0.35;
    // lights[0] = light(vec3(0, 0.2, 1.0), #FFFFFF, 1.0);
    // lights[0] = light(vec3(0, 0.2, 1.0), 0.5 + 0.5 * cos(TWO_PI * (pos + vec3(0, 0.33, 0.67))), 1.0);
    // lights[1] = light(vec3(0.4, 0.4, 1.0), 0.5 + 0.5 * cos(TWO_PI * (lightPosScale * pos + vec3(0, 0.1, 0.2))), 1.0);
    // lights[2] = light(vec3(0.4, 0, 1.0), 0.5 + 0.5 * cos(TWO_PI * (lightPosScale * pos + vec3(0.2, 0.3, 0.4))), 1.0);
    lights[0] = light(vec3(0, 0.2, 1.0), #FFFFFF, 1.0);
    lights[1] = light(vec3(0.4, 0.4, 1.0), #FF7777, 1.0);
    lights[2] = light(vec3(0.4, 0, 1.0), #77FFFF, 1.0);

    if (t.x>0.) {
      vec3 color = vec3(0.0);

      // Material Types
      float isBlack = isMaterialSmooth(t.y, 1.);
      float isFloor = isMaterialSmooth(t.y, 0.);
      float isNeon = 1. - isBlack;

      vec3 nor = getNormal2(pos, 0.0001 * t.x);
      // float bumpsScale = 0.01;
      // float bumpIntensity = 1.0;
      // nor += bumpIntensity * vec3(
      //     cnoise3(bumpsScale * 490.0 * mPos),
      //     cnoise3(bumpsScale * 670.0 * mPos + 234.634),
      //     cnoise3(bumpsScale * 310.0 * mPos + 23.4634));
      // nor = normalize(nor);
      gNor = nor;

      vec3 ref = reflect(rayDirection, nor);
      ref = normalize(ref);

      gRd = rayDirection;

      // Basic Diffusion
      vec3 diffuseColor = baseColor(pos, nor, rayDirection, t.y, t.w);

      float occ = calcAO(pos, nor);
      float amb = saturate(0.5 + 0.5 * nor.y);
      float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 0.5;
      float specCo = 0.50;

      float specAll = 0.0;

      vec3 directLighting = vec3(0);
      for (int i = 0; i < NUM_OF_LIGHTS; i++) {
        vec3 lightPos = lights[i].position;
        float diffMin = 0.6;
        float dif = max(diffMin, diffuse(nor, normalize(lightPos)));
        float spec = pow(clamp( dot(ref, normalize(lightPos)), 0., 1. ), 256.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        float shadowMin = 0.3;
        float sha = max(shadowMin, softshadow(pos, normalize(lightPos), 0.001, 4.75));
        dif *= sha;

        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        lin += specCo * spec * (1. - fre);
        specAll += specCo * spec * (1. - fre);

        // Ambient
        lin += 0.7 * amb * diffuseColor;

        float distIntensity = 1.0; // lights[i].intensity / pow(length(lightPos - gPos), 2.0);
        color +=
          saturate((dif * distIntensity) * lights[i].color * diffuseColor)
          + saturate(lights[i].intensity * mix(lights[i].color, vec3(1), 0.1) * lin * mix(diffuseColor, vec3(1), 0.4));

        vec3 fromLight = rayOrigin - lightPos;
        float lightMasked = 1. - smoothstep(t.x, t.x + 0.001, length(fromLight));
        float lightAngle = pow(max(0., dot(-rayDirection, normalize(fromLight))), 512.0);
        // directLighting +=
        //     lightMasked
        //   * mix(lights[i].color, vec3(1), 0.9 * lightAngle)
        //   * lightAngle;
      }

      color *= 1.0 / float(NUM_OF_LIGHTS);
      color += 1.0 * vec3(pow(specAll, 8.0));

      // vec3 reflectColor = vec3(0);
      // vec3 reflectionRd = reflect(rayDirection, nor);
      // reflectColor += 0.3 * reflection(pos, reflectionRd);
      // color += reflectColor;

      // vec3 refractColor = vec3(0);
      // vec3 refractionRd = refract(rayDirection, nor, 1.5);
      // refractColor += textures(refractionRd);
      // color += refractColor;

      // vec3 dispersionColor = dispersionStep1(nor, rayDirection, n2, n1);
      // vec3 dispersionColor = dispersion(nor, rayDirection, n2, n1);
      // color += 1.0 * dispersionColor;
      // color = mix(color, color + dispersionColor, ncnoise3(1.5 * pos));
      // color = pow(color, vec3(1.2));

      // Fog
      float d = max(0.0, t.x);
      color = mix(background, color, saturate((fogMaxDistance - d) / fogMaxDistance));
      color *= exp(-d * 0.005);

      // color += directLighting * exp(-d * 0.0005);

      // Inner Glow
      // color += 0.5 * innerGlow(5.0 * t.w);

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

      for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
        vec3 lightPos = lights[i].position;
        vec3 fromLight = rayOrigin - lightPos;
        float lightMasked = 1. - smoothstep(t.x, t.x + 0.001, length(fromLight));
        float lightAngle = pow(dot(-rayDirection, normalize(fromLight)), 512.0);
        color.rgb += lightMasked * mix(lights[i].color, vec3(1), lightAngle) * pow(dot(-rayDirection, normalize(fromLight)), 512.0);
      }

      // Radial Gradient
      // color = mix(vec4(vec3(0), 1.0), vec4(background, 1), saturate(pow((length(uv) - 0.25) * 1.6, 0.3)));

      // Glow
      // float i = saturate(t.z / (0.91 * float(maxSteps)));
      // vec3 glowColor = pow(#D93741, vec3(2.2));
      // color = mix(color, vec4(glowColor, 1.0), i);

      return color;
    }
}

void thrshld (inout float v, in float level, in float edge) {
  v = smoothstep(level, level + edge, v);
}

void bnd (inout float v, in float start, in float end, in float eps) {
  float v1 = v;
  thrshld(v1, start, eps);
  float v2 = v;
  thrshld(v2, end, -eps);

  v = v1 * v2;
}

float crcl (in vec2 uv, float r) {
  return length(uv) - r;
}

float sqr (in vec2 uv, float r) {
  vec2 absUv = abs(uv);
  float l = max(absUv.x, absUv.y);
  return l - r;
}

vec2 hash( vec2 p  ) {
  p = vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3)));
  return fract(sin(p)*18.5453);
}

float noiseDots (vec2 uv, float size) {
  vec2 q = uv;
  vec2 c = pMod2(q, vec2(size));

  float r = size * 0.2;
  r += 0.35 * r * cnoise2(c + 5.0 * q / size);

  return 1. - smoothstep(0., 0.002, length(q) - r);
}

vec4 tile (in vec2 uv) {
  const float scale = 0.25;
  uv *= scale;
  return vec4(
      step(0.125, cnoise2(scale * 9.0 * uv +   0.0)),
      step(0.125, cnoise2(scale * 9.4 * uv + 110.4)),
      step(0.125, cnoise2(scale * 9.4 * uv + 813.1)),
      step(0.125, cnoise2(scale * 9.4 * uv - 310.0)));
}

vec4 neighborsTile (in vec2 uv) {
  return vec4(
      tile(uv - vec2(1, 0)).z,
      tile(uv - vec2(0, 1)).w,
      tile(uv + vec2(1, 0)).x,
      tile(uv + vec2(0, 1)).y);
}

float edgeDist (in vec2 q, in float edge) {
  float d = 100.0;

  float capMask = smoothstep(0.0, edge, q.x);
  float axisD = abs(q.y);
  axisD = mix(100.0, axisD, capMask);
  d = min(d, axisD);

  float roundMask = length(q);
  d = min(d, roundMask);

  return d;
}
float edgeBand (in vec2 q, in float thickness, in float edge) {
  return smoothstep(thickness + edge, thickness, edgeDist(q, edge));
}

vec3 tileColor (vec2 q, vec4 sides, in float size, in float colorOffset) {
  vec2 absQ = abs(q);
  vec3 primaryColor = vec3(0);
  float mask = 0.;

  const float edge = 0.0001;
  const float thickness = 0.333333;

  float d = 100.0;
  // Left
  vec2 inputQ = q.xy * vec2(-1, 1);
  float edgeMask = edgeBand(inputQ, thickness * size, edge * size);
  mask += sides.x * edgeMask;
  float spaceD = edgeDist(inputQ, edge * size);
  float edgeD = mix(100., spaceD, sides.x * edgeMask);
  d = min(d, edgeD);
  // Right
  inputQ = q.xy * vec2( 1, 1);
  edgeMask = edgeBand(inputQ, thickness * size, edge * size);
  spaceD = edgeDist(inputQ, edge * size);
  edgeD = mix(100., spaceD, sides.z * edgeMask);
  d = min(d, edgeD);
  mask += sides.z * edgeMask;

  // Top
  inputQ = q.yx * vec2(-1, 1);
  edgeMask = edgeBand(inputQ, thickness * size, edge * size);
  spaceD = edgeDist(inputQ, edge * size);
  edgeD = mix(100., spaceD, sides.y * edgeMask);
  d = min(d, edgeD);
  mask += sides.y * edgeMask;
  // Bottom
  inputQ = q.yx * vec2( 1, 1);
  edgeMask = edgeBand(inputQ, thickness * size, edge * size);
  spaceD = edgeDist(inputQ, edge * size);
  edgeD = mix(100., spaceD, sides.w * edgeMask);
  d = min(d, edgeD);
  mask += sides.w * edgeMask;

  // Colors
  float colorI = d / (thickness * size);

  float numColors = 4.0;
  // colorI = floor(colorI * numColors) / numColors;
  colorI += colorOffset;
  primaryColor = vec3(0.51 + 0.5 * cos(TWO_PI * colorI * 2.5 - 4.0 * cosT));

  return saturate(mask) * primaryColor;
}

vec3 grid (in vec2 uv, in float size, in float colorOffset) {
  vec3 color = vec3(0);

  vec2 q = uv;
  vec2 c = pMod2(q, vec2(size));

  // Show Borders
  // vec2 absQ = abs(q);
  // float borderD = max(absQ.x, absQ.y) - 0.4 * size;
  // color = mix(color, vec3(0, 0, 1), smoothstep(0., 0.001, borderD));

  vec4 sides = tile(c);
  float isPrimary = floor(mod(dot(c, vec2(1)), 2.));

  vec3 primaryColor = tileColor(q, sides, size, colorOffset);
  vec3 secondaryColor = tileColor(q, neighborsTile(c), size, colorOffset);

  color = mix(color, primaryColor, isPrimary * smoothstep(0., 0.01, length(primaryColor)));
  color = mix(secondaryColor, color, isPrimary);

  return color;
}

vec2 cellBoxes (in vec2 q, in vec2 c, in float cellSize) {

#define MOVING_SQR 1
#ifdef MOVING_SQR
  float moveRate = 0.05;
#else
  float moveRate = 0.;
#endif

  vec2 d = vec2(200, 1);
  for (int i = 0; i < 6; i++) {
    vec2 qB = q + 1.3 * cellSize * vec2(
        noise(1.32345 * vec2(i) + c + 2.0 * moveRate * time + 0.),
        noise(2.92345 * vec2(i) + c + 1.0 * moveRate * time + 123.));
    float r = 0.2 * cellSize;
    vec2 absQ = abs(qB);
    float sqrD = max(absQ.x, absQ.y);
    float m = smoothstep(0., 0.01, sqrD - r);
    float n = smoothstep(0.001, 0.0, sqrD - 0.95 * r);

    vec2 crossQ = qB * rotMat2(PI * 0.25);
    absQ = abs(crossQ);
    sqrD = max(absQ.x, absQ.y);
    float axisD = min(absQ.x, absQ.y) - 0.05 * r;
    float axis = smoothstep(0.001, 0.0, axisD);
    float crossR = r * 0.5;
    float cross = axis * smoothstep(0.001, 0., sqrD - crossR);
    cross = 1. - cross;
    n *= cross;

    float z = 1.; // noise(vec2(i));

    vec2 b = vec2(m * 100. + z, n);
    d = dMin(d, b);
  }

  return d;
}

vec3 two_dimensional (in vec2 uv, in float generalT) {
  vec3 color = vec3(1);

  vec2 q = uv;

  const float gridScale = 20.0;
  const float invGridScale = 0.05;
  const float stretch = 0.20;

  q += vec2(
      2.0 *invGridScale,
      2.0 *invGridScale * stretch)
    * smoothstep(0.3, 1., sin(cosT - length(q)));

  const float size = 0.01;
  const float halfsize = size * 0.5;

  float inShape = 0.;

  vec2 fontUV = 0.5 * (uv + 1.);
  fontUV.y = 1. - fontUV.y;

  inShape = texture2D(textTex, fontUV).r;
  inShape = step(0.3, inShape);

  // Grid
  q *= gridScale;

  vec2 axis = mix(vec2(1, stretch), vec2(stretch, 1), inShape);
  vec2 c = pMod2(q, axis);

  float v = mod(dot(c, vec2(1)), 2.);

  v = saturate(v);
  const vec3 foregroundColor = vec3(1);
  color = mix(background, foregroundColor, v);

  return color;
}

vec3 two_dimensional (in vec2 uv) {
  return two_dimensional(uv, modT);
}

vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  return vec4(two_dimensional(uv), 1);
  vec4 t = march(ro, rd, 0.20);
  return shade(ro, rd, t, uv);
}

void main() {
    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;
    background = getBackground(uv);

    float gRAngle = TWO_PI * mod(time, totalT) / totalT;
    float gRc = cos(gRAngle);
    float gRs = sin(gRAngle);
    globalRot = mat3(
      gRc, 0.0, -gRs,
      0.0, 1.0,  0.0,
      gRs, 0.0,  gRc);
    float glRc = cos(-gRAngle);
    float glRs = sin(-gRAngle);
    globalLRot = mat3(
      glRc, 0.0, -glRs,
      0.0, 1.0,  0.0,
      glRs, 0.0,  glRc);

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec2 R = resolution * 2.;

    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
#ifdef ORTHO
            vec3 rd = vec3(0, 0, -1);
#else
            vec3 rd = getRayDirection(vec2(
                  float(x) / R.y + uv.x,
                  float(y) / R.y + uv.y),
                  projectionMatrix);
#endif
            rd = (vec4(rd, 1.) * cameraMatrix).xyz;
            rd = normalize(rd);
#ifdef ORTHO
            vec2 ndc = (gl_FragCoord.xy + 0.0 * vec2(x, y)) / resolution.xy * 2.0 - 1.0;
            float w = 2.0;
            float h = w / (resolution.x / resolution.y);
            ro += (vec4(
                  ndc * vec2(w * 0.5, h * 0.5),
                  0, 1) * cameraMatrix).xyz;
#endif
            color += saturate(sample(ro, rd, uv));
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
#ifdef ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = getRayDirection(uv, projectionMatrix);
#endif
    rd = (vec4(rd, 1.) * cameraMatrix).xyz;
    rd = normalize(rd);
#ifdef ORTHO
    vec2 ndc = gl_FragCoord.xy / resolution.xy * 2.0 - 1.0;
    float w = 2.0;
    float h = w / (resolution.x / resolution.y);
    ro += (vec4(
          ndc * vec2(w * 0.5, h * 0.5),
          0, 1) * cameraMatrix).xyz;
#endif
    gl_FragColor = saturate(sample(ro, rd, uv));
    #endif

    // gamma
    gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(0.454545));

    // Gradient effect
    // float brightness = length(gl_FragColor.rgb);
    // vec2 angle = normalize(vec2(0.0, 1.0));
    // gl_FragColor.rgb *= mix(
    //   vec3(1),
    //   mix(
    //     #FF1111,
    //     #00aaaa,
    //     saturate(-0.25 + dot(angle, uv.xy)))
    //   , 0.3);
    // gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(1.0 - 0.3 * brightness));
    // gl_FragColor.rgb *= 1.1;

    // Go to white as it gets brighter
    // float brightness = length(gl_FragColor.rgb);
    // gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(1.0 - 0.4 * brightness));

    // vec2 absUV = abs(uv);
    // float vignette = smoothstep(0.65, 1.4, max(absUV.x, absUV.y));
    // vignette *= vignette;
    // gl_FragColor.a += vignette;
    // vignette = 1.0 - vignette;
    // gl_FragColor.rgb *= vec3(vignette);
}
