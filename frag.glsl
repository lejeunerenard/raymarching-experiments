#extension GL_OES_standard_derivatives : enable
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)
#define saturate(x) clamp(x, 0.0, 1.0)

// #define debugMapCalls
// #define debugMapMaxed
// #define SS 2
// #define ORTHO 1
// #define NO_MATERIALS 1

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

uniform vec3 colors1;
uniform vec3 colors2;

uniform float d;

// KIFS
uniform mat4 kifsM;
uniform mat4 kifsM2;
uniform float scale;
uniform vec3 offset;
uniform float rot;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 1024
#define maxDistance 10.0
#define fogMaxDistance 5.25

#define slowTime time * 0.2
// v3
// #define slowTime time * 0.06666667

vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);
const float totalT = 6.0;
float modT = mod(time, totalT);
float norT = modT / totalT;
float cosT = TWO_PI / totalT * modT;
const float edge = 0.0025;
const float thickness = 0.05;

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
  const float scale = 2.0;

  q = vec2(
        iqFBM(p + vec2(0.0, 0.0)),
        iqFBM(p + vec2(7.2, 34.5)));

  s = vec2(
        iqFBM(p + scale * q + vec2(93.9, 234.0)),
        iqFBM(p + scale * q + vec2(3.2, 123.0)));

  r = vec2(
        iqFBM(p + scale * s + vec2(23.9, 74.0)),
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
float sdBox( vec4 p, vec4 b ) {
  vec4 d = abs(p) - b;
  return min(max(d.x,max(d.y,max(d.z, d.w))),0.0) + length(max(d,0.0));
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
float sdTorus82( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length8(q)-t.y;
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

#define Iterations 8
#pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
#pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM, Iterations=Iterations)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
#pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=5, kifsM=kifsM, trapCalc=trapCalc)
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

vec2 dMax (vec2 d1, vec2 d2) {
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
#pragma glslify: bounce= require(glsl-easings/bounce-in-out)
#pragma glslify: cubic = require(glsl-easings/cubic-in-out)
#pragma glslify: cubicOut = require(glsl-easings/cubic-out)
#pragma glslify: cubicIn = require(glsl-easings/cubic-in)
#pragma glslify: circ = require(glsl-easings/circular-in-out)
#pragma glslify: circIn = require(glsl-easings/circular-in)
#pragma glslify: circOut = require(glsl-easings/circular-out)
#pragma glslify: expo = require(glsl-easings/exponential-in-out)
#pragma glslify: expoIn = require(glsl-easings/exponential-in)
#pragma glslify: expoOut = require(glsl-easings/exponential-out)
#pragma glslify: elastic = require(glsl-easings/elastic-in-out)
#pragma glslify: sine = require(glsl-easings/sine-in-out)
#pragma glslify: sineOut = require(glsl-easings/sine-out)
#pragma glslify: quart = require(glsl-easings/quadratic-in-out)
#pragma glslify: quint = require(glsl-easings/quintic-in-out)
#pragma glslify: quintIn = require(glsl-easings/quintic-in)
#pragma glslify: quintOut = require(glsl-easings/quintic-out)
// #pragma glslify: elasticInOut = require(glsl-easings/elastic-in-out)
#pragma glslify: elasticOut = require(glsl-easings/elastic-out)
// #pragma glslify: elasticIn = require(glsl-easings/elastic-in)

#pragma glslify: voronoi = require(./voronoi, edge=edge, thickness=thickness, mask=sqrMask)
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
vec3 mPos2 = vec3(0);
mat3 mRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);

mat3 rotOrtho (in float t) {
  const vec3 rotAxis = vec3(0, 1, 0);
  return rotationMatrix(rotAxis, 1.5 * PI * (0.5 + 0.5 * cos(t)));
}

float getLayer (in float t) {
  float l = 0.;
  l += step(0.3, t);
  l += step(0.6, t);

  return l;
}

// Create multiple copies of an object - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
vec2 opRepLim( in vec2 p, in float s, in vec2 lim ) {
  return p-s*clamp(floor(p/s + 0.5),-lim,lim);
}

const float height = 0.2;
const float size = 0.1;
vec3 map (in vec3 p, in float dT) {
  vec3 d = vec3(maxDistance, 0, 0);

  p *= rotationMatrix(vec3(0, 0, 1), -0.110 * PI);
  vec3 q = p;

  float t = mod(cosT + PI, TWO_PI);

  const float warpScale = 0.65;
  const float r = 0.50;

  // Grid polar
  q = vec3(
      atan(q.z, q.x),
      length(q.zx) - 1.25 * r,
      q.y);

  const float size = 0.0625 * 2. * r;
  q.yz *= rotMat2(cosT + q.x + 0.125 * PI * cos(cosT + q.x));
  q.yz = opRepLim(q.yz, size, vec2(1, 2));

  mPos = q;
  vec3 s = vec3(sdCylinder(q.yxz, vec3(vec2(0), 0.2 * size)), 0, 0);
  d = dMin(d, s);

  q = p;
  q.y += 0.125 * r * sin(cosT);
  vec3 planet = vec3(length(q) - 0.65 * r, 1, 0);
  planet.x -= 0.0125 * cellular(4. * q);
  d = dMin(d, planet);

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
vec3 getNormal2 (in vec3 p, in float eps, in float generalT) {
  vec2 e = vec2(1.,0.) * eps;
  return normalize(vec3(
    map(p + e.xyy, generalT).x - map(p - e.xyy, generalT).x,
    map(p + e.yxy, generalT).x - map(p - e.yxy, generalT).x,
    map(p + e.yyx, generalT).x - map(p - e.yyx, generalT).x));
}
vec3 getNormal2 (in vec3 p, in float eps) {
  return getNormal2(p, eps, 0.);
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
const float n2 = 1.50;
const float amount = 0.1;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float spread = 1.; // saturate(1.0 - 1.0 * dot(-rd, gNor));
  // float n = smoothstep(0., 1.0, sin(150.0 * rd.x + 0.01 * noise(433.0 * rd)));

  float startPoint = 0.0;

  vec3 spaceScaling = 0.2 * vec3(0.734, 1.14, 0.2);
  float n = ncnoise3(spaceScaling * rd + startPoint);
  n = smoothstep(0.0, 0.80, n);

  /* vec3 spaceScaling = vec3(0.8); */
  /* float n = vfbmWarp(spaceScaling * rd + startPoint); */
  /* n = smoothstep(0.525, 0.80, n); */

  /* vec3 spaceScaling = vec3(9.8); */
  /* float n = vfbm4(spaceScaling * rd + startPoint); */
  /* n = smoothstep(0.125, 0.85, n); */

  /* float n = smoothstep(0.9, 1.0, sin(TWO_PI * (dot(vec2(8), rd.xz) + 2.0 * cnoise3(1.5 * rd)) + time)); */

  /* float n = cnoise3(3.5 * rd); */
  /* n = smoothstep(-0.1, 0.9, n); */

  // float n = 0.6 + 0.4 * sin(dot(vec3(PI), sin(3.18 * rd + sin(1.38465 * rd.yzx))));

  color = vec3(n * spread);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd, in float ior) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=amount, time=time, norT=norT)

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

#pragma glslify: dispersionStep1 = require(./glsl-dispersion, scene=secondRefraction, amount=amount, time=time, norT=norT)

vec3 baseColor (in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap, in float t) {
  vec3 color = vec3(0.5);

  vec3 dI = vec3(0.2); // vec3(angle1C);
  dI += vec3(0.159155, vec2(0.1)) * mPos;
  dI += 0.1 * nor;
  dI += 0.2 * pow(dot(nor, -rd), 2.);
  dI += norT; // angle1C;
  color = 0.5 + 0.5 * cos( TWO_PI * (dI + vec3(0, 0.33, 0.67)) );
  color += 0.2 * (0.5 + 0.5 * cos( TWO_PI * (nor + vec3(0, 0.33, 0.67)) ));

  color = mix(color, vec3(0.1), isMaterialSmooth(m, 1.));

#ifdef NO_MATERIALS
  color = vec3(0.5);
#endif
  return color;
}

#pragma glslify: reflection = require(./reflection, getNormal=getNormal2, diffuseColor=baseColor, map=map, maxDistance=maxDistance, epsilon=epsilon, maxSteps=512, getBackground=getBackground)

const vec3 glowColor = pow(#ED4F2C, vec3(2.2));

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)

vec4 shade ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv, in float generalT ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    gPos = pos;

    // Declare lights
    struct light {
      vec3 position;
      vec3 color;
      float intensity;
    };

    const int NUM_OF_LIGHTS = 3;
    const float repNUM_OF_LIGHTS = 0.333333;

    light lights[NUM_OF_LIGHTS];

    // vec2 lightPosRef = vec2(0.95, 0.2);
    // mat2 lightPosRefInc = rotMat2(TWO_PI * repNUM_OF_LIGHTS);

    // for (int i = 0; i < NUM_OF_LIGHTS; i++) {
    //   vec3 lightColor = hsb2rgb(vec3(float(i) * 1.1 * repNUM_OF_LIGHTS, 1., 1));
    //   lights[i] = light(vec3(lightPosRef, 0.25), lightColor, 2.0);
    //   lightPosRef *= lightPosRefInc;
    // }

    // lights[0] = light(normalize(vec3(  0.15, 0.25, 1.0)), #FFFFFF, 1.0);
    lights[0] = light(vec3(-1.0, 0.5,  0.5), #FFFFFF, 1.0);
    lights[1] = light(vec3(-1.0, 1.0,  1.0), #FFFFFF, 1.0);
    lights[2] = light(vec3( 0.0, 0.0,  1.0), #FFFFFF, 1.0);

    float backgroundMask = 1.;
    // Allow anything in top right corner
    // backgroundMask = max(backgroundMask, smoothstep(0., edge, dot(uv, vec2(1)) + 0.05));

    if (t.x>0. && backgroundMask > 0.) {
      vec3 color = vec3(0.0);

      // Material Types
      float isFloor = isMaterialSmooth(t.y, 1.);

      // Normals
      vec3 nor = getNormal2(pos, 0.005 * t.x, generalT);
      // float bumpsScale = 4.75;
      // float bumpIntensity = 0.4;
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
      vec3 diffuseColor = baseColor(pos, nor, rayDirection, t.y, t.w, generalT);

      float occ = calcAO(pos, nor);
      float amb = saturate(0.5 + 0.5 * nor.y);
      float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 0.95;
      float specCo = 0.4;

      float specAll = 0.0;

      vec3 directLighting = vec3(0);
      for (int i = 0; i < NUM_OF_LIGHTS; i++) {
        vec3 lightPos = lights[i].position; // * globalLRot;
        const float diffMin = 0.7;
        float dif = max(diffMin, diffuse(nor, normalize(lightPos)));
        float spec = pow(clamp( dot(ref, normalize(lightPos)), 0., 1. ), 128.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        const float shadowMin = 0.90;
        float sha = max(shadowMin, softshadow(pos, normalize(lightPos), 0.001, 4.75));
        dif *= sha;

        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        specAll += specCo * spec; // * (1. - fre);

        // Ambient
        lin += 0.00 * amb * diffuseColor;
        // dif += 0.000 * amb;

        float distIntensity = 1.; // lights[i].intensity / pow(length(lightPos - gPos), 1.0);
        distIntensity = saturate(distIntensity);
        color +=
          saturate((dif * distIntensity) * lights[i].color * diffuseColor)
          + saturate(distIntensity * mix(lights[i].color, vec3(1), 0.1) * lin * mix(diffuseColor, vec3(1), 0.4));

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

      vec3 reflectColor = vec3(0);
      vec3 reflectionRd = reflect(rayDirection, nor);
      reflectColor += mix(0.025, 0.75, isMaterialSmooth(t.y, 1.)) * reflection(pos, reflectionRd);
      color += reflectColor;

      /* vec3 refractColor = vec3(0); */
      /* vec3 refractionRd = refract(rayDirection, nor, 1.5); */
      /* refractColor += 0.05 * textures(refractionRd); */
      /* color += refractColor; */

#ifndef NO_MATERIALS
      vec3 dispersionColor = dispersionStep1(nor, rayDirection, n2, n1);
      // dispersionColor = textures(rayDirection);
      // vec3 dispersionColor = dispersion(nor, rayDirection, n2, n1);

      dispersionColor *= isMaterialSmooth(t.y, 0.);

      dispersionColor *= pow(saturate(dot(nor, -rayDirection)), 1.5);

      color += saturate(dispersionColor);

      // color = pow(color, vec3(1.5));
#endif

      // color = diffuseColor;

      // Fog
      float d = max(0.0, t.x);
      color = mix(background, color, saturate(pow(clamp(fogMaxDistance - d, 0., fogMaxDistance), 2.) / fogMaxDistance));
      color *= saturate(exp(-d * 0.05));

      // color += directLighting * exp(-d * 0.0005);

      // Inner Glow
      // color += 0.5 * innerGlow(5.0 * t.w);

      // Debugging
#ifdef NO_MATERIALS
      color = diffuseColor;
#endif

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
      /* float i = saturate(t.z / (1.0 * float(maxSteps))); */
      /* vec3 glowColor = vec3(1); */
      /* const float stopPoint = 0.975; */
      /* color = mix(color, vec4(glowColor, 1.0), smoothstep(stopPoint, stopPoint + edge, i)); */

      return color;
    }
}
// Without Delta T
vec4 shade ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
  return shade ( rayOrigin, rayDirection, t, uv , 0.);
}

float myFBM (in vec2 q) {
  const float scale = 1.4;

  float f = 0.;
  f += 0.5000 * snoise2(q); q *= scale + 0.100;
  f += 0.2500 * snoise2(q); q *= scale + 0.050;
  f += 0.1250 * snoise2(q); q *= scale - 0.095;
  f += 0.0625 * snoise2(q); q *= scale + 0.150;

  return f;
}

float myFBMWarp (in vec2 q, out vec2 s, out vec2 p, out vec2 r) {
  const float scale = 1.1;

  s = vec2(
      myFBM(q + vec2(0.234, 35.234)),
      myFBM(q + vec2(34.24, 93.234)));

  p = vec2(
      myFBM(q + scale * s + vec2(0.234, 35.234)),
      myFBM(q + scale * s + vec2(34.24, 93.234)));

  r = vec2(
      myFBM(q + scale * p + vec2(0.234, 35.234)),
      myFBM(q + scale * p + vec2(34.24, 93.234)));

  return myFBM(q + scale * r);
}

float myFBMWarp (in vec2 q) {
  vec2 r = vec2(0);
  vec2 s = vec2(0);
  vec2 p = vec2(0);

  return myFBMWarp(q, s, p, r);
}

// IQ's 2D isosceles triangle
// source: https://www.shadertoy.com/view/MldcD7
float sdTriangleIsosceles( in vec2 p, in vec2 q ) {
    p.x = abs(p.x);
    vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
    float s = -sign( q.y );
    vec2 d = min( vec2( dot(a,a), s*(p.x*q.y-p.y*q.x) ),
                  vec2( dot(b,b), s*(p.y-q.y)  ));
    return -sqrt(d.x)*sign(d.y);
}

// IQ's 2D Uneven capsule
// source: https://www.shadertoy.com/view/4lcBWn
float sdUnevenCapsule( vec2 p, float r1, float r2, float h ) {
    p.x = abs(p.x);
    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(p,vec2(-b,a));
    if( k < 0.0 ) return length(p) - r1;
    if( k > a*h ) return length(p-vec2(0.0,h)) - r2;
    return dot(p, vec2(a,b) ) - r1;
}

// IQ's line sdf
// source: https://www.shadertoy.com/view/lsXGz8
float sdLine( in vec2 p, in vec2 a, in vec2 b ) {
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

float drumstick (in vec2 q, in float size, in float biteMask) {
  float n = 0.;

  q *= 1.90; // Zoom to fit

  q.y += 0.35 * size;

  // Meat
  float meatScale = 0.5;
  float meatSmallEnd = meatScale * 0.25 * size;
  float meatHeight = meatScale * 0.80 * size;
  float meat = sdUnevenCapsule(q, meatScale * 0.5 * size, meatSmallEnd, meatHeight);
  meat = smoothstep(edge, 0., meat);

  // Meat to bone crop
  float meatToBoneCrop = q.y - meatHeight + 0.03125 * size * sin(TWO_PI * 97. * q.x);
  // meatToBoneCrop = smoothstep(edge, 0., meatToBoneCrop);
  meat = min(meat, smoothstep(edge, 0., meatToBoneCrop));

  // Meat bite
  vec2 meatBiteQ = q - vec2(0.215 * size, 0.015625 * size);
  float bite = length(meatBiteQ) - (0.190 * size + 0.015625 * cnoise2(103. * meatBiteQ));
  bite = smoothstep(0.125 * edge, 0., bite);
  meat = min(meat, 1. - biteMask * bite);

  n = max(n, meat);

  // Bone
  float boneEnd = 0.85 * size;
  float boneThickness = 0.05 * size;
  float bone = sdLine(q, vec2(0, meatHeight * 0.5), vec2(0, boneEnd)) - boneThickness;
  bone = smoothstep(edge, 0., bone);
  bone = min(bone, smoothstep(edge, 0., - (meatToBoneCrop - 0.075 * size)));
  n = max(n, bone);

  // Bone nubs
  float nubBigEnd = boneThickness * 1.25;
  const float nubAngle = 0.175 * PI;
  vec2 nub1Q = q - vec2(0, boneEnd);
  nub1Q *= rotMat2(-nubAngle);
  float nub1 = sdUnevenCapsule(nub1Q, boneThickness, nubBigEnd, 0.1 * size);
  nub1 = smoothstep(edge, 0., nub1);
  n = max(n, nub1);

  vec2 nub2Q = q - vec2(0, boneEnd);
  nub2Q *= rotMat2(nubAngle);
  float nub2 = sdUnevenCapsule(nub2Q, boneThickness, nubBigEnd, 0.1 * size);
  nub2 = smoothstep(edge, 0., nub2);
  n = max(n, nub2);

  return n;
}

vec3 two_dimensional (in vec2 uv, in float generalT) {
  vec3 color = vec3(1);

  vec2 q = uv;

  // Global Timing
  float t = mod(generalT, 1.);

  // Sizing
  const float size = 0.15;
  const float side = size * 0.5;
  const float triAngleH = sqrt(0.5 * side * side);

  float n = 0.;

  // Grid space
  vec2 preModQ = q;
  vec2 c = pMod2(q, vec2(size));

  // Drumstick
  // float angle = cosT - 0.5 * length(c) + 0.0 * snoise2(c * 4.2343);
  float angle = cosT - 0.25 * dot(c, vec2(1.0, 0.5));

  q *= rotMat2(-1.5 * PI + 0.5 * PI * sin(angle));

  float showBite = step(0.6, cnoise2(24.34 * c + 0.434));
  n = max(n, drumstick(q, size, showBite));

  // Foreground/background to color
  color = vec3(n) + vec3(0.2, 0.13, 0.10);
  /* color = pow(#FCF7D5, vec3(2.2)); */
  /* vec3 lineColor = mix(#8E1EFF, #4110E8, saturate(0.6 * uv.x)); */
  /* color = mix(color, lineColor, n); */

  // DEBUG
  /* color.x += 0.25 * smoothstep(0., edge, abs(q.x) - 0.5 * size); */
  /* color.y += 0.25 * smoothstep(0., edge, abs(q.y) - 0.5 * size); */

  return color.rgb;
}

vec3 two_dimensional (in vec2 uv) {
  return two_dimensional(uv, modT);
}

vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  // return vec4(two_dimensional(uv, norT), 1);

  vec4 color = vec4(0);
  float time = norT;
  vec4 t = march(ro, rd, time);
  vec4 layer = shade(ro, rd, t, uv, time);
  return layer;
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
    // vec2 angle = normalize(vec2(1.0, 1.0));
    // gl_FragColor.rgb *= mix(
    //   vec3(1),
    //   mix(
    //     #FF1111,
    //     #00aaaa,
    //     saturate(0.25 + dot(angle, uv.xy)))
    //   , 0.15);
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
