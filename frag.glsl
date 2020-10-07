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
#define maxSteps 2048
#define maxDistance 10.0
#define fogMaxDistance 11.0

#define slowTime time * 0.2
// v3
// #define slowTime time * 0.06666667

vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);
#pragma glslify: import(./time)
const float edge = 0.0025;
const float thickness = 0.05;

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: cnoise4 = require(glsl-noise/classic/4d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
#pragma glslify: snoise4 = require(glsl-noise/simplex/4d)
//#pragma glslify: pnoise3 = require(glsl-noise/periodic/3d)
#pragma glslify: vmax = require(./hg_sdf/vmax)

#define combine(v1, v2, t, p) mix(v1, v2, t/p)

#pragma glslify: rotationMatrix = require(./rotation-matrix3)
#pragma glslify: rotationMatrix4 = require(./rotation-matrix4)

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
float sdCapsule( vec4 p, vec4 a, vec4 b, float r ) {
    vec4 pa = p - a, ba = b - a;
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

float dot2( vec3 v ) { return dot(v,v); }
float udQuad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 dc = d - c; vec3 pc = p - c;
  vec3 ad = a - d; vec3 pd = p - d;
  vec3 nor = cross( ba, ad );

  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(dc,nor),pc)) +
     sign(dot(cross(ad,nor),pd))<3.0)
     ?
     min( min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
     dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float sdBox( vec2 p, vec2 b ) {
  vec2 d = abs(p) - b;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}
float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}
float sdBox( vec4 p, vec4 b ) {
  vec4 d = abs(p) - b;
  return min(max(d.x,max(d.y,max(d.z, d.w))),0.0) + length(max(d,0.0));
}

float lpBox ( vec4 p, vec4 b ) {
  vec4 d = abs(p) - b;
  // float metricD = max(d.x,max(d.y,max(d.z, d.w)));
  // float metricD = min(d.x,min(d.y,min(d.z, d.w)));
  float metricD = dot(d, vec4(1));
  // float pr = 0.25;
  // float metricD = pow(pow(d.x, pr) + pow(d.y, pr) + pow(d.z, pr) + pow(d.w, pr), 1. / pr);
  return min(metricD, 0.0) + length(max(d,0.0));
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

// IQ's line sdf
// source: https://www.shadertoy.com/view/lsXGz8
float sdLine( in vec2 p, in vec2 a, in vec2 b ) {
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

float sdLine( in vec3 p, in vec3 a, in vec3 b ) {
    vec3 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

#define Iterations 9
#pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
#pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM, Iterations=Iterations)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
#pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=5, kifsM=kifsM, trapCalc=trapCalc)
#pragma glslify: dodecahedronFold = require(./folds/dodecahedron-fold, Iterations=1, kifsM=kifsM)

//
// #pragma glslify: fold = require(./folds)
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

// HG_SDF
float smin(float a, float b, float k){
    float f = clamp(0.5 + 0.5 * ((a - b) / k), 0., 1.);
    return (1. - f) * a + f  * b - f * (1. - f) * k;
}

float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
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
#pragma glslify: pMod4 = require(./modulo/p-mod4.glsl)
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

// Source: https://www.shadertoy.com/view/MdcXzn
const float X_REPEAT_DIST = 0.90;
const float Z_REPEAT_DIST = 1.80;
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

// Create multiple copies of an object - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
vec2 opRepLim( in vec2 p, in float s, in vec2 lim ) {
  return p-s*clamp(floor(p/s + 0.5),-lim,lim);
}

// Create multiple copies of an object - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
vec3 opRepLim( in vec3 p, in float s, in vec3 lim ) {
  return p-s*clamp(floor(p/s + 0.5),-lim,lim);
}

// source: https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
// float opRevolution( in vec3 p, in sdf2d primitive, float o )
// {
//     vec2 q = vec2( length(p.xz) - o, p.y );
//     return primitive(q)
// }

float triangleWave (in float t) {
  return 2. * abs(mod(t, 1.) - 0.5);
}

vec2 triangleWave (in vec2 t) {
  return 2. * abs(mod(t, 1.) - 0.5);
}

vec3 triangleWave (in vec3 t) {
  return 2. * abs(mod(t, 1.) - 0.5);
}

vec4 triangleWave (in vec4 t) {
  return 2. * abs(mod(t, 1.) - 0.5);
}

float lengthP(in vec2 q, in float p) {
  return pow(dot(pow(q, vec2(p)), vec2(1)), 1.0 / p);
}

float lengthP(in vec3 q, in float p) {
  return pow(dot(pow(q, vec3(p)), vec3(1)), 1.0 / p);
}

float lengthP(in vec4 q, in float p) {
  return pow(dot(pow(q, vec4(p)), vec4(1)), 1.0 / p);
}

// Inverse stereographic projection of p,
// p4 lies onto the unit 3-sphere centered at 0.
// - mla https://www.shadertoy.com/view/lsGyzm
vec4 inverseStereographic(vec3 p, out float k) {
    k = 2.0/(1.0+dot(p,p));
    return vec4(k*p,k-1.0);
}

// tdhooper's adjustment to make inverse stereographic project to be lipschitz
// continuous.
// source: https://www.shadertoy.com/view/wsfGDS
//
// Distances get warped by the stereographic projection, this applies
// some hacky adjustments which makes them lipschitz continuous.

// The numbers have been hand picked by comparing our 4D torus SDF to
// a usual 3D torus of the same size, see DEBUG.

// vec3 d
//   SDF to fix, this should be applied after the last step of
//   modelling on the torus.

// vec3 k
//   stereographic scale factor

float fixDistance(float d, float k) {
    float sn = sign(d);
    d = abs(d);
    d = d / k * 1.82;
    d += 1.;
    d = pow(d, .5);
    d -= 1.;
    d *= 5./3.;
    d *= sn;
    return d;
}

float lengthMin(in vec2 p) {
  vec2 absP = abs(p);
  return min(absP.x, absP.y);
}

float lengthMax(in vec2 p) {
  vec2 absP = abs(p);
  return max(absP.x, absP.y);
}

float lengthDot(in vec2 p) {
  vec2 absP = abs(p);
  return dot(absP, vec2(0.6));
}


float lengthCustom (in vec2 p, in float power) {
  return length(p);
}

float fTorus(vec4 p4) {

    // Torus distance
    // We want the inside and outside to look the same, so use the
    // inverted outside for the inside.
    const float pPow = 0.75;
    float d1 = lengthCustom(p4.xy, pPow) / lengthCustom(p4.zw, pPow) - 1.;
    float d2 = lengthCustom(p4.zw, pPow) / lengthCustom(p4.xy, pPow) - 1.;
    float d = d1 < 0. ? -d1 : d2;

    // Because of the projection, distances aren't lipschitz continuous,
    // so scale down the distance at the most warped point - the inside
    // edge of the torus such that it is 1:1 with the domain.
    d /= PI;

    return d;
}

const float repetitions = 12.;
vec4 pieSpace (in vec3 p, in float relativeC) {
  float angle = relativeC * TWO_PI / repetitions;
  p.xz *= rotMat2(angle);
  float c = pModPolar(p.xz, repetitions);
  p.xz *= rotMat2(-angle);
  return vec4(p, c);
}

float r = 1.50;
float sdHollowBox (in vec3 q, in vec3 r, in float thickness) {
  float b = sdBox(q, r);

  // crop inners
  vec3 cropR = r - thickness;
  vec3 cropQ = q;
  float crop = sdBox(cropQ, vec3(30, cropR.y, cropR.z));

  if (abs(q.y) > abs(q.x)) {
    cropQ.yx = cropQ.xy;
  }
  crop = min(crop, sdBox(cropQ, vec3(30, cropR.x, cropR.z)));

  if (abs(q.z) > abs(q.x)) {
    cropQ.zx = cropQ.xz;
  }
  crop = min(crop, sdBox(cropQ, vec3(30, cropR.x, cropR.y)));

  return max(b, -crop);
}

float gridBump ( in vec3 q, float size ) {
  vec3 c = pMod3(q, vec3(size));

  vec3 absQ = abs(q);
  return vmax(absQ) / size;
}

vec3 sphericalCoords (in vec3 q) {
  float a = atan(q.z, q.x);
  float arc = atan(q.y, length(q.xz));
  return vec3(
      a,
      arc,
      length(q));
}

vec3 map (in vec3 p, in float dT, in float universe) {
  vec3 d = vec3(maxDistance, 0, 0);
  float minD = 0.;

  // p *= globalLRot;

  vec3 q = p;
  vec4 z = vec4(q, q.x);

  const float size = 0.1;
  float t = mod(dT, 1.);

  float warpScale = 1.5;

  // Warp
  // vec3 wQ = q;
  vec4 wQ = z;

  wQ += warpScale * 0.100000 * cos( 5.638 * wQ.yzwx + cosT );
  wQ.xyz = twist(wQ.xzy, 1.00 * wQ.z);
  wQ += warpScale * 0.050000 * cos(11.237 * wQ.yzxw + cosT );
  wQ += warpScale * 0.025000 * cos(23.123 * wQ.yzxw + cosT );
  // wQ.xzy = twist(wQ.xyz, 2.00 * wQ.y);
  wQ.xzy = twist(wQ.xyz, 3.00 * wQ.y);
  wQ += warpScale * 0.012500 * cos(37.123 * wQ.yzxw + cosT );

  // q = wQ.xyz;
  z = wQ;

  mPos = z.xyz;
  vec3 o = vec3(lpBox(z, vec4(r)), 0, 0);
  // vec3 o = vec3(sdBox(z, vec4(r)), 0, 0);
  // o.x += 0.001 * cellular(3. * z.xyz);
  d = dMin(d, o);

  d.x *= 0.0125;

  return d;
}

vec3 map (in vec3 p, in float dT) {
  return map(p, dT, 0.);
}

vec3 map (in vec3 p) {
  return map(p, 0., 0.);
}

vec4 march (in vec3 rayOrigin, in vec3 rayDirection, in float deltaT, in float universe) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  const float deltaTDelta = 0.000;
// #define OVERSTEP 1
#ifdef OVERSTEP
  // Source: https://www.shadertoy.com/view/MdcXzn
  const int halfMax = (maxSteps / 2);
  for( int i = 0; i < halfMax; i++ ) {
       vec3 d = map(rayOrigin + rayDirection * t, deltaT, universe);

       if( d.x < epsilon * 100.0 ) break;
       t += d.x;
       if (t > maxDistance) break;
      deltaT += deltaTDelta;
   }

   t -= Z_REPEAT_DIST * 0.5;

   for( int i = 0; i < halfMax; i++ ) {
       vec3 d = map(rayOrigin + rayDirection * t, deltaT, universe);

       if( d.x<epsilon ) return vec4(t + d.x, d.y, float(i), d.z);

       t += min(d.x, Z_REPEAT_DIST * 0.2);
      deltaT += deltaTDelta;
   }
#else
  for (int i = 0; i < maxSteps; i++) {
    vec3 d = map(rayOrigin + rayDirection * t, deltaT, universe);
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

vec4 march (in vec3 rayOrigin, in vec3 rayDirection, in float deltaT) {
  return march (rayOrigin, rayDirection, deltaT, 0.);
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

// --- Patterns ---
// #pragma glslify: checker = require(glsl-checker)
// #pragma glslify: herringBone = require(./patterns/herring-bone)

#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: rgb2hsv = require(./rgb2hsv.glsl)
#pragma glslify: hsb2rgb = require(./color-map/hsb2rgb)

float n1 = 1.;
float n2 = 2.4;
const float amount = 0.25;

float gM = 0.;
vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float dNR = dot(-rd, gNor);

  float spread = 1.; // saturate(1.0 - 1.0 * pow(dNR, 9.));
  // // float n = smoothstep(0., 1.0, sin(150.0 * rd.x + 0.01 * noise(433.0 * rd)));

  // float startPoint = 0.0;

  // vec3 spaceScaling = 0.4 * vec3(0.734, 1.14, 0.2);
  // float n = ncnoise3(spaceScaling * rd + startPoint);
  // n = smoothstep(0.0, 0.80, n);

  // vec3 spaceScaling = vec3(0.8);
  // float n = vfbmWarp(spaceScaling * rd + startPoint);
  // n = smoothstep(0.525, 0.80, n);

  // vec3 spaceScaling = vec3(9.8);
  // float n = vfbm4(spaceScaling * rd + startPoint);
  // n = smoothstep(0.125, 0.85, n);

  // float n = smoothstep(0.9, 1.0, sin(TWO_PI * (dot(vec2(8), rd.xz) + 2.0 * cnoise3(1.5 * rd)) + time));

  // float n = cnoise3(8.5 * rd);
  // n = smoothstep(-0.1, 0.9, n);

  // float n = 0.6 + 0.4 * sin(dot(vec3(PI), sin(3.18 * rd + sin(1.38465 * rd.yzx))));

  vec3 dI = vec3(dNR);
  dI += 0.2 * snoise3(0.1 * rd);
  dI += 0.3 * pow(dNR, 3.);

  dI *= angle3C;
  dI += offset.x;

  // color = vec3(0.098039, 0.960784, 0.960784) + vec3(0.2, 0.4, 0.2) * cos( TWO_PI * (vec3(2, 1, 1) * dI + vec3(0, 0.25, 0.25)) );
  color = 0.5 + 0.5 * cos( TWO_PI * ( dI + vec3(0, 0.33, 0.67) ) );

  color *= spread;

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

vec4 shade ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv, in float generalT );

float phaseHerringBone (in float c) {
  return c * 0.124 + 2. * (1. + mod(c, 2.)) * cosT;
}

#pragma glslify: herringBone = require(./patterns/herring-bone, phase=phaseHerringBone)

vec3 baseColor (in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap, in float t) {
  vec3 color = vec3(1.25);

  float dNR = dot(nor, -rd);
  vec3 dI = vec3(dNR);

  dI += 0.2 * pos;
  dI += 0.3 * pow(dNR, 3.);

  dI *= angle1C;
  dI += angle2C;

  color = 0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0,-0.33, 0.67)));
  dI += 0.1 * nor;
  color.rg += 0.5 + 0.5 * cos(TWO_PI * (dI.xy + vec2(0, 0.33)));

  color *= 0.65;

  // color = vec3(
  //     color.r + color.b,
  //     color.r + color.g,
  //     color.g + color.b);
  // color *= 0.435;

  gM = m;

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

    // float dNR = dot()
    // lights[0] = light(normalize(vec3(  0.15, 0.25, 1.0)), #FFFFFF, 1.0);
    lights[0] = light(vec3(0.0, 0.81, 1.0), #FFFFFF, 1.0);
    lights[1] = light(vec3(0.5, 1.2,  1.0), #FFFFFF, 1.0);
    lights[2] = light(vec3(-0.1, 1.0, -1.0), #FFFFFF, 1.0);
    // lights[3] = light(vec3( 0.3, 0.8, -0.4), #FFFFFF, 1.0);
    // lights[4] = light(vec3(-0.4, -.2, -1.0), #FFFFFF, 1.0);

    const float universe = 0.;
    background = getBackground(uv, universe);

    float backgroundMask = 1.;
    // Allow anything in top right corner
    // backgroundMask = max(backgroundMask, smoothstep(0., edge, dot(uv, vec2(1)) + 0.05));

    if (t.x>0. && backgroundMask > 0.) {
      vec3 color = vec3(0.0);

      // Normals
      vec3 nor = getNormal2(pos, 0.005 * t.x, generalT);
      // float bumpsScale = 5.75;
      // float bumpIntensity = 0.25;
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

      // Material Types

      float occ = calcAO(pos, nor);
      float amb = saturate(0.5 + 0.5 * nor.y);
      float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 1.0;
      float specCo = 0.5;

      float specAll = 0.0;

      vec3 directLighting = vec3(0);
      for (int i = 0; i < NUM_OF_LIGHTS; i++) {
        vec3 lightPos = lights[i].position;
        float diffMin = 0.7;
        float dif = max(diffMin, diffuse(nor, normalize(lightPos)));
        float spec = pow(clamp( dot(ref, normalize(lightPos)), 0., 1. ), 64.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        float shadowMin = 0.5;
        float sha = max(shadowMin, softshadow(pos, normalize(lightPos), 0.001, 4.75));
        dif *= sha;

        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        specAll += specCo * spec;
        lin += pow(specCo * spec, 4.);

        // Ambient
        lin += 0.0 * amb * diffuseColor;
        // dif += 0.000 * amb;

        float distIntensity = 1.; // lights[i].intensity / pow(length(lightPos - gPos), 1.0);
        distIntensity = saturate(distIntensity);
        color +=
          saturate((dif * distIntensity) * lights[i].color * diffuseColor)
          + saturate(distIntensity * mix(lights[i].color, vec3(1), 0.0) * lin * mix(diffuseColor, vec3(1), 1.0));

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
      reflectColor += 0.05 * reflection(pos, reflectionRd);
      color += reflectColor;

      // vec3 refractColor = vec3(0);
      // vec3 refractionRd = refract(rayDirection, nor, 1.5);
      // refractColor += 0.15 * textures(refractionRd);
      // color += refractColor;

#ifndef NO_MATERIALS

      // vec3 dispersionColor = dispersionStep1(nor, normalize(rayDirection), n2, n1);
      vec3 dispersionColor = dispersion(nor, rayDirection, n2, n1);

      float dispersionI = 1.; // dot(nor, -rayDirection);
      dispersionColor *= dispersionI;

      color += saturate(dispersionColor);
      // color = saturate(dispersionColor);

#endif

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

      // // Post processing coloring
      // float dNR = dot(gNor, -rayDirection);

      // vec3 midColor = #E68CA1 * mix(dNR, 1., 0.9);
      // const vec3 highlightColor = saturate(#E68CA1 * 1.5);
      // const vec3 shadowColor = #8C4858;

      // vec3 shaded = color;
      // color = mix(midColor, highlightColor, step(0.55, shaded.x));
      // color = mix(color, shadowColor, step(0.7, 1. - shaded.x));

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
      float i = saturate(t.z / (1.0 * float(maxSteps)));
      vec3 glowColor = vec3(0);
      const float stopPoint = 0.5;
      color = mix(color, vec4(glowColor, 1.0), smoothstep(stopPoint, stopPoint + edge, i));

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

float chopSquareChisel (vec2 q, float stage1, float stage2, float stage3, float stage4) {
  float chop = 1.;

  chop *= 1. - step(1., stage1) * step(0., -q.x) * step(0., -q.y);
  chop *= 1. - step(1., stage2) * step(0., -q.x) * step(0.,  q.y);
  chop *= 1. - step(1., stage3) * step(0.,  q.x) * step(0.,  q.y);

  return chop;
}

float chopSquareSub (vec2 q, float stage1, float stage2, float stage3, float stage4) {
  float chop = 1.;

  chop *= 1. - step(1., stage1) * step(0.,  q.x) * step(0.,  q.y);
  chop *= 1. - step(1., stage2) * step(0., -q.x) * step(0.,  q.y);
  chop *= 1. - step(1., stage3) * step(0., -q.x) * step(0., -q.y);

  return chop;
}

vec3 getRingColor (in float angle) {
  return 0.5 + 0.5 * cos(TWO_PI * (vec3(1) * angle + vec3(-0.5, 0.3, 0.6)));
}

vec3 getRing (in float dir, in vec2 q, in vec3 color, in float ringWidth, in float offset) {
  // Color
  vec2 pol = vec2(q.x, q.y);
  pol.x /= PI;
  pol.x += 1.;
  pol.x *= 0.5;
  pol.x += dir * 15. * (0.5 + 0.5 * sin(0.5 * cosT - PI * 0.5));
  pol.x += dir * PI * 0.5;
  pol.x += -0.123 * offset;

  // Offset for edges
  float edgeStart = 0.9 * 0.5 * ringWidth;
  float isEdge = smoothstep(edgeStart, edgeStart + edge, abs(pol.y));

  // Offset by a third rotation
  pol.x += 0.218 * isEdge * (1. - 2. * step(0., pol.y));

  vec3 ringColor = getRingColor(pol.x);

  // Ring Mask
  float isRing = smoothstep(edge, 0., abs(pol.y) - 0.5 * ringWidth);
  return mix(color, ringColor, isRing);
}

float triangleMotion(in vec2 q, in float r, in float t) {
  float n = 0.;

  const float bigR = 0.3;

  vec2 offsetPos = vec2(bigR, 0);
  offsetPos = mix(offsetPos, bigR * vec2(cos(2.0 * PI / 3.0), sin(2.0 * PI / 3.0)), saturate(t*3.0));
  offsetPos = mix(offsetPos, bigR * vec2(cos(2.0 * PI * 2. / 3.0), sin(2.0 * PI * 2. / 3.0)), saturate(t*3.0 - 1.0));
  offsetPos = mix(offsetPos, bigR * vec2(1., 0.), saturate(t*3.0 - 2.0));

  n = length(q - offsetPos) - r;
  n = smoothstep(edge, 0., n);

  return n;
}

float swirlMask (in vec2 q, in float nOffset) {
  float a = atan(q.y, q.x);
  vec2 pol = vec2(a / TWO_PI + 0.5, length(q));

  vec2 scaleN = 0.8 * scale * vec2(49, 10);
  float n1 = fbmWarp(scaleN * pol + nOffset);
  float n2 = fbmWarp(scaleN * vec2(1. - pol.x, pol.y) + nOffset);
  return mix(n1, n2, saturate((pol.x - 0.8) / 0.2));
}

vec2 eclipseQ (in vec2 abR, in float t) {
  return abR * vec2(cos(t), sin(t));
}

// IQ's ellipse 2D sdf
// source: https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float sdEllipse( in vec2 p, in vec2 ab )
{
    p = abs(p); if( p.x > p.y ) {p=p.yx;ab=ab.yx;}
    float l = ab.y*ab.y - ab.x*ab.x;
    float m = ab.x*p.x/l;      float m2 = m*m; 
    float n = ab.y*p.y/l;      float n2 = n*n; 
    float c = (m2+n2-1.0)/3.0; float c3 = c*c*c;
    float q = c3 + m2*n2*2.0;
    float d = c3 + m2*n2;
    float g = m + m*n2;
    float co;
    if( d<0.0 )
    {
        float h = acos(q/c3)/3.0;
        float s = cos(h);
        float t = sin(h)*sqrt(3.0);
        float rx = sqrt( -c*(s + t + 2.0) + m2 );
        float ry = sqrt( -c*(s - t + 2.0) + m2 );
        co = (ry+sign(l)*rx+abs(g)/(rx*ry)- m)/2.0;
    }
    else
    {
        float h = 2.0*m*n*sqrt( d );
        float s = sign(q+h)*pow(abs(q+h), 1.0/3.0);
        float u = sign(q-h)*pow(abs(q-h), 1.0/3.0);
        float rx = -s - u - c*4.0 + 2.0*m2;
        float ry = (s - u)*sqrt(3.0);
        float rm = sqrt( rx*rx + ry*ry );
        co = (ry/sqrt(rm-rx)+2.0*g/rm-m)/2.0;
    }
    vec2 r = ab * vec2(co, sqrt(1.0-co*co));
    return length(r-p) * sign(p.y-r.y);
}

float starSpeed (in float i, in float numStars) {
  // return 1.;
  // return 1. / numStars;
  return 4. / numStars;
}

float starPhase (in float i, in float numStars) {
  return TWO_PI / numStars; // / (i + 2.);
}

float numStarsCalc (in float i) {
  return (i + 1.) * 8.;
}

float starEllipseStar (in vec2 q, in vec2 abR, in float starR, in float t) {
  float n = 1.;

  vec2 eclQ = eclipseQ(abR, t);
  vec2 relQ = q - eclQ;
  float angle = atan(relQ.y, relQ.x) + PI;
  angle /= TWO_PI; // normalize
  float normT = mod(t + PI, TWO_PI) / TWO_PI;

  // Ring
  float ring = sdEllipse(q, abR);
  ring = abs(ring);
  ring = smoothstep(0.50 * edge, 0., ring);
  // Make ring a tail
  ring *= pow(1. - saturate(mod(normT - angle, 1.) / 0.5), 0.5);

  // ring *= mod(t, TWO_PI) / TWO_PI;
  n *= 1. - saturate(10. * d) * ring;

  // Star
  float star = length(relQ) - starR;
  star = smoothstep(edge, 0., star);
  n *= 1. - star;

  return n;
}

float starEllipse (in vec2 q, in vec2 abR, in float starR, in float t, in float i) {
  float n = 1.;

  // Ring
  float ring = sdEllipse(q, abR);
  ring = abs(ring);
  ring = smoothstep(0.50 * edge, 0., ring);
  // n *= 1. - saturate(d) * ring;

  vec2 eclQ = eclipseQ(abR, -t);

  // i += 1.;
  float speed = 1.;
  float phase = 0.;
  float noiseTune = 1000. * angle3C; // 912.423;
  float noiseTuneAmp = PI * angle2C;
  float numStars = 0.;
  float rotAmount = PI * -0.036;

  abR = angle1C * abR;
#define numStarsInOrbit 5.
  numStars = numStarsInOrbit; // Fix because I cant use dynamic variables in loop conditionals // numStarsCalc(i);
  speed = starSpeed(i, numStars);
  phase = starPhase(i, numStars);
  for (float j = 0.; j < numStarsInOrbit; j++) {
    n *= starEllipseStar(q - eclQ, abR, starR, speed * t + j * phase + noiseTuneAmp * noise(vec2((i + 1.) * noiseTune )));
  }

  return n;
}

float starEllipse (in vec2 q, in vec2 abR, in float starR, in float t) {
  return starEllipse(q, abR, starR, t, 1.);
}

vec3 gradient (in float i) {
  const float increment = 0.3333;

  float lowerEdge = 0.;

  vec3 color = #9CF0C0;
  color = mix(color, #D3A8F0, smoothstep(lowerEdge, lowerEdge + increment, i));
  lowerEdge += increment;
  color = mix(color, #9FF08F, smoothstep(lowerEdge, lowerEdge + increment, i));
  lowerEdge += increment;
  color = mix(color, #F08178, smoothstep(lowerEdge, lowerEdge + increment, i));
  lowerEdge += increment;
  color = mix(color, #E1F084, smoothstep(lowerEdge, lowerEdge + increment, i));
  lowerEdge += increment;

  return color;
}

float halfGridSplit (in vec2 q, in float scale, in float t) {
  const float size = 0.2;
  // Initial R is fudged so the slightly bigger r of the two spheres overlapping
  // w/ a smooth min will equal the radius of one sphere by the end.
  float initialR = 0.15 * size * 0.65;
  float r = mix(initialR, 0.15 * size * scale, t);

  float rowOdd = mod(floor((q.y + 0.5 * size) / size), 2.);
  q.x += 1.0 * size * rowOdd;

  vec2 c = pMod2(q, vec2(2. * size, size));

  float offset = 0.075 * length(c);
  float separation = size * 0.5 * smoothstep(offset, offset + 0.5, t);

  float s1 = length(q - vec2(separation, 0)) - r;
  float s2 = length(q + vec2(separation, 0)) - r;
  float n = fOpUnionRound(s1, s2, 2. * r);
  n = smoothstep(edge, 0., n);

  return n;
}

vec2 cellQ;
float cellImage (in vec2 q, in float size, in vec2 offset) {
  offset *= 0.1;

  const float warpScale = 0.40;

  float r = size * 0.023438;

  offset *= rotMat2(0.25 * PI);

  float a = atan(offset.y, offset.x);
  vec2 shift = vec2(0);
  shift += warpScale * 0.1000 * cos( 2. * (offset.yx + shift.yx) + cosT );
  shift *= warpScale * 2.5 * sin(cosT - 3. * length(offset + shift));
  shift += warpScale * 0.0500 * cos( 3. * (offset.yx + shift.yx) + cosT );
  shift += warpScale * 0.0250 * cos( 4. * (offset.yx + shift.yx) + cosT );
  // shift 

  q += shift;
  cellQ = q;
  float s = length(q) - r;
  return s;
}

vec3 layerColor (in vec2 q, in float phase) {
  return 0.5 + 0.5 * cos(TWO_PI * (0.123 * vec3(q, 0.) + phase + vec3(0, 0.33, 0.67)));
}

// vec3 pieMap (in vec3 q, in float c) {
//   vec3 d = vec3(maxDistance, 0, 0);
//   vec2 uv = q.xy;

//   float angleSegment = TWO_PI * 0.142857;
//   float numberOfPositionsShifted = 2.;
//   uv *= rotMat2(-numberOfPositionsShifted * angleSegment * norT);

//   const float groundR = 0.25;

//   float a = atan(uv.y, uv.x);
//   float h = uv.x - (groundR + bounceOffset(c, angleSegment, numberOfPositionsShifted, a));

//   vec3 t = vec3(abs(h) - 0.0125, 1, -a);
//   // d = dMin(d, t);
//   vec2 ballUV = uv;
//   float ballR = 0.05;

//   ballUV.x -= groundR + bounceOffset(c, angleSegment, numberOfPositionsShifted);
//   vec3 b = vec3(length(ballUV) - ballR, 0, 0);
//   d = dMin(d, b);

//   float groundThickness = 0.05;
//   vec3 g = vec3(abs(length(uv) - (groundR - 0. * ballR - groundThickness)), 0, 0);
//   d = dMin(d, g);

//   return d;
// }

// #pragma glslify: ringSpace = require(./pie-slice, map=pieMap, repetitions=7, cosT=cosT, )

const float size = 0.035;
float map (in vec2 q, in vec2 c) {
  float oddColumn = mod(c.x, 2.);
  float oddRow = mod(c.y, 2.);
  float evenColumn = 1. - mod(c.x + oddRow, 2.);
  float evenRow = 1. - mod(c.y + oddColumn + evenColumn, 2.);

  q.x += 0.5 * size * oddRow;
  float t = mod(norT, 1.);

  // float phase = dot(c, vec2(1));
  // float phase = -length(c);
  // Square metric
  // vec2 absC = abs(c);
  // float phase = -max(absC.x, absC.y);
  // float r = mix(0.95, 1.249, (0.5 + 0.5 * sin(cosT + phase))) * size;
  // float r = angle1C * size;

  const float warpScale = 0.;

  vec2 absC = abs(c);
  // float warpPhase = dot(c, vec2(0.10));
  // float warpPhase = -0.123 * length(c);
  float warpPhase = -2.413 * max(absC.x, absC.y);
  q += warpScale * 0.1000 * cos(19. * q.yx + cosT + warpPhase);
  q += warpScale * 0.0500 * cos(61. * q.yx + cosT + warpPhase);
  q += warpScale * 0.0250 * cos(65. * q.yx + cosT + warpPhase);

  // Random positions
  // q += 0.7 * size * vec2(
  //     snoise2(0.2 * sin(cosT - 2.9 * length(c)) + c * 0.123),
  //     snoise2(0.2 * sin(cosT - 1.9 * length(c)) + c * 8.123 + 0.777778));

  // Circle
  float d = length(q);
  d = smoothstep(0., edge, d - size * 0.25);

  // // Hexagon
  // float d = sdHexPrism(vec3(q, 0), vec2(0., 1));

  // // Square
  // vec2 absQ = abs(q);
  // float d = max(absQ.x, absQ.y);

  return d;
}

vec2 getCoords ( in float r, in float i, in float num ) {
  float angleIncrement = TWO_PI / num;

  return r * vec2(
      cos(i * angleIncrement),
      sin(i * angleIncrement));
}

float lineTrap (in vec2 point) {
  vec2 pointOnLine = 2. * (colors1.xy - vec2(0.5));
  vec2 lineNor = normalize(vec2(0.4, 0.8));

  vec2 polToP = pointOnLine - point;
  return length(polToP - dot(polToP, lineNor) * lineNor);
}

#pragma glslify: neighborGrid = require(./modulo/neighbor-grid, map=map, maxDistance=maxDistance, numberOfNeighbors=2.)
vec3 two_dimensional (in vec2 uv, in float generalT) {
  vec3 color = vec3(0);

  vec2 q = uv;
  // q.yx = q.xy;

  q *= scale;
  q += offset.xy;
  float warpScale = 1.;

  // Global Timing
  float t = mod(generalT + 0.0, 1.);

  // Fractal Warp
  float minD = maxDistance;
  float avgD = 0.;
  const int iterations = 50;
  float iteration = 0.;
  for (int i = 0; i < iterations; i++) {
    float fI = float(i);

    // q.x = abs(q.x);
    // q *= angle2C;;
    // // q.xy = abs(q.xy);
    // q.x += angle3C;
    // q *= rotMat2(offset.x * PI);

    // Julia set
    vec2 c = vec2(angle1C, angle2C); // vec2(0.4, 0.6);
    c += 0.01 * cos(cosT + vec2(0, 0.5 * PI) + 0. * uv);
    // q^2 = q.x^2 + 2. * q.x * q.y * i + - q.y^2
    q = vec2(
        q.x * q.x - q.y * q.y,
        2. * q.x * q.y);
    q += c;

    float dis = dot(q, q);
    if (dis > 8.) break;

    // float d = length(q);
    float d = dot(q, q);
    // float pr = 5.5;
    // float d = pow(dot(pow(q, vec2(pr)), vec2(1)), 1. / pr);

    // float d = lineTrap(q);

    avgD += d;
    minD = min(minD, d);
    iteration += 1.;
  }

  avgD /= float(iterations);

  // Final shape
  float n = length(q);
  // n = smoothstep(edge, 0., n);

  // n *= 0.4;

  // Option 1
  float option1I = angle3C * minD;
  // float option1I = iteration / float(iterations);
  // vec3 option1 = 0.5 + 0.5 * cos( TWO_PI * (option1I + vec3(0, 0.33, 0.66)) );
  vec3 option1 = vec3(0.5 + 0.5 * sin(TWO_PI * option1I));

  // Option 2
  // float option2I = 0.5 + 0.5 * sin(offset.y * TWO_PI * n);
  // float option2I = iteration / float(iterations);
  // float option2I = n;
  float option2I = offset.z * minD;
  // vec3 option2 = vec3(smoothstep(0., edge, sin(offset.y * TWO_PI * n)));
  vec3 option2 = 0.5 + 0.5 * cos( TWO_PI * (option2I + vec3(0, 0.33, 0.67)) );
  // vec3 option2 = option2I * mix(#026FF2, vec3(1), pow(option2I, 3.0));

  // color = mix(option2, option1, norT);
  color = option1;
  // color = option2;

  float mask = smoothstep(0., edge, n - r);
  // color = mix(color, background, mask);

  if (iteration != float(iterations)) {
    color = option2;
  }

  return color.rgb;
}

vec3 two_dimensional (in vec2 uv) {
  return two_dimensional(uv, norT);
}

vec3 screenBlend (in vec3 a, in vec3 b) {
  return 1. - (1. - a) * (1. - b);
}

vec3 overlay (in vec3 a, in vec3 b) {
  vec3 colorOut = vec3(0);
  vec3 screen = screenBlend(a, b);

  colorOut.x = mix(2. * a.x * b.x, screen.x, step(0.5, a.x));
  colorOut.y = mix(2. * a.y * b.y, screen.y, step(0.5, a.y));
  colorOut.z = mix(2. * a.z * b.z, screen.z, step(0.5, a.z));

  return colorOut;
}

vec3 softLight1 (in vec3 a, in vec3 b) {
  return (1. - 2. * b) * a * a + 2. * b * a;
}

vec3 softLight2 (in vec3 a, in vec3 b) {
  return pow(a, pow(vec3(2.), 2. * (0.5 - b)));
}

vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  return vec4(two_dimensional(uv, norT), 1);

  // vec3 color = vec3(0);

  // const int slices = 5;
  // for (int i = 0; i < slices; i++) {
  //   float fI = float(i);
  //   // vec3 layerColor = 0.5 + 0.5 * cos(TWO_PI * (fI / float(slices) + vec3(0, 0.33, 0.67)));
  //   vec3 layerColor = vec3(
  //       saturate(mod(fI + 0., 3.)),
  //       saturate(mod(fI + 1., 3.)),
  //       saturate(mod(fI + 2., 3.))
  //   );

  //   vec3 dI = vec3(fI);
  //   dI += 0.4 * uv.x;
  //   dI += 0.3 * dot(uv, vec2(1));
  //   // dI += 0.2 * snoise2(3. * uv);
  //   layerColor += 1.0 * (0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.33, 0.67))));
  //   layerColor *= 0.6;

  //   // Add black layer as first layer
  //   // layerColor *= step(0.5, fI);

  //   // layerColor = pow(layerColor, vec3(4 + slices));

  //   float mask = two_dimensional(uv, cosT + 1.148 * fI / float(slices)).x;
  //   layerColor *= mask;
  //   // if (i == 0) {
  //   //   color = layerColor;
  //   // } else {
  //   //   color = overlay(color, layerColor);
  //   // }
  //   // color *= layerColor;

  //   // vec3 layerColorA = softLight2(color, layerColor);
  //   // vec3 layerColorB = color * layerColor;
  //   // layerColor = layerColorA + layerColorB;
  //   // layerColor *= 0.85;

  //   // layerColor = overlay(color, layerColor);
  //   // layerColor = screenBlend(color, layerColor);
  //   // color = mix(color, layerColor, 0.3);

  //   // Add
  //   color += layerColor;

  //   // Multiply
  //   // color *= layerColor;

  //   // Pseudo Multiply
  //   // color = mix(color, color * layerColor, mask);
  // }

  // color = pow(color, vec3(1.75));
  // color /= float(slices);

  // // float posX = mod(18. * abs(uv.x), 3.);
  // // color = vec3(
  // //     saturate(mod(posX + 0., 3.)),
  // //     saturate(mod(posX + 1., 3.)),
  // //     saturate(mod(posX + 2., 3.))
  // //     );

  // return vec4(color, 1.);

  float time = 2. * norT;
  vec4 t = march(ro, rd, time);
  vec4 layer = shade(ro, rd, t, uv, time);
  return layer;
}

void main() {
  // Update timing
  modT = mod(time, totalT);
  norT = modT / totalT;
  cosT = TWO_PI / totalT * modT;

    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;

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
