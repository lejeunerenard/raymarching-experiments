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
// #define DOF 1

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
uniform sampler2D sdf2DTexture;

uniform sampler2D uninitTex;

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
#define maxSteps 512
#define maxDistance 10.0
#define fogMaxDistance 5.0

#define slowTime time * 0.2
// v3
// #define slowTime time * 0.06666667

vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);
vec3 dRd = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);
#pragma glslify: import(./time)
const float edge = 0.0025;
const float thickness = 0.01;

// Dispersion parameters
float n1 = 1.;
float n2 = 1.5;
const float amount = 0.05;

// Dof
float doFDistance = length(cameraRo) - 0.275;

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: vmax = require(./hg_sdf/vmax)

float vmin (in vec2 t) {
  return min(t.x, t.y);
}

float vmin (in vec3 t) {
  return min(t.x, min(t.y, t.z));
}

#define combine(v1, v2, t, p) mix(v1, v2, t/p)

float range (in float start, in float stop, in float t) {
  return saturate((t - start) / (stop - start));
}

// Math
#pragma glslify: rotationMatrix = require(./rotation-matrix3)
#pragma glslify: rotationMatrix4 = require(./rotation-matrix4)

vec4 qSquare (in vec4 q) {
  return vec4(q.x*q.x - q.y*q.y - q.z*q.z - q.w*q.w, 2.0*q.x*q.yzw);
}

float qLength2 (in vec4 q) {
  return dot(q, q);
}

vec4 qCube ( in vec4 q ) {
  vec4  q2 = q*q;
  return vec4(q.x  *(    q2.x - 3.0*q2.y - 3.0*q2.z - 3.0*q2.w), 
      q.yzw*(3.0*q2.x -     q2.y -     q2.z -     q2.w));
}

// return: [0, 1]
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

// Noise
#pragma glslify: cnoise4 = require(glsl-noise/classic/4d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
#pragma glslify: snoise4 = require(glsl-noise/simplex/4d)

//#pragma glslify: pnoise3 = require(glsl-noise/periodic/3d)
float ncnoise2(in vec2 x) {
  return smoothstep(-1.00, 1.00, cnoise2(x));
}
float ncnoise3(in vec3 x) {
  return smoothstep(-1.00, 1.00, cnoise3(x));
}

// 3D noise function (IQ)
// float noise(vec3 p) {
//   vec3 ip=floor(p);
//     p-=ip;
//     vec3 s=vec3(7,157,113);
//     vec4 h=vec4(0.,s.yz,s.y+s.z)+dot(ip,s);
//     p=p*p*(3.-2.*p);
//     h=mix(fract(sin(h)*43758.5),fract(sin(h+s.x)*43758.5),p.x);
//     h.xy=mix(h.xz,h.yw,p.y);
//     return mix(h.x,h.y,p.z);
// }
// source: https://www.shadertoy.com/view/lsl3RH
float noise( in vec2 x ) {
  return sin(1.5*x.x)*sin(1.5*x.y);
}
float noise( in vec3 x ) {
  return sin(1.5*x.x)*sin(1.5*x.y)*sin(1.5*x.z);
}

// Source: https://www.shadertoy.com/view/fsyGD3
float h21 (vec2 a) {
  return fract(sin(dot(a.xy,vec2(12.9898,78.233)))*43758.5453123);
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

// IQ's 2D Line Segment SDF
// Source: http://iquilezles.untergrund.net/www/articles/distfunctions2d/distfunctions2d.htm
float sdSegment( in vec2 p, in vec2 a, in vec2 b ) {
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
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
float sdTorus28( vec3 p, vec2 t )
{
  vec2 q = vec2(length8(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

// Maximetric in radius direction & euclidean round the radius loci
float sdTorus2M( vec3 p, vec3 t )
{
  vec2 d = abs(p.xz) - t.xy;
  vec2 q = vec2(min(vmax(d), 0.) + length(max(d, 0.)),p.y);
  return length(q)-t.z;
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

// IQ's arc SDF
float sdArc( in vec2 p, in vec2 sca, in vec2 scb, in float ra, float rb )
{
    p *= mat2(sca.x,sca.y,-sca.y,sca.x);
    p.x = abs(p.x);
    float k = (scb.y*p.x>scb.x*p.y) ? dot(p,scb) : length(p);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
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
#pragma glslify: halfTetraFold = require(./folds/half-tetrahedral)
#pragma glslify: tetraFold = require(./folds/tetrahedral)
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
vec2 dSMax (vec2 d1, vec2 d2, in float r) {
  float h = saturate(0.5 + 0.5 * (d1.x - d2.x) / r);
  float d = mix(d2.x, d1.x, h) + h * ( 1.0 - h ) * r;
  return vec2(d, (d1.x < d2.x) ? d1.y : d2.y);
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
#pragma glslify: gyroid = require(./model/gyroid)

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
#pragma glslify: sineIn = require(glsl-easings/sine-in)
#pragma glslify: quart = require(glsl-easings/quadratic-in-out)
#pragma glslify: quartIn = require(glsl-easings/quadratic-in)
#pragma glslify: quartOut = require(glsl-easings/quadratic-out)
#pragma glslify: quint = require(glsl-easings/quintic-in-out)
#pragma glslify: quintIn = require(glsl-easings/quintic-in)
#pragma glslify: quintOut = require(glsl-easings/quintic-out)
// #pragma glslify: elasticInOut = require(glsl-easings/elastic-in-out)
#pragma glslify: elasticOut = require(glsl-easings/elastic-out)
// #pragma glslify: elasticIn = require(glsl-easings/elastic-in)


// vec3 versions
vec3 expo (in vec3 x) {
  return vec3(
      expo(x.x),
      expo(x.y),
      expo(x.z)
      );
}
vec3 quad (in vec3 x) {
  return vec3(
      quad(x.x),
      quad(x.y),
      quad(x.z)
      );
}

vec3 expoWave (in vec3 q) {
  q = triangleWave(q);
  return -1. + 2. * expo(q);
}

#pragma glslify: voronoi = require(./voronoi, edge=edge, thickness=thickness, mask=sqrMask)
#pragma glslify: sdFBM = require(./model/sdf-fbm, REPS=7, smoothness=0.6, clipFactor=0.1, smax=smax, smin=smin);
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
// TODO figure out if this outputs from [0, 1]
float sigmoid ( in float x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

vec3 sigmoid ( in vec3 x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

#pragma glslify: gyroidTriangle = require(./model/gyroid-trianglewave, triangleWave=triangleWave)
#pragma glslify: gyroidExpo = require(./model/gyroid-trianglewave, triangleWave=expoWave)

// Smooth polar mod by Paulo Falcao
// source: https://www.shadertoy.com/view/NdS3Dh
vec2 smoothPModPolar(in vec2 q, in float repetitions, in float smoothness, in float correction, in float displacement) {
  repetitions *= 0.5;
  float k = length(q);
  float x = asin(sin(atan(q.x, q.y) * repetitions) * (1. - smoothness)) * k;
  float ds = k * repetitions;
  float y = mix(ds, 2. * ds - length(vec2(x, ds)), correction);
  return vec2(x / repetitions, y / repetitions - displacement);
}

vec3 opElogate ( in vec3 q, in vec3 h, out float correction ) {
   q = abs(q) - h;

   correction = min(vmax(q), 0.);

   return max(q, 0.);
}

void quadrantIndex (in vec2 q, inout float carryIndex) {
  // --- Quadrant index ---
  // Assuming starting with 1.

  // -- y split --
  //  + | + 
  // ---|---
  //  - | - 
  carryIndex *= sign(q.y);

  // -- x split --
  //  x2 | x1 
  // ----|----
  //  x2 | x1 
  carryIndex *= q.x > 0. ? 1. : 2.;

  //  So: 
  //   2 |  1 
  // ----|----
  //  -2 | -1 
  // This should be okay to apply multiple times to have nested quadrants
}

vec2 opElogateS ( in vec2 p, in vec2 h ) {
   vec2 q = p - clamp(p, -h, h);

   return q;
}

float opExtrude ( in vec3 q, in float d, in float h ) {
  vec2 w = vec2(d, abs(q.z) - h);
  return min(max(w.x, w.y), 0.) + length(max(w, 0.));
}

#define jTrap 14
vec2 julia (in vec4 z, in vec4 c, in float t) {
  // Fractal General setup
  float minD = 1e10;
  float dist2dq = 1.;
  float modulo2;

  float avgD = 0.;
  const int iterations = 200;
  float dropOutInteration = float(iterations);
  float iteration = 0.;

  const float warpScale = 0.5;

  for (int i = 0; i < iterations; i++) {
    float fI = float(i);

    // // General space pre-warp
    // z += warpScale * 0.1000 * triangleWave(3. * z.yzwx + t * TWO_PI);

    // // Kifs
    // z = abs(z);
    // // z *= angle2C;;
    // // z.xy = abs(z.xy);
    // // z.x += angle3C;
    // z *= kifsM;
    // z.xyz += offset;
    // // z.yzwx = z.xyzw;
    // // z *= rotMat2(offset.x * PI + 0.125 * PI * sin(TWO_PI * t) + 0.3);

    // Julia set
    // z³ power
    // z' = 3q² -> |z'|² = 9|z²|²
    dist2dq *= 9. * qLength2(qSquare(z));
    z = qCube(z);

    // // z² power
    // // z' = 2q -> |z'|² = 4|z|²
    // dist2dq = 2. * modulo2;
    // z = qSquare(z);

    z += c;
    modulo2 = qLength2(z);

    // // Mandelbrot set
    // vec2 c = uv;
    // z = cSquare(z);
    // z += c;

    // if (i > 2) {
      // float trap = length(z);
      // float trap = dot(z, z);
      // float pr = 5.5;
      // float d = pow(dot(pow(z, vec2(pr)), vec2(1)), 1. / pr);
      float trap = length(z.xy - vec2(0.669, -0.323) + 0.4 * sin(z.zw + PI)); // circle trap
      // trap = 0.5 + 0.5 * sin(TWO_PI * trap);
      trap = abs(trap);
      // trap -= 0.00625 * iteration / float(iterations);
      trap -= 0.082;

    // float trap = lineTrap(z);

      avgD += trap;
      minD = min(minD, trap);
    // }

    float dis = modulo2;
    if (dis > 256.) break;
    if (iteration >= dropOutInteration) break;

    iteration += 1.;
  }

  avgD /= float(dropOutInteration);

  // Fractal Hubbard-Douady potential distance estimation
  // SDF(z) = log|z|·|z|/|dz| : https://iquilezles.org/www/articles/distancefractals/distancefractals.htm
  float sdfD = 0.25 * log(modulo2) * sqrt(modulo2 / dist2dq);
  return vec2(sdfD, minD);
}
vec2 julia (in vec4 z, in vec4 c) {
  return julia(z, c, 0.);
}
vec2 julia (in vec4 z) {
  vec4 c = vec4(0);
  return julia(z, c, 0.);
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

float onion (in float d, in float thickness) {
  return abs(d) - thickness;
}

mat3 rotOrtho (in float t) {
  const vec3 rotAxis = vec3(0, 1, 0);
  return rotationMatrix(rotAxis, 1.5 * PI * (0.5 + 0.5 * cos(t)));
}

// Create multiple copies of an object - http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
float opRepLim( in float p, in float s, in float lim ) {
  return p-s*clamp(floor(p/s + 0.5),-lim,lim);
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
vec2 opRevolution( in vec3 p, float o ) {
    return vec2( length(p.xz) - o, p.y );
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

vec4 pieSpace (in vec3 p, in float relativeC, in float repetitions) {
  float angle = relativeC * TWO_PI / repetitions;
  p.xz *= rotMat2(angle);
  float c = pModPolar(p.xz, repetitions);
  p.xz *= rotMat2(-angle);
  return vec4(p, c);
}

vec4 pieSpace (in vec3 p, in float relativeC) {
  const float repetitions = 6.;
  float angle = relativeC * TWO_PI / repetitions;
  p.xz *= rotMat2(angle);
  float c = pModPolar(p.xz, repetitions);
  p.xz *= rotMat2(-angle);
  return vec4(p, c);
}

float r = 1.50;
float sdHollowBox (in vec3 q, in vec3 r, in float thickness) {
  const float cropHeight = 30.;

  float b = sdBox(q, r);

  // crop inners
  vec3 cropR = r - thickness;
  vec3 cropQ = q;
  float crop = sdBox(cropQ, vec3(cropHeight, cropR.y, cropR.z));

  crop = min(crop, sdBox(cropQ, vec3(cropR.x, cropHeight, cropR.z)));

  crop = min(crop, sdBox(cropQ, vec3(cropR.x, cropR.y, cropHeight)));

  return max(b, -crop);
}

float sdHollowBox (in vec4 q, in vec4 r, in float thickness) {
  const float cropHeight = 30.;

  float b = sdBox(q, r);

  // crop inners
  vec4 cropR = r - thickness;
  vec4 cropQ = q;
  float crop = sdBox(cropQ, vec4(cropHeight, cropR.y, cropR.z, cropR.w));

  crop = min(crop, sdBox(cropQ, vec4(cropR.x, cropHeight, cropR.z, cropR.w)));

  crop = min(crop, sdBox(cropQ, vec4(cropR.x, cropR.y, cropHeight, cropR.w)));

  crop = min(crop, sdBox(cropQ, vec4(cropR.x, cropR.y, cropR.z, cropHeight)));

  return max(b, -crop);
}

vec2 sdBoxWEdges (in vec3 q, in vec3 r, in float thickness) {
  float b = sdBox(q, r);
  float m = 0.;

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

  if (-crop > b) {
    m = 1.;
  }
  return vec2(b, m);
}

float sdCross (in vec2 q, in float thickness) {
  return vmin(abs(q)) - thickness;
}

float sdCappedCross (in vec2 q, in vec2 r, in float thickness) {
  float cross = sdCross(q, thickness);
  return max(cross, sdBox(q, r));
}

float sdCappedCross (in vec2 q, in float r, in float thickness) {
  return sdCappedCross(q, vec2(r), thickness);
}

//--------------------------------------------------------------------------------
// ray-sphere intersection
// http://iquilezles.org/www/articles/intersectors/intersectors.htm
//--------------------------------------------------------------------------------
vec2 iSphere( in vec3 ro, in vec3 rd, in float rad ) {
  float b = dot( ro, rd );
  float c = dot( ro, ro ) - rad*rad;
  float h = b*b - c;
  if( h<0.0 ) return vec2(-1.0);
  h = sqrt(h);
  return vec2(-b-h, -b+h );
}

float gridBump ( in vec3 q, float size ) {
  vec3 c = pMod3(q, vec3(size));

  vec3 absQ = abs(q);
  return vmax(absQ) / size;
}

vec2 polarCoords (in vec2 q) {
  return vec2(atan(q.y, q.x), length(q));
}

vec2 polarToCartesian (in vec2 q) {
  float r = q.y;
  float a = q.x;
  return vec2(r, 0) * rotMat2(a);
}

vec3 sphericalCoords (in vec3 q) {
  float a = atan(q.z, q.x);
  float arc = atan(q.y, length(q.xz));
  return vec3(
      a,
      arc,
      length(q));
}

vec3 mobiusOriginPos (in vec3 q, in float bigR, in float yOffset, in float angle) {
  // q *= rotationMatrix(vec3(0, 0, -1), angle);
  q.xy *= rotMat2(-angle);
  q.x -= bigR;
  q.xz *= rotMat2(0.5 * angle);
  q.x -= yOffset;

  return q;
}

vec3 rotTorus (in vec3 q, in float angle, in float r) {
  float a = atan(q.z, q.x);
  q.xz *= rotMat2(-a);
  q.x -= r; // Move to center
  q.xy *= rotMat2(angle);
  q.x += r; // Move to Ring radius
  q.xz *= rotMat2(a); // Rotate back

  return q;
}

vec3 rotPlane (in vec3 q, in float ind, in float t) {
  float angle = 0.5 * PI * t * step(0.9, ind);
  q *= rotationMatrix(vec3(1, 0, 0), angle);
  return q;
}

vec2 foldAcross45s (in vec2 q) {
  q *= rotMat2(0.25 * PI);

  q = abs(q);

  q *= rotMat2(-0.25 * PI);

  // // Mirror around y=x
  // if (q.y >= q.x) {
  //   q.xy = q.yx;
  // }

  return q;
}

vec3 foldAcross45s (in vec3 q) {
  q = abs(q);

  // Mirror around y=x
  if (q.y >= q.x) {
    q.xy = q.yx;
  }
  // Mirror around z=x
  if (q.z >= q.x) {
    q.xz = q.zx;
  }

  return q;
}

vec2 lissajous (in float bigA, in float bigB, in float a, in float b, in float delta, in float t) {
  return vec2(bigA, bigB) * sin(vec2(a, b) * t + vec2(delta, 0));
}

vec2 split (in vec2 q, inout float mask, in float angle, in float gap, in float start) {
  vec2 axis = vec2(1, 0);
  axis *= rotMat2(angle);

  // Shape - extrude
  vec2 shapeQ = q;

  shapeQ += axis * start;

  shapeQ *= rotMat2(-angle);
  shapeQ = opElogateS( shapeQ, vec2(1, 0) * gap);
  shapeQ *= rotMat2( angle);

  shapeQ -= axis * start;

  // Mask
  float newMask = -abs(dot(q, axis) + start) + gap;
  mask = max(mask, newMask);
  return shapeQ;
}

const float numBreaks = 5.;
vec3 splitParams (in float i, in float t) {
  const float nullBuffer = 0.10;
  const float splitLength = (1. - 2. * nullBuffer) / numBreaks;

  // Distinct times per crack
  float localT = range( nullBuffer + splitLength * i, nullBuffer + splitLength * (i + 1.), t * (1. - smoothstep(0.8, 1., t)));

  // Slight overlap in time
  // float localT = range( nullBuffer + splitLength * 0.45 * i, nullBuffer + splitLength * (i + 2.0), triangleWave(t + 0.5) );
    // Collapse together
    // * (1. - range(0.7 + 0.002 * i, 0.8 + 0.001 * i, t));

  float gapAmount = 0.220 * (0.95 * mod(133.2830 * i, 1.0) + 0.05);
  float gap = gapAmount * quart(localT);

  float angle = snoise2(vec2(3.157143 * i) + 1.482349) * PI;

  float start = 0.15 * snoise2(vec2(10.8123 * i));

  return vec3(angle, gap, start);
}

const vec2 gSize = vec2(0.015);
float microGrid ( in vec2 q ) {
  vec2 cMini = pMod2(q, vec2(gSize * 0.10));

  return mod(dot(cMini, vec2(1)), 2.);
}

vec2 gC = vec2(0);
float localCosT = cosT;
float localT = norT;
float second = maxDistance;
vec2 shape (in vec2 q, in vec2 c) {
  vec2 d = vec2(maxDistance, -1.);

  vec2 uv = q;

  // float dC = vmax(abs(c));
  float dC = dot(c, vec2(-1, 1));

  float odd = mod(dC, 2.);
  float even = 1. - odd;

  // vec2 size = vec2(0.85, 0.15);

  // // Assume [0,1] range per dimension
  // vec2 bigSize = vec2(2);
  // vec2 bigC = floor(abs(c) / bigSize);
  // vec2 miniC = mod(c, bigSize);

  // Create a copy so there is no cross talk in neighborGrid
  float locallocalT = localT;
  // locallocalT = angle1C;
  // locallocalT -= 0.03 * length(c);
  // locallocalT -= 0.07 * vmax(abs(0.4 * c));
  // locallocalT -= 0.07 * vmax(vec2(0.4, 0.3) * c);
  // locallocalT -= atan(c.y, c.x) / PI;
  locallocalT += 0.0125 * dC;
  // locallocalT += 0.125 * snoise2(0.05 * (c + gC) + vec2(19.7, 113.1273));
  // locallocalT += 0.02 * odd;
  // locallocalT += 2.00 * q.x;
  // NOTE Flip time offset if there are gaps
  // Might fix some of the gaps caused by the time offset
  // A hack but getting closer to a general solution
  // const float clip = 0.0125;
  // locallocalT = clamp(locallocalT, clip, 1. - clip);

  float t = mod(locallocalT, 1.);
  // t = expo(t);
  float localCosT = TWO_PI * t;

  // Local C that transitions from one cell to another
  float shift = 0.;
  vec2 shiftDir = vec2(1, 1);

  vec2 localC = mix(c, c + shift * shiftDir, t);

  // // Vanilla cell coordinate
  // vec2 localC = c;

  float localNorT = 0.5 + 0.5 * cos(localCosT);
  float warpScale = 0.45 * expo(localNorT);

  vec2 size = gSize;
  vec2 r = 0.225 * size;

  // q.x += 0.5 * localNorT * size.x * mod(localC.y, 2.);

  // Make grid look like random placement
  float nT = 0.5 + 0.5 * sin(localCosT); // 0.5; // triangleWave(t);
  q += 1.5 * localNorT * size.x * mix(
      0.2 * vec2(1, -1) * snoise2(2.417 * localC + 73.17123),
      vec2(1) * snoise2(8.863 * localC + 2.37),
      nT);

  // float side = step(abs(c.y), abs(c.x));
  // q.x += sign(c.x) * side * size.x * (0.5 + 0.5 * cos(localCosT));

  // q.x += 1.1 * size.x * (0.5 + 0.5 * cos(localCosT));

  // q.x += t * size.x * mod((shift * shiftDir).y, 2.);

  // q.x += size.x * (1. - 2. * mod(c.y, 2.)) * (0.5 + 0.5 * cos(localCosT + 0.2 * c.x));

  t -= 0.45;
  t = mod(t, 1.);
  q += 3. * size * shiftDir * t;

  // vec2 center = vec2(size.x * (c + gC));
  // center += size.x * warpScale * 0.10000 * cos( 3.17823 * center.yx + localCosT + vec2(9.2378));
  // center += size.x * warpScale * 0.05000 * cos( 7.91230 * center.yx + localCosT + vec2(-10.2378));
  // center *= rotMat2(0.005 * PI * cos(localCosT - length(0.1 * c)));
  // center += size.x * warpScale * 0.02500 * cos(13.71347 * center.yx + localCosT);
  // center -= size.x * c;
  // q += center;

  // // Cosine warp
  // float warpScale2 = warpScale * 0.5;
  // q += vec2(-1, 1) * warpScale2 * 0.10000 * cos( 2. * vec2(-1, 1) * q.yx + localCosT + 0.71283 * gC);
  // q += vec2(-1, 1) * warpScale2 * 0.05000 * cos( 3. * vec2(-1, 1) * q.yx + localCosT + 0.91283 * gC);
  // q += vec2(-1, 1) * warpScale2 * 0.02500 * cos( 5. * vec2(-1, 1) * q.yx + localCosT + 1.11283 * gC);
  // q += vec2(-1, 1) * warpScale2 * 0.01250 * cos( 7. * vec2(-1, 1) * q.yx + localCosT - 0.71283 * gC);
  // q += vec2(-1, 1) * warpScale2 * 0.00625 * cos(11. * vec2(-1, 1) * q.yx + localCosT + 0.31283 * gC);

  // q *= rotMat2(0.5 * PI * cos(localCosT));

  // c = floor((q + 0.5 * size) / size);

  // q.x += 0.333 * size.x * mod(c.y, 3.);
  // c = pMod2(q, size);

  // q -= shiftDir * shift * size * t;

  // // Rotate randomly
  // q *= rotMat2(1.0 * PI * snoise2(0.263 * localC));

  float internalD = length(q) - r.x;
  // float internalD = abs(q.y);
  // internalD = max(internalD, abs(q.x) - 0.7 * vmax(size));
  // internalD = min(internalD, abs(q.x));
  // internalD = max(internalD, sdBox(q, 0.5 * size));

  // float internalD = abs(dot(q, vec2(-1, 1)));
  // internalD = max(internalD, sdBox(q, vec2(0.5 * size)));
  // float internalD = vmax(abs(q));
  // float internalD = dot(abs(q), vec2(1));
  // float internalD = sdBox(q, r);
  // internalD = abs(internalD) - 0.05 * vmax(r);

  // vec2 absQ = abs(q);
  // float internalD = min(absQ.x, absQ.y);
  // float crossMask = sdBox(q, vec2(0.35 * size));
  // internalD = max(internalD, crossMask);

  // float internalD = sdHexPrism(vec3(q, 0), vec2(r, 1.));

  // // Arrow Up
  // vec2 arrowQ = q;
  // arrowQ.y += abs(arrowQ.x);
  // float internalD = abs(arrowQ.y) - 0.1 * size;
  // internalD = max(internalD, sdBox(q, vec2(r)));

  // // -- 2D SDF Texture --
  // // Convert from center 0 to center 0.5
  // q += 0.5;
  // float internalD = texture2D(sdf2DTexture, q).r;
  // // Unpack from [0, 1] to [-1, 1]
  // internalD -= 0.5;
  // internalD *= 2.;

  // float internalD = sdBox(q, r);

  vec2 o = vec2(internalD, 0.);
  // vec2 o = vec2(internalD - 0.03 * size.x, 0.);
  // float o = microGrid(q);
  d = dMin(d, o);

  // // Outline
  // const float adjustment = 0.0;
  // d = abs(d - adjustment) - r * 0.025;

  // Mask
  float mask = 0.;
  // mask = step(0., dot(abs(c), vec2(1)) - 12.);
  // mask = step(0., vmax(abs(c)) - 12.);
  // mask = step(0., abs(sdBox(c, vec2(26))) - 8.);
  // mask = step(0., sdBox(q, size * vec2(1, 5)));
  // mask = step(0., abs(length(c) - 4.) - 2.));
  // mask = step(0., length(c) - 35.);
  // mask = step(0., abs(c.y - 25.) - 8.); // Mask below a line
  // Convert circle into torus
  // mask = step(0., abs(length(c) - 16.) - 10.);

  // Apply mask
  d.x = mix(d.x, maxDistance, mask);

  return d;
}

// NOTE Don't fully understand whether this is correct or how it works
vec2 circleInversion (in vec2 q) {
  return q / dot(q, q);
  // q is point in the circle & a is the inverse
  // dot(q, a) = r * r
  // q.x * a.x + q.y * a.y = r * r // i don't know what the invert of a dot product is...
}

#pragma glslify: neighborGrid = require(./modulo/neighbor-grid, map=shape, maxDistance=maxDistance, numberOfNeighbors=3.)

float thingy (in vec2 q, in float t) {
  float d = maxDistance;

  vec2 uv = q;

  localCosT = TWO_PI * t;
  localT = t;

  float bigR = r * 1.5;

  for (float i = 0.; i < 3.; i++) {
    vec2 localQ = q;
    localQ += lissajous(bigR, bigR, 1., 2., PI * 0.5, localCosT + TWO_PI * 0.333 * i);
    float b = length(localQ) - r;
    d = min(d, b);
  }

  // // Outline
  // const float adjustment = 0.0;
  // d = abs(d - adjustment) - 0.02 / rollingScale;

  float stop = angle3C;
  // d = smoothstep(stop, 0.3 * edge + stop, d);
  // d = 1. - d;

  return d;
}

float sdHalfDome (in vec3 q, in float r, in float thickness) {
  float d = length(q) - r;
  d = onion(d, thickness);
  d = max(d, q.x);
  d *= 0.125;
  return d;
}

// A box hollow on one side
float sdBin (in vec3 q, in vec3 r, in float thickness) {
  float b = sdBox(q, r);

  // crop
  vec2 innerR = r.xy - thickness;
  float longR = 2.0 * r.z;
  float crop = sdBox(q - vec3(0, 0, r.z + thickness), vec3(innerR, longR));
  b = max(b, -crop);

  return b;
}

float arrowUpTexture (in vec2 q, in float size) {
  float r = 0.4 * size;

  // Arrow Up texture
  vec2 arrowQ = q;
  vec2 c = pMod2(arrowQ, vec2(size));
  vec2 localQ = arrowQ;
  arrowQ.y += abs(arrowQ.x);
  arrowQ.y -= size * 0.25;
  float internalD = abs(arrowQ.y) - 0.1 * size;
  return max(internalD, sdBox(localQ, vec2(r)));
}

vec2 piston (in vec3 q, in float r, in float t) {
  vec2 d = vec2(maxDistance, -1.);

  const float headBodyRatio = 0.25;
  const float invHeadBodyRatio = 1. - headBodyRatio;

  // Assume pointed in the +x direction
  vec3 headQ = q;
  headQ.x += invHeadBodyRatio * r + t * 2. * r;

  float head = sdBox(headQ, 0.95 * vec3(headBodyRatio * r, r, r));
  d = dMin(d, vec2(head, 0));

  // Shaft
  float shaftR = 0.2 * r;
  float shaftLength = t * r;

  vec3 shaftQ = q.yxz;
  shaftQ.y += shaftLength + 2. * (invHeadBodyRatio - 0.5) * r;
  float shaft = sdCappedCylinder(shaftQ, vec2(shaftR, shaftLength));
  d = dMin(d, vec2(shaft, 2));

  vec3 bodyQ = q;
  bodyQ.x -= headBodyRatio * r;
  float body = sdBox(bodyQ, vec3(invHeadBodyRatio * r, r, r));
  d = dMin(d, vec2(body, 1));

  return d;
}

float dotTexture (in vec2 q, in float size) {
  float r = 0.4 * size;

  // Arrow Up texture
  vec2 arrowQ = q;
  vec2 c = pMod2(arrowQ, vec2(size));
  vec2 localQ = arrowQ;
  float internalD = abs(length(localQ) - 0.4 * size) - 0.05 * size;
  return max(internalD, sdBox(localQ, vec2(r)));
}

// IQ's 2D isosceles triangle
// source: https://www.shadertoy.com/view/MldcD7
float sdTriangleIsosceles( in vec2 p, in vec2 q ) {
  p.x = abs(p.x);
  vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(p.x*q.y-p.y*q.x),k*(p.y-q.y)  );
  return sqrt(d)*sign(s);
}

// IQ's 2D triangle
// source: https://www.shadertoy.com/view/XsXSz4
float sdTriangle( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 ) {
  vec2 e0 = p1 - p0;
  vec2 e1 = p2 - p1;
  vec2 e2 = p0 - p2;

  vec2 v0 = p - p0;
  vec2 v1 = p - p1;
  vec2 v2 = p - p2;

  vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );
  vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
  vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );

  float s = e0.x*e2.y - e0.y*e2.x;
  vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),
        vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),
      vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));

  return -sqrt(d.x)*sign(d.y);
}

float gear (in vec3 p, in float r, in float thickness, in float thinness, in float teeth) {
  vec3 q = p;
  float d = maxDistance;
  float b = sdCappedCylinder(q, vec2(r, thinness));
  d = min(d, b);

  pModPolar(q.xz, teeth);
  float toothThickness = thickness * 0.4;
  q.x -= r + 1.0 * toothThickness;
  float cog = sdBox(q, vec3(toothThickness, thinness, toothThickness));
  d = min(d, cog);

  q = p;
  float cutOut = sdCappedCylinder(q, vec2(r - thickness, 1.));
  d = max(d, -cutOut);

  pModPolar(q.xz, 6.);
  float spoke = sdBox(q, vec3(r, thinness, thickness));
  d = min(d, spoke);

  return d;
}

float axialStar (in vec3 q, in float r, in float thickness) {
  // Create plus sign
  if (abs(q.z) > abs(q.x)) {
    q.zx = q.xz;
  }

  return sdBox(q, vec3(r, vec2(thickness)));
}

vec2 conveyerBelt (in vec3 q, in vec3 beltDims, in float thickness, in float t) {
  vec2 d = vec2(maxDistance, -1);

  float beltSpacing = beltDims.y;
  float ridgeSpacing = 0.15 * beltDims.x;

  vec3 beltTopQ = vec3(q);
  beltTopQ.y = abs(beltTopQ.y) - beltSpacing;
  vec2 beltTop = vec2(sdBox(beltTopQ, vec3(beltDims.x, thickness, beltDims.z)), 0);

  float ridgeDepth = 0.25 * thickness;
  vec3 ridgesQ = beltTopQ;

  ridgesQ.x += sign(q.y) * ridgeSpacing * t;

  pMod1(ridgesQ.x, ridgeSpacing);
  ridgesQ.y -= thickness;
  float ridges = sdBox(ridgesQ, vec3(ridgeDepth, ridgeDepth, beltDims.z * 2.));
  beltTop.x = max(beltTop.x, - ridges);
  d.xy = dMin(d.xy, beltTop);

  // axis
  float rollerR = beltSpacing - thickness;
  vec3 axisQ = q;

  pMod1(axisQ.x, 3.5 * rollerR);
  axisQ.zy = axisQ.yz;

  vec2 roller = vec2(sdCappedCylinder(axisQ, vec2(rollerR, beltDims.z * 0.96)), 1);
  d.xy = dMin(d.xy, roller);
  vec2 axis = vec2(sdCappedCylinder(axisQ, vec2(rollerR * 0.5, beltDims.z * 1.02)), 2);
  d.xy = dMin(d.xy, axis);

  return d;
}

#pragma glslify: loopNoise = require(./loop-noise, noise=snoise3)

float crystal (in vec3 q, in float r, in vec3 h, in float angle) {
  float d = maxDistance;

  float o = 0.;
  q = opElogate(q, h, o);

  float b = icosahedral(q, 42., r);
  d = min(d, b);

  q *= rotationMatrix(vec3(1), angle);
  float crop = dodecahedral(q, 42., 0.95 * r);
  d = max(d, crop);

  return d;
}

float crossGyroid (in vec3 p, in float thickness) {
  vec3 xross = cross(sin(p), cos(p.yzx));
  float gyroid = dot(xross, xross);
  return abs(gyroid) - thickness;
}

vec4 componentShift (in vec4 q) {
  return q.yzwx;
}

vec3 componentShift (in vec3 q) {
  return q.yzx;
}

vec2 componentShift (in vec2 q) {
  return q.yx;
}

float gyroid (in vec4 p, in float thickness) {
  float gyroid = dot(sin(p), cos(componentShift(p)));
  return abs(gyroid) - thickness;
}

float squiggle ( in vec2 q, in float r, in float thickness ) {
  float d = maxDistance;

  float a1 = 0.25 * PI;
  float a2 = 0.25 * PI;
  float n = sdArc(q, vec2(sin(a1), cos(a1)), vec2(sin(a2), cos(a2)), r, thickness);
  d = min(d, n);

  for (int i = 0; i < 27; i++) {
    a1 -= -1. * PI;
    q.y -= 2. * r;
    n = sdArc(q, vec2(sin(a1), cos(a1)), vec2(sin(a2), cos(a2)), r, thickness);
    d = min(d, n);

    a1 -= 1.00 * PI;
    q.x -= 2. * r;
    n = sdArc(q, vec2(sin(a1), cos(a1)), vec2(sin(a2), cos(a2)), r, thickness);
    d = min(d, n);
  }

  return d;
}

float squiggle (in vec3 q, in float r, in float thickness) {
  float d = maxDistance;

  vec3 originalQ = q;

  float c = pMod1(q.x, 4. * r);
  float b = sdTorus(q.xzy, vec2(r, thickness));
  b = max(b, max(0., -q.y));
  d = min(d, b);

  q = originalQ;
  q.x += 2. * r;
  pMod1(q.x, 4. * r);
  b = sdTorus(q.xzy, vec2(r, thickness));
  b = max(b, max(0., q.y));
  d = min(d, b);

  return d;
}


#define splitTime(x, t) range(prevT, prevT + x, t); prevT += x;

float tileTime (in vec2 c, in float t) {
  vec2 cTDelay = -0.25 * vec2(0.02, 0.05);

  return dot(c, cTDelay);
}

float tile (in vec3 q, in vec2 c, in float r, in vec2 size, in float t) {
  float d = maxDistance;

  float originalR = r;

  q.xz += c * size;
  // q.z += 2. * r * t; // Recenter over time

  float prevT = 0.;
  float growT = splitTime(0.25, t);
  growT = expo(growT);
  float fallT = splitTime(0.25, t);
  fallT = expoIn(fallT);

  float splitThreshold = step(prevT, t);
  float splitT = splitTime(0.5, t);
  splitT = expoOut(splitT);

  float h = r * (1. + growT);
  h = mix(h, r, splitThreshold);

  // Squeeze when stretching up
  r -= 0.2 * originalR * growT;

  // Relax
  h = mix(h, 0.95 * r, fallT);
  r = mix(r, 3. * originalR, fallT);

  // Split
  r = mix(r, originalR, splitThreshold);
  h = mix(h, originalR, splitT);

  vec3 boxQ = q;
  boxQ.y -= h; // Offset growing height

  float b = sdBox(boxQ, vec3(r, h, r));
  d = min(d, b);

  vec3 direction = vec3(0);
  vec3 splitBoxQ = boxQ;

  float splitB = maxDistance;

  direction = vec3( 0, 0, 1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3( 0, 0,-1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3( 1, 0, 0);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3(-1, 0, 0);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3( 1, 0, 1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3(-1, 0, 1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3( 1, 0,-1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  direction = vec3(-1, 0,-1);
  splitBoxQ = boxQ;
  splitBoxQ -= 2. * r * direction * splitThreshold;
  splitBoxQ -= 2. * r * direction * splitT;

  splitB = sdBox(splitBoxQ, vec3(r, h, r));
  d = min(d, splitB);

  return d;
}

#pragma glslify: subdivide = require(./modulo/subdivide.glsl, vmin=vmin, noise=h21)

float gR = 0.25;
bool isDispersion = false;
bool isSoftShadow = false;
vec3 map (in vec3 p, in float dT, in float universe) {
  vec3 d = vec3(maxDistance, 0, 0);
  vec2 minD = vec2(1e19, 0);

  float t = mod(dT, 1.);
  float localCosT = TWO_PI * t;
  float r = gR;
  vec2 size = r * vec2(4.75);

  // Positioning adjustments

  // // -- Pseudo Camera Movement --
  // // Wobble Tilt
  // const float tilt = 0.08 * PI;
  // p *= rotationMatrix(vec3(1, 0, 0), 0.25 * tilt * cos(localCosT));
  // p *= rotationMatrix(vec3(0, 1, 0), 0.2 * tilt * sin(localCosT - 0.2 * PI));

  // p *= globalRot;

  vec3 q = p;

  float warpScale = 0.25;
  float warpFrequency = 0.2;
  float rollingScale = 1.;

  // Warp
  vec3 preWarpQ = q;
  vec3 wQ = q.xyz;

  // vec4 wQ = vec4(q.xyz, 0);

#define distortT localCosT

  // float worldScale = 1.0;
  // wQ *= worldScale;

  float phasePeriod = 0.5 * (0.5 + 0.5 * cos(dot(wQ, vec3(1)) + localCosT));
  vec3 warpPhase = TWO_PI * phasePeriod * vec3(0., 0.33, 0.67) + 0.9;
  // vec4 warpPhase = TWO_PI * phasePeriod * vec4(0., 0.25, 0.5, 0.75) + 0.9;

  const float warpPhaseAmp = 0.9;

  wQ += 0.100000 * warpScale * cos( 2.182 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  wQ += 0.050000 * warpScale * cos( 5.732 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  warpPhase += warpPhaseAmp * componentShift(wQ);
  wQ.xyz = twist(wQ.xzy, 1. * wQ.z + 0.2 * PI * cos(localCosT + 0.9 * wQ.z));
  wQ += 0.025000 * warpScale * cos( 9.123 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  wQ += 0.012500 * warpScale * cos(13.923 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  warpPhase += warpPhaseAmp * componentShift(wQ);
  // wQ.xyz = twist(wQ.xzy, 0.35 * wQ.x + 0.305 * sin(localCosT + wQ.x));
  wQ += 0.006250 * warpScale * cos(17.369 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  wQ += 0.003125 * warpScale * cos(19.937 * warpFrequency * componentShift(wQ) + distortT + warpPhase);
  warpPhase += warpPhaseAmp * componentShift(wQ);
  wQ += 0.001125 * warpScale * cos(23.937 * warpFrequency * componentShift(wQ) + distortT + warpPhase);

  // Commit warp
  q = wQ.xyz;
  mPos = q;

  r *= 1. - saturate(-q.z / 2.5);

  vec3 b = vec3(r - length(q.xy), 0, 0);
  d = dMin(d, b);

  // // Fractal Scale compensation
  // d.x /= rollingScale;

  // // Scale compensation
  // d.x /= worldScale;

  // // Under step
  // d.x *= 0.70;

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
  float minM = maxDistance;

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
    maxI = float(i);
    if (d.x < epsilon) return vec4(t + d.x, d.y, maxI, d.z);
    t += d.x;
    trap = min(trap, d.z);
    if (d.y == 1.) { // Select a material id to track
      minM = min(minM, t);
    }
    if (t > maxDistance) break;
    deltaT += deltaTDelta;
  }
#endif
  return vec4(-1., minM, maxI, trap);
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
#pragma glslify: simpleshadow = require(./lighting/simple-shadow.glsl, map=map)
#pragma glslify: calcAO = require(./ao, map=map)

// --- Patterns ---
// #pragma glslify: checker = require(glsl-checker)
// #pragma glslify: herringBone = require(./patterns/herring-bone)

#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: rgb2hsv = require(./rgb2hsv.glsl)
#pragma glslify: hsb2rgb = require(./color-map/hsb2rgb)

float gM = 0.;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float dNR = dot(-rd, gNor);

  float spread = saturate(1.0 - 1.0 * pow(dNR, 9.));
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
  dI += 0.3 * snoise3(0.3 * mPos);
  dI += 0.5 * pow(dNR, 5.);

  dI *= angle1C;
  dI += angle2C;
  // dI += norT;

  // dI += gC.z;

  // dI += 0.2 * gPos;

  // dI *= 0.3;

  // -- Colors --
  color = 0.5 + 0.5 * cos( TWO_PI * ( vec3(1) * dI + vec3(0, 0.33, 0.67) ) );
  // color = mix(#FF0000, #00FFFF, 0.5 + 0.5 * sin(TWO_PI * (dI)));

  // // - Rotated Components -
  // float angle = 1.0 * gPos.y;
  // mat3 rot = rotationMatrix(vec3(1), angle);

  // color = vec3(0);
  // color += vec3(1, 0, 1) * rot * snoise3(0.2 * dI);
  // color += vec3(1, 0, 1) * rot * -dNR;
  // color += vec3(1, 1, 0) * rot * cos(TWO_PI * (snoise3(0.3 * gPos) + vec3(0, 0.33, 0.67)));

  // color *= 0.6;

  // color = vec3(n);

  color *= spread;

  // color = getBackground(rd.xy, 0.);

  // // Identity scene color
  // color = vec3(1);

  // color *= 1.3;

  // // Unbounded
  // return color;

  // Clamped
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

  for (int i = 0; i < 64; i++) {
    float d = dispersionMap(rayOrigin + rayDirection * t, norT).x;
    if (d >= 0.0) break;
    d = min(d, -0.0001);

    t += abs(d);
  }
  return t;
}

vec3 secondRefraction (in vec3 rd, in float ior) {
  float d = 0.0;

  // Perturbate initial ray direction
  float angle = ncnoise3(rd);
  rd *= rotationMatrix(vec3(1, 4, 2), smoothstep(0.5, 0.6, angle));

  #if 1
  d = dispersionMarch(rd);
  #else
  const int samples = 12;
  for (int i = 0; i < samples; i++) {
    vec3 lightDir = rd;

    // Add noise to walk to create imperfections
    lightDir += 0.1 * vec3(
        noise(gPos),
        noise(1.3 * gPos + 340.0),
        noise(-1.9 * gPos + 640.0));

    d += dispersionMarch(lightDir);
  }
  d /= float(samples);
  #endif

  vec3 reflectionPoint = gPos - gNor * 0.1 + rd * d;
  vec3 reflectionPointNor = getNormal2(reflectionPoint, 0.001, norT);
  dNor = reflectionPointNor;
  dRd = rd;
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
  vec3 color = vec3(0.);

  vec2 pol = vec2(
      atan(mPos.y, mPos.x),
      mPos.z);
  pol.x /= PI;

  float baseSize = 0.025;
  float lSize = 12.;
  vec2 size = baseSize * vec2(1. + floor(max(0., -(pol.y + 0.5 * 12. * baseSize) / (12. * baseSize))), 12.);
  vec2 r = vec2(0.1, 0.4) * size;

  vec2 c = pMod2(pol, size);

  pol.y += 0.15 * size.y * snoise2(1.238 * c);

  float b = sdBox(pol, r);
  float n = step(b, 0.);

  color = vec3(n);

  return color;

  // float n = dot(mPos.xyz, vec3(1));
  // n *= TWO_PI;
  // n *= 40.;
  // n = sin(n);
  // n += 0.6;
  // n = smoothstep(0., edge, n);
  // n *= 1.4;
  // return vec3(n);

  float dNR = dot(nor, -rd);
  vec3 dI = 0.3 * vec3(dot(nor, vec3(1)));
  // dI += 2. * pow(dNR, 2.);
  dI.xy += 0.1 * fragCoord.xy;

  dI += 0.0125 * pos;
  dI += 0.5 * snoise3(0.3 * mPos);

  dI *= angle1C;
  dI += angle2C;

  color = vec3(0.5) + vec3(0.5) * cos(TWO_PI * (vec3(1) * dI + vec3(0.0, 0.33, 0.67)));

  // float angle = 20.13 * PI + 0.8 * pos.y;
  // mat3 rot = rotationMatrix(vec3(1), angle);

  // // color = vec3(0);
  // color += vec3(1, 0, 0) * rot * cos(dI);
  // color += vec3(0, 1, 0) * rot * dNR;
  // color += vec3(0, 0, 1) * rot * 2. * snoise3(0.3 * pos);

  // color *= 0.95;

  // // -- Holo --
  // vec3 beforeColor = color;

  // const float numSteps = 30.;
  // const float stepSize = 0.20;
  // const float holoIoR = 1.0;
  // vec3 holoQ = pos;
  // // vec3 holoRd = refract(nor, rd, holoIoR);
  // vec3 holoRd = rd;
  // holoRd += 0.2 * refract(nor, rd, holoIoR);
  // for (float i = 0.; i < numSteps; i++) {
  //   vec3 holoPos = holoQ + i * stepSize * holoRd;
  //   // float inclusion = snoise3(0.5 * holoPos);
  //   vec3 s = vec3(0);
  //   float inclusion = fbmWarp(0.25 * vec3(1, 2, 1) * holoPos, s);
  //   // float inclusion = snoise3(1.225 * vec3(1) * holoPos + length(holoPos));
  //   // inclusion = step(0.25, inclusion);

  //   float colorSpeed = 0.5;

  //   // // Basic layer
  //   // dI = vec3(0.05 * i);

  //   // Single Noise
  //   dI = vec3(snoise3(vec3(1, 2, 2) * holoPos + 0.10 * i));

  //   // // Multi-noise
  //   // vec3 nQ = holoPos;
  //   // // nQ *= rotationMatrix(vec3(1), length(holoPos));

  //   // dI = vec3(
  //   //     snoise3(length(nQ) + colorSpeed * vec3(0.8, 1, 1.1) * nQ + 0.05 * i),
  //   //     snoise3(length(nQ) + colorSpeed * vec3(0.99, 1, 0.7) * nQ + 0.08 * i),
  //   //     snoise3(length(nQ) + colorSpeed * vec3(1, 0.2, 1.1) * nQ + 0.15 * i));

  //   // dI += s;

  //   dI *= angle1C;
  //   dI += angle2C;

  //   vec3 layerColor = vec3(0.5) + vec3(0.4, 0.5, 0.5) * cos(TWO_PI * (vec3(2, 1, 0.8) * dI + vec3(0,0.2, 0.4)));
  //   // color += saturate(inclusion) * layerColor;
  //   // color = mix(color, layerColor, saturate(inclusion));
  //   // color = mix(pow(color, vec3(2.2)), pow(layerColor, vec3(2.2)), saturate(inclusion));
  //   color = pow(color, vec3(0.454545));
  // }

  // color /= pow(numSteps, 0.30);
  // color *= 1.2;
  // color /= numSteps;

  // color.r = pow(color.r, 0.7);
  // color.b = pow(color.b, 0.7);

  // color += 0.4 * beforeColor;

  color += 0.4 * background;

  // color = mix(color, vec3(0.5), 0.2);
  // color = mix(color, vec3(1), 0.4);

  // color *= 0.5 + 0.5 * dNR;

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

    const float universe = 0.;
    background = getBackground(uv, universe);

    // Declare lights
    struct light {
      vec3 position;
      vec3 color;
      float intensity;
      float size;
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

    // // Test light
    // lights[0] = light(vec3(0.01,  1.0, 0.1), #FFFFFF, 1.0, 32.);

    lights[0] = light(2. * vec3(-0.9, 0.4, 1.0), #FFDCDC, 1.0, 32.0);
    lights[1] = light(2. * vec3( 0.6, 0.7, 0.8), #DCFFFF, 1.0, 16.0);
    lights[2] = light(2. * vec3( 0.2, 0.8,-1.3), #FFFFFF, 2., 16.0);

    float m = step(0., sin(TWO_PI * (0.25 * fragCoord.x + generalT)));

    float backgroundMask = 1.;
    // Allow anything in top right corner
    // backgroundMask = max(backgroundMask, smoothstep(0., edge, dot(uv, vec2(1)) + 0.05));

    if (t.x>0. && backgroundMask > 0.) {
      vec3 color = vec3(0.0);

      // Normals
      vec3 nor = getNormal2(pos, 0.001 * t.x, generalT);
      // float bumpsScale = 1.0;
      // float bumpIntensity = 0.075;
      // nor += bumpIntensity * vec3(
      //     snoise3(bumpsScale * 490.0 * mPos),
      //     snoise3(bumpsScale * 670.0 * mPos + 234.634),
      //     snoise3(bumpsScale * 310.0 * mPos + 23.4634));
      // nor -= 0.125 * cellular(5. * mPos);

      // // Cellular bump map
      // nor += 0.3 * (0.5 + 0.5 * dot(nor, rayDirection)) * cellular(vec3(9, 1, 9) * mPos);

      nor = normalize(nor);
      gNor = nor;

      vec3 ref = reflect(rayDirection, nor);
      ref = normalize(ref);

      gRd = rayDirection;

      // Basic Diffusion
      vec3 diffuseColor = baseColor(pos, nor, rayDirection, t.y, t.w, generalT);

      // Material Types
      float isFloor = isMaterialSmooth(t.y, 1.);

      float occ = calcAO(pos, nor, generalT);
      float amb = saturate(0.5 + 0.5 * nor.y);
      float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 0.0;
      float specCo = 0.0;

      vec3 specAll = vec3(0.0);

      // Shadow minimums
      float diffMin = 1.;
      float shadowMin = 1.;

      vec3 directLighting = vec3(0);
      for (int i = 0; i < NUM_OF_LIGHTS; i++) {
        vec3 lightPos = lights[i].position;
        vec3 nLightPos = normalize(lightPos);
        vec3 lightRd = normalize(lightPos - pos);

        float dif = max(diffMin, diffuse(nor, nLightPos));

        // // "Dithered" shadows
        // vec2 ditherQ = fragCoord.xy;
        // const float ditherSize = 0.0065;
        // vec2 ditherC = pMod2(ditherQ, vec2(ditherSize));
        // // float dither = length(ditherQ) - ditherSize * 0.25;
        // float dither = vmax(abs(ditherQ)) - ditherSize * 0.125;
        // float ditherAmount = 0.3 + 0.7 * range(0., 0.5 * ditherSize, dither);
        // dif = mix(1., ditherAmount, 1. - step(0.1, diffuse(nor, nLightPos)));

        float spec = pow(saturate(dot(ref, nLightPos)), 128.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        isSoftShadow = true;
        float sha = max(shadowMin, softshadow(pos, lightRd, 0.001, 1.0, lights[i].size, generalT));
        isSoftShadow = false;
        dif *= sha;

        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * occ * sha;
        lin += fre; // Commit Fresnel
        specAll += mix(lights[i].color, vec3(1), 0.2) * specCo * spec * sha;

        // // Ambient
        // lin += mix(0.0750, 0., isFloor) * amb * diffuseColor;
        // dif += mix(0.0750, 0., isFloor) * amb;

        float distIntensity = 1.; // lights[i].intensity / pow(length(lightPos - gPos), 0.55);
        distIntensity = saturate(distIntensity);
        color +=
          (dif * distIntensity) * lights[i].color * diffuseColor
          + distIntensity * mix(lights[i].color, vec3(1), 0.2) * lin * mix(diffuseColor, vec3(1), 1.0);

        // // Debug Light(s)
        // color += sha; // Soft shadows only
        // // color += dif; // Lambert shadows only

        // // -- Add in light flare --
        // vec3 fromLight = rayOrigin - lightPos;
        // float lightMasked = 1. - smoothstep(t.x, t.x + 0.001, length(fromLight));
        // float lightAngle = pow(max(0., dot(-rayDirection, normalize(fromLight))), 512.0);
        // directLighting +=
        //     lightMasked
        //   * mix(lights[i].color, vec3(1), 0.9 * lightAngle)
        //   * lightAngle;
      }

      color *= 1.0 / float(NUM_OF_LIGHTS);
      color += 1.0 * pow(specAll, vec3(8.0));

      // // Reflect scene
      // vec3 reflectColor = vec3(0);
      // vec3 reflectionRd = reflect(rayDirection, nor);
      // reflectColor += 0.3 * mix(diffuseColor, vec3(1), 1.0) * reflection(pos, reflectionRd, generalT);
      // color += reflectColor;

      // vec3 refractColor = vec3(0);
      // vec3 refractionRd = refract(rayDirection, nor, 1.5);
      // refractColor += 0.10 * textures(refractionRd);
      // color += refractColor;

#ifndef NO_MATERIALS

// -- Dispersion --
// #define useDispersion 1

#ifdef useDispersion
      // Set Global(s)
      dNor = gNor;

      isDispersion = true; // Set mode to dispersion

      // vec3 dispersionColor = dispersionStep1(nor, normalize(rayDirection), n2, n1);
      vec3 dispersionColor = dispersion(nor, rayDirection, n2, n1);

      isDispersion = false; // Unset dispersion mode

      float dispersionI = 1.0 * pow(0. + dot(dNor, -gRd), 2.);
      // float dispersionI = 1.0;
      dispersionI *= 2.;

      dispersionColor *= dispersionI;

      // Dispersion color post processing
      // dispersionColor.r = pow(dispersionColor.r, 0.7);
      // dispersionColor.b = pow(dispersionColor.b, 0.7);
      // dispersionColor.g = pow(dispersionColor.g, 0.8);

      // dispersionColor *= 0.9;

      color += saturate(dispersionColor);
      // color = mix(color, dispersionColor, saturate(pow(dot(dNor, -gRd), 1.5)));
      // color = saturate(dispersionColor);
      // color = vec3(dispersionI);
#endif

#endif

      // // Fog
      // float d = max(0.0, t.x);
      // color = mix(background, color, saturate(
      //       pow(clamp(fogMaxDistance - d, 0., fogMaxDistance), 1.2)
      //       / fogMaxDistance
      // ));
      // color *= saturate(exp(-d * 0.025));
      // color = mix(background, color, saturate(exp(-d * 0.05)));

      // color += directLighting * exp(-d * 0.0005);

      // Inner Glow
      // color += 0.5 * innerGlow(5.0 * t.w);

      // // Fade to background
      // color = mix(color, background, saturate(isFloor * pow(0.5 * t.w, 1.1)));

      // color = diffuseColor;

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

      // // Cut off shading
      // color = vec3(step(0.3, length(color)));

      // // Levels adjustment of sorts
      // color = mix(color, sigmoid(color), 0.7);

      // Max Call based material
      float stepIndex = t.z / float(maxSteps);
      stepIndex = pow(stepIndex, 0.85);

      // color += 0.5 + 0.5 * cos(TWO_PI * (2. * stepIndex + vec3(0, 0.33, 0.67) + 0.125));
      // color += 0.2 * vec3(stepIndex);

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

      // // -- Direct lighting --
      // for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
      //   vec3 lightPos = lights[i].position;
      //   vec3 fromLight = rayOrigin - lightPos;
      //   float lightMasked = 1. - smoothstep(t.x, t.x + 0.001, length(fromLight));
      //   float lightAngle = pow(dot(-rayDirection, normalize(fromLight)), 512.0);
      //   color.rgb += lightMasked * mix(lights[i].color, vec3(1), lightAngle) * pow(dot(-rayDirection, normalize(fromLight)), 512.0);
      // }

      // // Cartoon outline
      // // Requires trap be the distance even when the object is missed
      // // Doesn't detect edges not on the background.
      // float outlineStop = 0.00125;
      // vec3 outlineColor = vec3(0);
      // float outlineT = t.w;
      // outlineT = smoothstep(outlineStop, 0.125 * edge + outlineStop, outlineT);
      // outlineT = 1. - outlineT;
      // color = mix(color, vec4(outlineColor, 1), outlineT);

      // Radial Gradient
      // color = mix(vec4(vec3(0), 1.0), vec4(background, 1), saturate(pow((length(uv) - 0.25) * 1.6, 0.3)));

      // // Glow
      // float stepScaleAdjust = 0.03;
      // vec2 a = polarCoords(uv);

      // // t.z += 1.20 * snoise2(2123. * uv);

      // vec3 s = vec3(0);
      // a.x *= 3.;
      // t.z += 3.20 * fbmWarp(vec3(a, 0.5 * cos(-length(uv) + cosT)), s);
      // // t.z += length(s);
      // // t.z += cellular(3. * s);

      // float i = saturate(t.z / (stepScaleAdjust * float(maxSteps)));
      // // float i = 1. - saturate(pow(2.0 * t.w, 0.25));
      // vec3 glowColor = vec3(1, 0.9, 0);
      // // const float stopPoint = 0.5;
      // // i = smoothstep(stopPoint, stopPoint + edge, i);
      // i = pow(i, 1.25);
      // color = mix(color, vec4(glowColor, 1.0), i);

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

// Credit: IQ
// source: https://www.shadertoy.com/view/3l23RK
// c is the sin/cos of the angle. r is the radius
float sdPie( in vec2 p, in vec2 c, in float r )
{
  p.x = abs(p.x);
  float l = length(p) - r;
  float m = length(p - c*clamp(dot(p,c),0.0,r) );
  return max(l,m*sign(c.y*p.x-c.x*p.y));
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

float cLength2 (in vec2 q) {
  return dot(q,q);
}

vec2 cSquare (in vec2 q) {
  return vec2(
      q.x * q.x - q.y * q.y,  // real
      2. * q.x * q.y);        // complex
}

vec2 cCube (in vec2 q) {
  // q³ = (q.x² + 2. * q.x * q.y * i + - q.y²) * (q.x + q.y * i)
  // q³ = (  q.x³ + q.y * q.x² * i
  //        + 2. * q.x² * q.y * i + - 2. * q.x * q.y²
  //        + - q.x * q.y² - q.y³ * i)
  // q³ =  q.x³ + - 2. * q.x * q.y² + - q.x * q.y²
  //        + (q.y * q.x² + 2. * q.x² * q.y - q.y³) * i 
  return vec2(
      q.x * q.x * q.x - 3. * q.x * q.y * q.y,   // real
      3. * q.x * q.x * q.y - q.y * q.y * q.y);  // complex
}

float star (in vec2 q, in float r, in float n) {
  float d = maxDistance;
  float angle = TWO_PI / n;
  float c = pModPolar(q, n);
  q.y = abs(q.y);

  q *= rotMat2(0.5 * angle);

  float line = abs(q.x - r);
  d = min(d, line);
  return d;
}

float rotateStarAngle (in float r0, in float r1, in float n) {
  //  theta = top angle
  // b|\ a
  //  |/ c
  float b = r0;
  float c = r1;

  float theta = PI / n;
  float a = b * cos(theta) - sqrt(c * c - b * b * pow(sin(theta), 2.));

  //  theta = top angle
  // b|\ a
  //  |- <- d
  //  |/ c
  // sin(theta) = d / a
  float d = a * sin(theta);

  //  theta = top angle
  // b|\ a
  //  |- <- d
  //  |/ c
  //  phi = bottom angle
  float phi = asin(d / c);

  // hmm closer..

  // debug
  // phi = theta;
  // phi = a;
  // phi = d;

  // what. doubling it is close to perfect (might even be perfect). I dont know
  // why. definitely breaks when it isn't n = 8
  return phi;
}

vec2 factorCenter (in vec2 q, in float num, in float r) {
  pModPolar(q, num);

  // okay so r needs to be such that it creates a side with length 2 * r.
  //    /b __  
  // x / |\  x 
  //  /r | \ | 
  // //_\|_a\  
  // a  2r   
  //
  // θ = b / 2
  // sin(θ) = r / x
  // x = r / sin(θ)
  q.x -= r / sin(TWO_PI / num * 0.5);

  return q;
}

float factor (in vec2 q, in float num, in float r) {
  float d = maxDistance;

  q = factorCenter(q, num, r);

  float c = length(q) - r;
  d = min(d, c);

  return d;
}

vec3 factorColor (in float layerI, in float totalLayers, in vec2 uv) {
  vec3 color = vec3(0);

  float percent = layerI / totalLayers;

  // percent *= angle1C;
  // percent += angle2C;

  // color = 0.5 + 0.5 * cos(TWO_PI * (0.2 * vec3(uv, 0) + percent + vec3(0, 0.33, 0.67)));
  // color += 0.4 * (0.5 + 0.5 * cos(TWO_PI * (0.2 * vec3(uv, 0) + color + vec3(0, 0.33, 0.67))));
  // color *= 0.8;

  // Greyscale
  color = vec3(percent);

  color = pow(color, vec3(1.5));

  return color;
}

const int NUM_CIRCLES = 10;
// Okay now to figure out how to get the last angle to perfectly complete the border
// okay. i think i figured it out.
float getCircleBorderAngle (in float i, in float t) {
  float angleInc = TWO_PI / float(NUM_CIRCLES);

  return i * angleInc + angleInc * (0.1 * snoise2(2.91283 * vec2(i)) + 0.075 * cos(t + (0.12 + 0.07 * noise(vec2(i))) * PI * i));
}

float getCircleBorderRadius (in float fI, in float baseR, in float t) {
  float r0 = 0.0612062695 + 0.0045 * cos(t);

  float r = r0;
  for (float i = 1.; i < float(NUM_CIRCLES); i++) {
    if (i > fI) break;
    float prevI = float(mod(i - 1., float(NUM_CIRCLES)));

    vec2 center = vec2(baseR, 0) * rotMat2(getCircleBorderAngle(i, t));
    vec2 other = vec2(baseR, 0) * rotMat2(getCircleBorderAngle(prevI, t));

    r = length(center - other) - r;
  }

  return r;
}

float circleEdgeCircle (in vec2 q, in float localCosT) {
// #define debugCircleEdge 1

  float d = maxDistance;

  const float baseR = 0.25;

#ifdef debugCircleEdge
  // Circle the edge is based on
  float m = length(q) - baseR;
  m = abs(m);
  d = min(d, m);
#endif

  float angleInc = TWO_PI / float(NUM_CIRCLES);

  float accAngle = 0.;
  float prevR = 0.;
  vec2 prevCenter = vec2(0);

  float subs = maxDistance;

  for (int i = 0; i < NUM_CIRCLES - 1; i++) {
    float fI = float(i);
    float prevI = float(mod(fI - 1., float(NUM_CIRCLES)));

    float angle = getCircleBorderAngle(fI, localCosT);
    vec2 center = vec2(baseR, 0) * rotMat2(angle);
    vec2 otherCenter = vec2(baseR, 0) * rotMat2(getCircleBorderAngle(prevI ,localCosT));

    float r = getCircleBorderRadius(fI, baseR ,localCosT);

    vec2 localQ = q - center;

    // Add "fill"
    float fill = sdTriangle(q, vec2(0), center, otherCenter);
#ifdef debugCircleEdge
    fill = abs(fill);
#endif
    if (fI != 0.) {
      // d = min(d, fill);
    }

    float o = length(localQ) - r;
#ifdef debugCircleEdge
    o = abs(o);
#endif
    // if (mod(fI, 2.) == 0.) {
      d = min(d, o);
    // } else {
// #ifndef debugCircleEdge
    //   subs = min(subs, o);
// #else
    //   d = min(d, o);
// #endif
    // }

    prevCenter = center;
    prevR = r;
  }

  // -- Add last circle --
  float fI = float(NUM_CIRCLES - 1);
  vec2 firstCenter = vec2(baseR, 0) * rotMat2(getCircleBorderAngle(0., localCosT));
  float r0 = getCircleBorderRadius(0., baseR, localCosT);
  vec2 secondToLastCenter = prevCenter;
  float rSecondToLast = prevR;
  float remainderSpace = length(firstCenter - secondToLastCenter) - (r0 + rSecondToLast);
  vec2 middle = secondToLastCenter
    + normalize(firstCenter - secondToLastCenter) * (remainderSpace * 0.5 + rSecondToLast);
  vec2 center = normalize(middle) * baseR;

  float r = length(firstCenter - center) - r0;

  vec2 localQ = q - center;

  // Add "fill"
  float fill = sdTriangle(q, vec2(0), center, secondToLastCenter);
#ifdef debugCircleEdge
  fill = abs(fill);
#endif
  // d = min(d, fill);

  // To first center
  fill = sdTriangle(q, vec2(0), center, firstCenter);
#ifdef debugCircleEdge
  fill = abs(fill);
#endif
  // d = min(d, fill);

  float o = length(localQ) - r;
#ifdef debugCircleEdge
  o = abs(o);
#endif
  // if (mod(fI, 2.) == 0.) {
    d = min(d, o);
  // } else {
#ifndef debugCircleEdge
    // subs = min(subs, o);
#else
    // d = min(d, o);
#endif
    // d = max(d, -o);
  // }

#ifndef debugCircleEdge
  d = max(d, -subs);
#endif

  return d;
}

float uShape (in vec2 q, in float r, in float size) {
  // "U" shape
  r = r * 0.707107; // sqrt(2) * 0.5
  float o = length(q) - r;

  // Add straight section
  vec2 straightQ = q;
  straightQ -= vec2(0.25 * size);
  straightQ *= rotMat2(0.25 * PI);
  float straight = sdBox(straightQ, vec2(r));
  o = min(o, straight);

  return o;
}

vec3 truchetSpace (in vec2 q, in float size, in vec2 seed) {
  vec2 c = pMod2(q, vec2(size));

  vec2 truchetQ = q;
  if (h21(c + seed) < 0.5) truchetQ.x *= -1.; // Flip, don't rotate
  if (dot(truchetQ, vec2(1)) <= 0.) truchetQ *= -1.;

  // TODO maybe remove as this is for shapes at the corner specifically
  truchetQ -= vec2(0.5 * size);

  // Angle
  float angle = atan(truchetQ.x, truchetQ.y);
  angle += PI;
  angle *= 2. / PI;
  float checker = mod(dot(c, vec2(1)), 2.);
  angle *= checker * 2. - 1.;

  return vec3(truchetQ, angle);
}

float roseCircles (in vec2 q, in float d, in float r) {
  vec2 sQ = (q);
  sQ *= rotMat2(cosT);
  sQ.x -= 0.5 * r;
  float s = length(sQ) - r;
  d = mix(d, 1. - d, step(0., s));

  sQ.x += 1.0 * r;
  s = length(sQ) - r;
  d = mix(d, 1. - d, step(0., s));

  sQ.xy += r * vec2(-0.5, 0.5);
  s = length(sQ) - r;
  d = mix(d, 1. - d, step(0., s));

  sQ.y += r * -1.0;
  s = length(sQ) - r;
  d = mix(d, 1. - d, step(0., s));

  s = length(q) - 1.5 * r;
  d = mix(d, 1. - d, step(0., s));

  return d;
}

float uiBoxCorner (in vec2 q, in vec2 bigBoxR) {
  vec2 boxR = vec2(3.0e-4, 0.009);

  q.x *= -1.;
  q -= vec2(-1, 1) * bigBoxR;

  // Reflect diagonally
  if (q.y > -q.x) {
    q.xy = -q.yx;
  }
  // Move "down" so full height of box is shown
  q.y += boxR.y;

  return sdBox(q, boxR);
}

float uiBoxCorners (in vec2 q, in vec2 bigBoxR) {
  q = abs(q.xy);
  return uiBoxCorner(q, bigBoxR);
}

float stripedBox (in vec2 q, in vec2 r, in float stripeFreq, in float stripeThick, in float angle) {
  vec2 stripeR = vec2(stripeFreq * stripeThick, 0.5);
  vec2 stripeQ = q;
  stripeQ *= rotMat2(angle);
  vec2 rotStripeQ = stripeQ;
  pMod1(stripeQ.x, stripeFreq);

  float s = sdBox(stripeQ, stripeR);
  float stripeMask = sdBox(q, r);
  return max(s, stripeMask);
}

float stripedBox (in vec2 q, in vec2 r, in float stripeFreq, in float stripeThick) {
  return stripedBox(q, r, stripeFreq, stripeThick, 0.25 * PI);
}

float get2DSDF (in vec2 q) {
  // Convert from center 0 to center 0.5
  q += 0.5;
  float sdf2D = texture2D(sdf2DTexture, q).r;
  // Unpack from [0, 1] to [-1, 1]
  sdf2D -= 0.5;
  sdf2D *= 2.;

  return sdf2D;
}

float bunch (in vec2 q, in float r, in float t , float cOffset, float cX) {
  const float tNum = 1.;

  float ySize = 2. * r;
  q.y += ySize * cOffset;
  float c = floor((q.y + 0.5 * ySize) / ySize);
  c += cOffset;
  q.y -= tNum * ySize * expo(t);
  c = pMod1(q.y, ySize);

  float virtualC = c + tNum * t;

  const float waveAmount = 0.;
  q.y -= waveAmount * r * cos(PI * (8. + pow(virtualC, 1.25)) * q.x);
  q.y -= waveAmount * r * snoise2(vec2(18. * q.x + 0.273 * cX, virtualC)); // * (1. - smoothstep(0.99, 1., abs(q.x)));

  return sdBox(q, vec2(r, 0.2 * r));
}

float ripple (in vec2 q, in float t) {
  float d = maxDistance;
  const float thickness = 0.01;
  float r = -0.0125;

  r += 1. * 0.5 * (0.45 * expoIn(range(0., 0.4, t)) + 1.55 * quartOut(range(0.4, 1., t)));

  float b = abs(length(q) - r) - thickness;
  d = min(d, b);

  q *= rotMat2(0.5 * PI * t);

  q = polarCoords(q);
  q.x /= PI;
  q.y -= r;
  q.y += 6. * thickness;

  pMod1(q.x, 2.5 * thickness);

  // // float s = length(q) - thickness;
  // float s = sdBox(q, vec2(thickness));
  // d = min(d, s);

  return d;
}

#define dotNum 8.
vec2 dotPosition (in float i, in float t, in float r) {
  const vec2 start = vec2(0);
  vec2 q = start;

  float maxI = floor((dotNum + 1.) * t) + 1.;
  float withinStep = expo((dotNum + 1.) * mod(t, 1./(dotNum + 1.)));

  float isFirst = expo(range(0.1/dotNum, 1./dotNum, t));

  vec2 offset = mix(0., r, isFirst) * cos(TWO_PI * saturate(i / maxI) + vec2(0, 0.5 * PI));
  vec2 prevOffset = mix(0., r, isFirst) * cos(TWO_PI * saturate(i / (maxI - 1.)) + vec2(0, 0.5 * PI));

  return (1. - quint(range((dotNum - 1.) / dotNum, 1., t))) * mix(prevOffset, offset, withinStep);
}

vec2 mUv = vec2(0);
vec4 two_dimensional (in vec2 uv, in float generalT) {
  vec4 color = vec4(vec3(0), 1.);
  vec2 d = vec2(1, -1);

  vec2 q = uv;

  // Global Timing
  // generalT = angle1C;
  float t = mod(generalT + 0.0, 1.0);
  localCosT = TWO_PI * t;
  localT = t;

  float warpScale = 1.00;
  float warpFrequency = 1.;

  vec2 r = vec2(0.01);
  vec2 size = vec2(2.5) * vmax(r);

  // -- Warp --
  vec2 wQ = q.xy;

  float warpT = localCosT;

  // vec2 c = floor((wQ + size*0.5)/size);

  // // Odd row offset
  // wQ.x += 0.5 * size.x * mod(c.y, 2.);

  vec2 c = vec2(0);

  // Fake "Isometric" perspective
  wQ.y *= 1.35;
  wQ *= rotMat2(0.2 * PI);

  // wQ += 0.100000 * warpScale * cos( 3.0 * warpFrequency * componentShift(wQ) + warpT );
  // wQ += 0.050000 * warpScale * cos( 9.0 * warpFrequency * componentShift(wQ) + warpT );
  // wQ += 0.050000 * warpScale * snoise2(1. * warpFrequency * componentShift(wQ));
  // wQ += 0.025000 * warpScale * cos(15.0 * warpFrequency * componentShift(wQ) + warpT );

  // wQ *= rotMat2(-0.5 * PI);

  // wQ = polarCoords(wQ);
  // wQ.x /= PI;
  // wQ.y -= 0.2;

  c = floor((wQ + size*0.5)/size);
  wQ = opRepLim(wQ, vmax(size), vec2(11));
  // c = pMod2(wQ, vec2(0.50));
  // gC = c;
  // c.y += cIshShift;

  q = wQ;
  mUv = q;

  // // Adjust R per cell
  // r *= 0.7 - 0.4 * snoise2(1.7238 * c);

  // -- Cell T --
  float cellT = t;

  // Center out
  cellT -= 0.025 * length(c);

  // // Coordinate offset
  // // cellT -= 0.080 * c.y;
  // cellT -= 0.020 * c.x;

  // // Vmax offset
  // cellT -= 0.1 * vmax(vec2(vmin(c), dot(c, vec2(-1, 1))));
  // cellT -= 0.20 * vmax(abs(c));

  // // Dot product offset
  // float dC = dot(c, vec2(1, -1));
  // cellT -= dC * 0.075;

  // // Noise offset
  // cellT -= 0.175 * snoise2(1.2 * c);

  // Rectify
  cellT = mod(cellT, 1.);

  // cellT = triangleWave(cellT);
  // cellT = range(0.0, 1., cellT);

  // // Scale radius from -x to 1 where x is 0.05
  // r *= cellT + 0.05 * (-1. + cellT);

  // -- Local Space offsets ---
  // // Shift by random noise
  // q += 0.4 * vmax(r) * vec2(
  //     snoise2(c + vec2( 0.0100,-0.9000)),
  //     snoise2(c + vec2(-9.7000, 2.7780)));

  // --- Distance field(s) ---
  // // Texture SDF
  // float sdf2D = get2DSDF(q);
  // vec2 o = vec2(sdf2D, 0);
  // d = dMin(d, o);

  // q *= rotMat2(PI * cos(localCosT + dot(c, vec2(0.15))));

  // vec2 b = vec2(sdBox(q, vec2(r)), 0);
  // d = dMin(d, b);

  // vec2 b = vec2(length(q - r * vec2(1)) - 1.8 * vmax(r), 0);
  //
  r -= 1.2 * r * quart(0.5 + 0.5 * cos(TWO_PI * cellT));
  vec2 b = vec2(sdBox(q, 0.9 * r), 0);
  d = dMin(d, b);

  // vec2 b = vec2(neighborGrid(q, gSize).x, 0);
  // d = dMin(d, b);

  // --- Mask ---
  float mask = 1.;
  vec2 maskQ = wQ;

  // vec2 maskSize = vec2(boxIshR, 2. * evaporateR);
  // mask = sdBox(c - vec2(0, maskSize.y - maskSize.x), maskSize);

  // mask = length(maskQ) - 0.40;
  // mask = sdBox(maskQ, vec2(r));
  // mask = abs(vmax(abs(maskQ)) - 0.3) - 0.1;

  // // mask = max(mask, -sdBox(maskQ, vec2(0.05, 2.)));
  // mask = smoothstep(0., edge, mask);
  // mask = 1. - mask;
  // // mask = 0.05 + 0.95 * mask;

  // --- Output ---
  float n = d.x;

  // // Repeat
  // n = sin(0.25 * TWO_PI * n);

  // // Outline
  // n = abs(n) - mix(-9. * thickness, thickness, pow(boxT, 0.5));

  // // Cyan glow
  // color.rgb = 0.8 * vec3(0, 1.0, 0.4) * mix(0., 1., saturate(1. - 1.8 * saturate(pow(saturate(n + 0.00), 0.125))));

  // Hard Edge
  n = smoothstep(0., 1.0 * edge, n - 0.0);

  // Invert
  n = 1. - n;

  // // Solid
  // color.rgb = vec3(1);

  // B&W
  color.rgb = vec3(n);

  // // B&W Repeating
  // color.rgb = vec3(0.5 + 0.5 * cos(TWO_PI * n));

  // // Simple Cosine Palette
  // float cosineIndex = d.y;
  // cosineIndex *= 0.17283;
  // cosineIndex += t;
  // cosineIndex += dot(uv, vec2(1));
  // color.rgb = saturate(n) * (0.5 + 0.5 * cos(TWO_PI * (cosineIndex + vec3(0, 0.33, 0.67))));
  // color.a = 1.;

  // // Mix
  // color.rgb = mix(vec3(0., 0.05, 0.05), vec3(1, .95, .95), n);

  // // JS colors
  // color.rgb = mix(colors1, colors2, n);

  // TODO This is too messy. Move it into functions / modules

  // // Cosine Palette
  // vec3 dI = vec3(n);
  // // dI += 0.125 * fallOff;
  // dI += dot(uv, vec2(0.4));
  // dI += 0.2 * cos(localCosT + dot(uv, vec2(0.2, -0.4)));
  // dI *= 0.75;
  // // color.rgb = mix(color.rgb, n * (0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.33, 0.67)))), isMaterialSmooth(d.y, 1.));
  // color.rgb = 0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.33, 0.67)));

  // // Stripes
  // const float numStripes = 8.;
  // vec2 axis = vec2(1, 0) * rotMat2(0.5 * PI * n);
  // float line = dot(q, axis);
  // line = sin(TWO_PI * numStripes * line);
  // line -= 0.95;
  // line = smoothstep(0., 2. * edge, line);
  // color.rgb = vec3(line);

  // // radial stripes
  // float angle = atan(q.y, q.x);
  // angle += 6. * n;
  // float line = angle;
  // line = sin(TWO_PI * 30. * line);
  // line = smoothstep(0., edge, line);
  // color.rgb = vec3(line);
  // color.rgb = mix(vec3(0), color.rgb, step(edge, n));

  // // Grid spinners?
  // const float baseGridSize = 0.10;
  // float gridSize = baseGridSize;
  // gridSize += 0.3 * gridSize * cos(localCosT - length(q));

  // vec2 c = pMod2(q, vec2(gridSize));
  // q *= rotMat2(localCosT + 10. * n - 0.05 * length(c) + 0.75 * PI * cnoise2(c));
  // // float line = abs(q.y) - 0.015625 * baseGridSize;
  // float line = sdBox(q, vec2(0.015625, 0.3) * baseGridSize);
  // line = smoothstep(1.0 * edge, 0., line);
  // // line *= step(0., -sdBox(c * rotMat2(0.25 * PI), vec2(5)));
  // line *= step(0., -sdBox(c, vec2(8)));
  // color.rgb = vec3(line);

  // // Grid crosses
  // float gridSize = 0.0175;
  // q = uv;
  // vec2 c = pMod2(q, vec2(gridSize));
  // q *= rotMat2(-localCosT + 12. * n - 0.05 * length(c));
  // float line = min(abs(q.x), abs(q.y)) - 0.125 * 0.015625 * gridSize;
  // // line = max(line, sdBox(q, vec2(0.25 * gridSize)));
  // line = smoothstep(edge, 0., line);
  // color.rgb = vec3(line);

  // // Grid circles
  // float gridSize = 0.02;
  // q = uv;
  // vec2 c = pMod2(q, vec2(gridSize));
  // float line = length(q) - 0.125 * gridSize * (-0.2 + 1.2 * range(0.1, 1., 0.5 + 0.5 * cos(PI * 1.123 * snoise2(3.0237 * c) - localCosT)));
  // line = smoothstep(0.5 * edge, 0., line);
  // color.rgb = vec3(line);

  // // Tint
  // color.rgb *= vec3(1, 0.9, 0.9);

  // // Darken negative distances
  // color.rgb = mix(color.rgb, vec3(0), 0.2 * smoothstep(0., 3. * edge, -n));

  color.rgb *= saturate(mask);
  // color.rgb *= color.a; // Don't leak color channels at expense of edges loosing color

  // color.rgb *= 1.3;

  return color;
}

vec4 two_dimensional (in vec2 uv) {
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

const vec3 sunPos = vec3(0.01,  0.,-1.2);
vec3 sunColor (in vec3 q) {
  const float darker = 0.8;
  return vec3(1, darker, darker);
  // return mix(vec3(darker, darker, 1), vec3(1, darker, darker), 0.5 * length(q.xy));
}

// #pragma glslify: godRays = require(./effects/god-rays.glsl, softshadow=softshadow, noise=h21, sunPos=sunPos, sunColor=sunColor)
#pragma glslify: godRays = require(./effects/god-rays-poisson.glsl, shadow=simpleshadow, noise=h21, volumeNoise=vfbmWarp, sunPos=sunPos, sunColor=sunColor)

// renderSceneLayer()
// This returns a rendered scene layer. Similar to `sample()` it takes:
// - `ro` : Ray origin vector
// - `rd` : Ray direction / Eye vector
// - `uv` : UV coordinate for a 'screen space'
// - `time` : The current time for the layer w/ a range of [0, 1)
// and returns a rgba color value for that coordinate of the scene.
vec4 renderSceneLayer (in vec3 ro, in vec3 rd, in vec2 uv, in float time) {

// #define is2D 1
#ifdef is2D
  // 2D
  vec4 layer = two_dimensional(uv, time);

#else
  // 3D
  vec4 t = march(ro, rd, time);
  vec4 layer = shade(ro, rd, t, uv, time);

  // // -- 3D : Effects --
  // layer = godRays(ro, rd, t, uv, layer, time);

#endif

  return layer;
}

vec4 renderSceneLayer (in vec3 ro, in vec3 rd, in vec2 uv) {
  return renderSceneLayer(ro, rd, uv, norT);
}

#pragma glslify: outline = require(./effects/outline.glsl, map=two_dimensional, steps=3.)

// sample()
// This function gets a textile / sample given:
// - `ro` : Ray origin vector
// - `rd` : Ray direction / Eye vector
// - `uv` : UV coordinate for a 'screen space'
// and returns a rgba color value for that sample.
vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  vec4 color = vec4(0, 0, 0, 1);

  // -- Single layer --
  return renderSceneLayer(ro, rd, uv);

  // // -- Single layer : Outline --
  // float layerOutline = outline(uv, angle3C);
  // // Hard Edge
  // layerOutline = smoothstep(0., 0.40 * edge, layerOutline - angle2C);

  // return vec4(vec3(1. - layerOutline), 1);

  // // -- Echoed Layers --
  // const float echoSlices = 8.;
  // for (float i = 0.; i < echoSlices; i++) {
  //   vec4 layerColor = renderSceneLayer(ro, rd, uv, norT - 0.04 * i);

  //   // // Outlined version
  //   // float layerOutline = outline(uv, angle3C, norT - 0.0075 * i);
  //   // // Hard Edge
  //   // layerOutline = smoothstep(0., 0.20 * edge, layerOutline - angle2C);
  //   // vec4 layerColor = vec4(vec3(1. - layerOutline), 1);

  //   // Echo Dimming
  //   // layerColor *= (1. - pow(i / (echoSlices + 1.), 0.125));
  //   layerColor.a *= (1. - pow(i / (echoSlices + 1.), 0.067));

  //   // Blend mode
  //   // Additive
  //   color += vec4(vec3(layerColor.a), 1) * layerColor;

  //   // color.rgb = mix(color.rgb, layerColor.rgb, layerColor.a);

  //   // -- Offsets --
  //   // Incremental offset
  //   uv.y += 0.0080;

  //   // Initial Offset
  //   uv.y += i == 0. ? 0.0075 : 0.;

  //   // uv.y += 0.0125 * i * loopNoise(vec3(norT, 0.0000 + 2. * uv), 0.3, 0.7);
  //   // uv.y += 0.012 * i * abs(snoise3(vec3(uv.y, sin(TWO_PI * norT + vec2(0, 0.5 * PI)))));
  // }

  // color.a = saturate(color.a);
  // // color.rgb = mix(vec3(1), color.rgb, color.a);
  // color.rgb += pow(1. - color.a, 1.3) * vec3(0);
  // color.a = 1.;

  // return color;

  // -- Color delay --
  const float slices = 8.;
  const float delayLength = 0.180;

  for (float i = 0.; i < slices; i++) {
    vec3 layerColor = vec3(0.);

    // -- Apply Scene as Mask --
    float layerT = norT
      - delayLength * i / slices;
    layerColor = renderSceneLayer(ro, rd, uv, layerT).rgb;

    // -- Get Layer Color --
    vec3 layerTint = vec3(0.); // 0.5 + 0.5 * cos(TWO_PI * (i / slices + vec3(0, 0.33, 0.67)));
    // vec3 layerTint = vec3(
    //     saturate(mod(i + 0., 3.)),
    //     saturate(mod(i + 1., 3.)),
    //     saturate(mod(i + 2., 3.))
    // );

    // Cosine Palette
    vec3 dI = vec3(i / slices);
    dI += dot(uv, vec2(0.7));
    // dI += 0.125 * snoise2(vec2(2, 1) * mUv);

    // dI *= 0.6;

    dI += 0.1 * cos(cosT + dot(uv, vec2(-1, 1))); // Vary over time & diagonal space

    // layerTint = 1.00 * (vec3(0.5) + vec3(0.5) * cos(TWO_PI * (vec3(0.5, 1, 1) * dI + vec3(0., 0.2, 0.3))));
    layerTint = 1.0 * (0.5 + 0.5 * cos(TWO_PI * (vec3(1) * dI + vec3(0, 0.33, 0.67))));
    // layerTint += 0.8 * (0.5 + 0.5 * cos(TWO_PI * (layerTint + pow(dI, vec3(2.)) + vec3(0, 0.4, 0.67))));
    // layerTint *= mix(vec3(1.0, 0.6, 0.60), vec3(1), 0.3);

    // // Solid Layer color
    // layerTint = vec3(5.0);

    // CYM
    // layerTint = vec3(0);
    // layerTint += vec3(0, 1, 1) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 1.0 * dI.x + vec3(0, 0.33, 0.67))));
    // layerTint += vec3(1, 0, 1) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 1.2 * dI.y + vec3(0, 0.33, 0.67))));
    // layerTint += vec3(1, 1, 0) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 0.8 * dI.z + vec3(0, 0.33, 0.67))));
    // layerTint *= 0.65;
    // layerTint *= vec3(1.0, 0.6, 0.60);

    // // -- Layer Post Processing --
    // layerTint *= colors1; // Tint w/ color 1
    // layerTint *= 1.5;

    // Add black layer as first layer
    // layerTint *= step(0.5, i);

    // // Darken exponentially based on total layers
    // layerTint = pow(layerTint, vec3(4. + slices));

    // -- Apply Layer Tint --
    layerColor *= layerTint;

    // -- Blend / Aggregate --
    // if (i == 0.) {
    //   color = layerColor;
    // } else {
    //   color = overlay(color, layerColor);
    // }
    // color *= layerColor;

    // vec3 layerColorA = softLight2(color, layerColor);
    // vec3 layerColorB = color * layerColor;
    // layerColor = layerColorA + layerColorB;
    // layerColor *= 0.85;

    // layerColor = overlay(color, layerColor);
    // layerColor = screenBlend(color, layerColor);
    // color = mix(color, layerColor, 0.3);

    // Add
    color.rgb += layerColor;

    // // Multiply
    // color.rgb *= layerColor;

    // // Pseudo Multiply
    // float mask = layerColor.x;
    // color.rgb = mix(color.rgb, color.rgb * layerColor, mask);
  }

  color.rgb = pow(color.rgb, vec3(1.750));
  color.rgb /= slices;

  // // Final layer
  // color.rgb += 1.0 * renderSceneLayer(ro, rd, uv, norT).rgb;

  // // Color manipulation
  // color.rgb = 1. - color.rgb;

  return color;
}

void main() {
  // Update timing
  modT = mod(time, totalT);
  norT = modT / totalT;
  cosT = TWO_PI / totalT * modT;

  const float orthoZoom = 0.5;

  vec3 ro = cameraRo + cOffset;

  vec2 uv = fragCoord.xy;

// #define pixelated
#ifdef pixelated
  // Pixelate UVs
  const float pixelSize = 1.0 * 0.009375;
  vec2 innerUV = uv;
  pMod2(innerUV, vec2(pixelSize));
  uv = floor((uv + pixelSize * 0.5) / pixelSize) * pixelSize;
#endif

    float gRAngle = -TWO_PI * mod(time, totalT) / totalT;
    float gRc = cos(gRAngle);
    float gRs = sin(gRAngle);
    globalRot = mat3(
      gRc, 0.0, -gRs,
      0.0, 1.0,  0.0,
      gRs, 0.0,  gRc);

#ifdef DOF
    const float dofCoeficient = 0.0175;
#endif

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec2 R = resolution * 2.;

    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 ssRo = ro;
#ifdef ORTHO
            vec3 rd = vec3(0, 0, -1);
#else
            vec3 rd = getRayDirection(vec2(
                  float(x) / R.y + uv.x,
                  float(y) / R.y + uv.y),
                  projectionMatrix);
#endif

            // jitter
            // rd.xy += 0.1 * noise(vec2(x, y));
            rd = (vec4(rd, 1.) * cameraMatrix).xyz;
            rd = normalize(rd);

#ifdef DOF
            // DoF
            // source: shadertoy.con/view/WtSfWK
            vec3 fp = ssRo + rd * doFDistance;
            ssRo.xy += dofCoeficient * vec2(
                cnoise2(238. * (uv + rd.xy) + 100. * cos(cosT) + 123. + 2384. * vec2(x, y)),
                cnoise2(323. * (uv + rd.xy) + 100. * cos(cosT) + 2034.123 * vec2(x, y)));
            rd = normalize(fp - ssRo);
#endif

#ifdef ORTHO
            vec2 ndc = (gl_FragCoord.xy + 0.0 * vec2(x, y)) / resolution.xy * 2.0 - 1.0;
            ndc *= orthoZoom;
            float w = 2.0;
            float h = w / (resolution.x / resolution.y);
            ssRo += (vec4(
                  ndc * vec2(w * 0.5, h * 0.5),
                  0, 1) * cameraMatrix).xyz;
#endif
            color += saturate(sample(ssRo, rd, uv));
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
#ifdef ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = getRayDirection(uv, projectionMatrix);
#endif

    // jitter
    // rd.xy += 0.001 * noise(2183. * uv);
    rd = (vec4(rd, 1.) * cameraMatrix).xyz;
    rd = normalize(rd);

#ifdef DOF
    // DoF
    // source: shadertoy.con/view/WtSfWK
    vec3 fp = ro + rd * doFDistance;
    ro.xy += dofCoeficient * vec2(
        cnoise2( 938.127 * (uv + rd.xy) + 100. * cos(cosT) + 123.),
        cnoise2(1323.379 * (uv + rd.xy) + 100. * cos(cosT) + 20034.123));
    rd = normalize(fp - ro);
#endif

#ifdef ORTHO
    vec2 ndc = gl_FragCoord.xy / resolution.xy * 2.0 - 1.0;
    ndc *= orthoZoom;
    float w = 2.0;
    float h = w / (resolution.x / resolution.y);
    ro += (vec4(
          ndc * vec2(w * 0.5, h * 0.5),
          0, 1) * cameraMatrix).xyz;
#endif
    gl_FragColor = saturate(sample(ro, rd, uv));
    #endif

#ifdef pixelated
    // Pixelated shading
    gl_FragColor.rgb *= 0.80 + 0.20 * (dot(innerUV, vec2(0.2, 1)) / pixelSize);
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
