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
#define maxSteps 256
#define maxDistance 10.0
#define fogMaxDistance 10.0

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
const float thickness = 0.01;

// Dispersion parameters
float n1 = 1.;
float n2 = 1.282;
const float amount = 0.05;

// Dof
float doFDistance = angle1C;

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
#pragma glslify: quartIn = require(glsl-easings/quadratic-in)
#pragma glslify: quartOut = require(glsl-easings/quadratic-out)
#pragma glslify: quint = require(glsl-easings/quintic-in-out)
#pragma glslify: quintIn = require(glsl-easings/quintic-in)
#pragma glslify: quintOut = require(glsl-easings/quintic-out)
// #pragma glslify: elasticInOut = require(glsl-easings/elastic-in-out)
#pragma glslify: elasticOut = require(glsl-easings/elastic-out)
// #pragma glslify: elasticIn = require(glsl-easings/elastic-in)

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
float sigmoid ( in float x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

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

vec3 sphericalCoords (in vec3 q) {
  float a = atan(q.z, q.x);
  float arc = atan(q.y, length(q.xz));
  return vec3(
      a,
      arc,
      length(q));
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

  float gapAmount = 0.0220 * (0.95 * mod(133.2830 * i, 1.0) + 0.05);
  float gap = gapAmount * quart(localT);

  float angle = snoise2(vec2(3.157143 * i) + 1.482349) * PI;

  float start = 0.15 * snoise2(vec2(10.8123 * i));

  return vec3(angle, gap, start);
}

const vec2 gSize = vec2(0.020);
float microGrid ( in vec2 q ) {
  vec2 cMini = pMod2(q, vec2(gSize * 0.10));

  return mod(dot(cMini, vec2(1)), 2.);
}

float localCosT = cosT;
float localT = norT;
float second = maxDistance;
vec2 shape (in vec2 q, in vec2 c) {
  vec2 d = vec2(maxDistance, -1.);
  // cool now I have material support

  vec2 uv = q;

  // float dC = vmax(abs(c));
  float dC = dot(c, vec2(-1, 2));

  float odd = mod(dC, 2.);
  float even = 1. - odd;

  const float warpScale = 0.;
  float size = gSize.y;

  // // Assume [0,1] range per dimension
  // vec2 bigSize = vec2(2);
  // vec2 bigC = floor(abs(c) / bigSize);
  // vec2 miniC = mod(c, bigSize);

  // Create a copy so there is no cross talk in neighborGrid
  float locallocalT = localT;
  // locallocalT = angle1C;
  // locallocalT -= 0.05 * length(c);
  locallocalT -= 0.0075 * dC;
  // locallocalT += 0.05 * odd;
  // NOTE Flip time offset if there are gaps
  // Might fix some of the gaps caused by the time offset
  // A hack but getting closer to a general solution
  // const float clip = 0.0125;
  // locallocalT = clamp(locallocalT, clip, 1. - clip);

  float t = mod(locallocalT, 1.);
  t = expo(t);
  float localCosT = TWO_PI * t;

  // Local C that transitions from one cell to another
  float shift = 8.;
  vec2 shiftDir = vec2(0, -1);

  vec2 localC = mix(c, c + shift * shiftDir, t);

  // // Vanilla cell coordinate
  // vec2 localC = c;

  float r = 0.025 * size;

  // // Make grid look like random placement
  // float nT = 0.5 + 0.5 * sin(localCosT); // 0.5; // triangleWave(t);
  // q += 0.25 * size * mix(
  //     vec2(1, -1) * snoise2(0.417 * localC + 23.17123),
  //     vec2(1) * snoise2(0.123 * localC),
  //     nT);

  // float side = step(abs(c.y), abs(c.x));
  // q.x += sign(c.x) * side * size * (0.5 + 0.5 * cos(localCosT));

  q -= shiftDir * shift * size * t;

  // // Rotate randomly
  // q *= rotMat2(1.0 * PI * snoise2(0.263 * localC));

  float internalD = length(q);
  // float internalD = abs(q.y);
  // internalD = max(internalD, abs(q.x) - 0.3 * size);
  // internalD = min(internalD, abs(q.x));
  // internalD = max(internalD, sdBox(q, vec2(0.3 * size, 0.1 * size)));

  // float internalD = abs(dot(q, vec2(-1, 1)));
  // internalD = max(internalD, sdBox(q, vec2(0.5 * size)));
  // float internalD = vmax(abs(q));
  // float internalD = dot(abs(q), vec2(1));
  // float internalD = sdBox(q, vec2(r));
  // vec2 absQ = abs(q);
  // float internalD = min(absQ.x, absQ.y);
  // float crossMask = sdBox(q, vec2(0.35 * size));
  // internalD = max(internalD, crossMask);

  // vec2 o = vec2(internalD, 0.);
  vec2 o = vec2(internalD - r, 0.);
  // float o = microGrid(q);
  d = dMin(d, o);

  // // Outline
  // const float adjustment = 0.0;
  // d = abs(d - adjustment) - r * 0.1000;

  // Mask
  // d = mix(d, maxDistance, step(0., dot(abs(c), vec2(1)) - 12.));
  // d = mix(d, maxDistance, step(0., vmax(abs(c)) - 12.));
  // d = mix(d, maxDistance, step(0., sdBox(c, vec2(14))));
  // d = mix(d, maxDistance, step(0., abs(length(c) - 4.) - 2.));
  // d = mix(d, maxDistance, step(0., length(c) - 15.));
  // // Convert circle into torus
  // d = mix(d, maxDistance, step(0., abs(length(c) - 13.) - 7.));

  return d;
}

// NOTE Don't fully understand whether this is correct or how it works
vec2 circleInversion (in vec2 q) {
  return q / dot(q, q);
  // q is point in the circle & a is the inverse
  // dot(q, a) = r * r
  // q.x * a.x + q.y * a.y = r * r // i don't know what the invert of a dot product is...
}

#pragma glslify: neighborGrid = require(./modulo/neighbor-grid, map=shape, maxDistance=maxDistance, numberOfNeighbors=8.)

float baseR = 0.4;
float thingy (in vec2 q, in float t) {
  float d = maxDistance;

  q *= 0.7;

  vec2 uv = q;

  localCosT = TWO_PI * t;
  localT = t;

  float thickness = 0.007;

  float r = 0.125;

  // Something is happening that looks circular
  float rollingScale = 1.;
  for (float i = 0.; i < 7.; i++) {
    q = circleInversion(q * 3.0) * 0.333;
    q = abs(q);
    // foldNd(q, vec2(angle1C, angle2C));
    q += offset.xy;
    q *= rotMat2(angle2C + 0.5 * localCosT);
    q *= scale;

    rollingScale *= scale;
  }

  // float o = dot(q, vec2(20)) * sdBox(q, vec2(r));
  float o = length(q) - r;
  o /= rollingScale;
  d = min(d, o);

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
  float longR = 2.5 * r.z;
  float crop = sdBox(q - vec3(0, 0, longR), vec3(innerR, longR));
  b = max(b, -crop);

  return b;
}

float gR = 0.6;
vec3 map (in vec3 p, in float dT, in float universe) {
  vec3 d = vec3(maxDistance, 0, 0);
  vec2 minD = vec2(1e19, 0);

  float t = mod(dT, 1.);
  float localCosT = TWO_PI * t;
  float r = gR;
  float size = 0.125;

  p *= globalRot;
  // p *= rotationMatrix(vec3(0, 1, 0), 0.15 * PI * cos(localCosT + 0.5 * p.x));

  vec3 q = p;

  float warpScale = 1.25;
  float warpFrequency = 1.0;
  float rollingScale = 1.;

  // Warp
  vec3 wQ = q.xyz;

  wQ += warpScale * 0.100000 * cos( 2. * wQ.xzx * warpFrequency + localCosT );
  wQ += warpScale * 0.050000 * cos( 7. * wQ.yzx * warpFrequency + localCosT );
  wQ.xzy = twist(wQ.xyz, (2.0 + 1.0 * cos(8. * wQ.y - localCosT)) * wQ.y / scale);
  wQ += warpScale * 0.025000 * cos(13. * wQ.xzx * warpFrequency + localCosT );
  wQ += warpScale * 0.012500 * cos(19. * wQ.yzx * warpFrequency + localCosT );
  wQ += warpScale * 0.006250 * cos(23. * wQ.yzx * warpFrequency + localCosT );

  // Commit warp
  q = wQ.xyz;

  mPos = q;

  // vec3 b = vec3(sdBox(q, vec3(r)), 0, 0);
  vec3 b = vec3(icosahedral(q, 52., r), 0, 0);
  b.x /= rollingScale;
  d = dMin(d, b);


  d.x *= 0.6;

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
    trap = min(trap, d.z);
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
  dI += 0.2 * pow(dNR, 3.);

  // dI += 0.25 * sin(TWO_PI * rd.x);
  dI *= angle1C;
  dI += angle2C;

  // dI += gPos;

  dI *= 0.3;

  // -- Colors --
  color = 0.5 + 0.5 * cos( TWO_PI * ( dI + vec3(0, 0.33, 0.67) ) );
  // color = mix(#FF0000, #00FFFF, 0.5 + 0.5 * sin(TWO_PI * (dI)));

  // color += 0.4 * (0.5 + 0.5 * cos( TWO_PI * ( color + dI + vec3(0, 0.1, 0.3) ) ));

  // color = vec3(n);

  color *= spread;

  // color = getBackground(rd.xy, 0.);

  // // Identity scene color
  // color = vec3(1);

  color *= 1.3;

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
  vec3 color = vec3(1);

  float dNR = dot(nor, -rd);
  vec3 dI = vec3(dNR);

  dI *= angle1C;
  dI += angle2C;

  color = 0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.3333, 0.67)));
  // color += 0.5 + 0.5 * cos(TWO_PI * (color + d * trap + vec3(0, 0.3333, 0.67)));
  color *= 1.3;

  color = mix(color, vec3(1), 0.5);

  return color;


  // return 2. * vec3(dNR);

  // dI += 0.5 * pow(dNR, 3.);
  // dI += 0.5 * snoise3(mPos);

  vec3 s = vec3(0);
  // dI += 0.5 * fbmWarp(18. * mPos, s);

  color = 0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.3333, 0.67)));
  color += 0.4 * (0.5 + 0.5 * cos(TWO_PI * (dI + vec3(0, 0.2, 0.4))));

  color *= 0.5;

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

    lights[0] = light(vec3(-0.2, 0.7, 1.0), 0.7 * #CCEEFF, 1.0);
    lights[1] = light(vec3(-0.5, 0.25, 1.0), #FFEEFF, 1.0);
    lights[2] = light(vec3(0.3, 0.3, 0.9), #FFFFEE, 1.0);

    const float universe = 0.;
    background = getBackground(uv, universe);

    float m = step(0., sin(TWO_PI * (0.25 * fragCoord.x + generalT)));

    float backgroundMask = 1.;
    // Allow anything in top right corner
    // backgroundMask = max(backgroundMask, smoothstep(0., edge, dot(uv, vec2(1)) + 0.05));

    if (t.x>0. && backgroundMask > 0.) {
      vec3 color = vec3(0.0);

      // Normals
      vec3 nor = getNormal2(pos, 0.001 * t.x, generalT);
      float bumpsScale = 1.8;
      float bumpIntensity = 0.105;
      nor += bumpIntensity * vec3(
          cnoise3(bumpsScale * 490.0 * mPos),
          cnoise3(bumpsScale * 670.0 * mPos + 234.634),
          cnoise3(bumpsScale * 310.0 * mPos + 23.4634));
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

      float freCo = 1.0;
      float specCo = 0.7;

      float specAll = 0.0;

      vec3 directLighting = vec3(0);
      for (int i = 0; i < NUM_OF_LIGHTS; i++) {
        vec3 lightPos = lights[i].position;
        // lightPos *= globalLRot; // Apply rotation
        vec3 nLightPos = normalize(lightPos);

        float diffMin = 0.0;
        float dif = max(diffMin, diffuse(nor, nLightPos));

        float spec = pow(clamp( dot(ref, nLightPos), 0., 1. ), 96.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        float shadowMin = 0.0;
        float sha = max(shadowMin, softshadow(pos, nLightPos, 0.01, 2.00, generalT));
        dif *= sha;

        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre; // Commit Fresnel
        specAll += specCo * spec;

        // // Ambient
        // lin += 0.100 * amb * diffuseColor;
        // dif += 0.100 * amb;

        float distIntensity = 1.; // lights[i].intensity / pow(length(lightPos - gPos), 1.0);
        distIntensity = saturate(distIntensity);
        color +=
          (dif * distIntensity) * lights[i].color * diffuseColor
          + distIntensity * mix(lights[i].color, vec3(1), 0.0) * lin * mix(diffuseColor, vec3(1), 1.0);

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
      // reflectColor += 0.10 * diffuseColor * reflection(pos, reflectionRd, generalT);
      // color += reflectColor;

      // vec3 refractColor = vec3(0);
      // vec3 refractionRd = refract(rayDirection, nor, 1.5);
      // refractColor += 0.10 * textures(refractionRd);
      // color += refractColor;

#ifndef NO_MATERIALS

      // vec3 dispersionColor = dispersionStep1(nor, normalize(rayDirection), n2, n1);
      vec3 dispersionColor = dispersion(nor, rayDirection, n2, n1);

      // float dispersionI = 1.0 * pow(1. - 1.0 * dot(nor, -rayDirection), 1.00);
      float dispersionI = 1.0;
      dispersionColor *= dispersionI;

      // dispersionColor.r = pow(dispersionColor.r, 0.7);

      color += saturate(dispersionColor);
      // color = saturate(dispersionColor);

#endif

      // // Fog
      // float d = max(0.0, t.x);
      // color = mix(background, color, saturate(pow(clamp(fogMaxDistance - d, 0., fogMaxDistance), 2.) / fogMaxDistance));
      // color *= saturate(exp(-d * 0.05));
      // // color = mix(background, color, saturate(exp(-d * 0.05)));

      // color += directLighting * exp(-d * 0.0005);

      // Inner Glow
      // color += 0.5 * innerGlow(5.0 * t.w);

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

      // // Cartoon outline
      // // Requires trap be the distance even when the object is missed
      // // Doesn't detect edges not on the background.
      // float outlineStop = 0.0005;
      // vec3 outlineColor = vec3(0);
      // float outlineT = t.w;
      // outlineT = smoothstep(outlineStop, 0.5 * edge + outlineStop, outlineT);
      // outlineT = 1. - outlineT;
      // color = mix(color, vec4(outlineColor, 1), outlineT);

      // Radial Gradient
      // color = mix(vec4(vec3(0), 1.0), vec4(background, 1), saturate(pow((length(uv) - 0.25) * 1.6, 0.3)));

      // // Glow
      // float stepScaleAdjust = 1.0;
      // float i = saturate(t.z / (stepScaleAdjust * float(maxSteps)));
      // vec3 glowColor = vec3(1);
      // const float stopPoint = 0.5;
      // // i = smoothstep(stopPoint, stopPoint + edge, i);
      // // i = pow(i, 0.90);
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

#pragma glslify: subdivide = require(./modulo/subdivide.glsl, vmin=vmin, noise=h21)

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

vec2 mUv = vec2(0);
vec3 two_dimensional (in vec2 uv, in float generalT) {
  vec3 color = vec3(0);
  vec2 d = vec2(maxDistance, -1);

  vec2 q = uv;

  // Global Timing
  // generalT = angle1C;
  float t = mod(generalT + 0.0, 1.0);
  localCosT = TWO_PI * t;
  localT = t;

  const float warpScale = 1.2;
  const vec2 size = gSize;
  float r = 0.25;
  vec2 seed = vec2(angle2C);

  vec2 wQ = q.xy;
  wQ *= rotMat2(-0.20 * PI);

  wQ += warpScale * 0.10000 * cos( 3. * vec2( 1, 1) * wQ.yx + 0. * localCosT + 0.1237);
  wQ += warpScale * 0.05000 * cos( 9. * vec2(-1, 1) * wQ.yx + 1. * localCosT + 1.937);
  wQ *= rotMat2(0.7 * PI + 0.5 * length(wQ) - 0.0125 * PI * cos(localCosT - 7.2 * length(wQ)));
  wQ += warpScale * 0.02500 * cos(16. * vec2( 1,-1) * wQ.yx + 1. * localCosT );
  wQ += warpScale * 0.01250 * cos(23. * vec2( 1, 1) * wQ.yx + 1. * localCosT );

  q = wQ;
  mUv = q;

  float x = dot(q, vec2(1));
  x = sin(TWO_PI * 30. * x);
  vec2 b = vec2(x, 0);
  d = dMin(d, b);

  // float bigMaskR = 0.35;
  // vec2 s = vec2(length(uv) - bigMaskR, 1.);
  // d = dMin(d, s);

  // float mask = maxDistance;
  // mask = smoothstep(0., 0.5 * edge, mask);
  // mask = 1. - mask;

  float n = d.x;

  // Hard Edge
  n = smoothstep(0., 1. * edge, n + 0.975);

  // Invert
  n = 1. - n;

  // // Solid
  // color = vec3(1);

  // B&W
  color = vec3(n);

  // // Mix
  // color = mix(vec3(0., 0.05, 0.05), vec3(1, .95, .95), n);

  // // JS colors
  // color = mix(colors1, colors2, n);

  // // Cosine Palette
  // vec3 dI = vec3(n);
  // // dI += 0.1238 * d.y;
  // color = 0.55 + 0.45 * cos(TWO_PI * (dI + vec3(0, 0.33, 0.67)));

  // // Stripes
  // const float numStripes = 60.;
  // vec2 axis = vec2(1, 0) * rotMat2(TWO_PI * n);
  // float line = dot(q, axis);
  // line = sin(TWO_PI * numStripes * line);
  // line = smoothstep(0., 2. * edge, line);
  // color = vec3(line);

  // // radial stripes
  // float angle = atan(q.y, q.x);
  // angle += 6. * n;
  // float line = angle;
  // line = sin(TWO_PI * 20. * line);
  // line = smoothstep(0., edge, line);
  // color = vec3(line);
  // color = mix(vec3(0), color, step(edge, n));

  // // Grid spinners?
  // float gridSize = 0.0175;
  // vec2 c = pMod2(q, vec2(gridSize));
  // q *= rotMat2(localCosT + 12. * n - 0.05 * length(c));
  // float line = abs(q.y) - 0.015625 * gridSize;
  // line = smoothstep(edge, 0., line);
  // color = vec3(line);

  // // Grid crosses
  // float gridSize = 0.0175;
  // q = uv;
  // vec2 c = pMod2(q, vec2(gridSize));
  // q *= rotMat2(-localCosT + 12. * n - 0.05 * length(c));
  // float line = min(abs(q.x), abs(q.y)) - 0.125 * 0.015625 * gridSize;
  // // line = max(line, sdBox(q, vec2(0.25 * gridSize)));
  // line = smoothstep(edge, 0., line);
  // color = vec3(line);

  // // Grid circles
  // float gridSize = 0.02;
  // q = uv;
  // vec2 c = pMod2(q, vec2(gridSize));
  // float line = length(q) - 0.125 * gridSize * (-0.2 + 1.2 * range(0.1, 1., 0.5 + 0.5 * cos(PI * 1.123 * snoise2(3.0237 * c) - localCosT)));
  // line = smoothstep(0.5 * edge, 0., line);
  // color = vec3(line);

  // // Tint
  // color *= vec3(1, 0.9, 0.9);

  // color *= mask;

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
  // return vec4(two_dimensional(uv, norT), 1);

  // vec3 color = vec3(0);

  // const int slices = 50;
  // for (int i = 0; i < slices; i++) {
  //   float fI = float(i);
  //   vec3 layerColor = vec3(0.); // 0.5 + 0.5 * cos(TWO_PI * (fI / float(slices) + vec3(0, 0.33, 0.67)));
  //   // vec3 layerColor = vec3(
  //   //     saturate(mod(fI + 0., 3.)),
  //   //     saturate(mod(fI + 1., 3.)),
  //   //     saturate(mod(fI + 2., 3.))
  //   // );

  //   vec3 dI = vec3(fI / float(slices));
  //   dI += 0.7 * dot(uv, vec2(1));
  //   // dI += 0.5 * snoise2(vec2(2, 1) * mUv);

  //   // layerColor = 1.00 * (vec3(0.5) + vec3(0.5) * cos(TWO_PI * (vec3(0.5, 1, 1) * dI + vec3(0., 0.2, 0.3))));
  //   layerColor = 1.0 * (0.5 + 0.5 * cos(TWO_PI * (vec3(1, 1, 1.5) * dI + vec3(0, 0.33, 0.67))));
  //   // layerColor += 0.8 * (0.5 + 0.5 * cos(TWO_PI * (layerColor + pow(dI, vec3(2.)) + vec3(0, 0.4, 0.67))));
  //   // layerColor *= mix(vec3(1.0, 0.6, 0.60), vec3(1), 0.3);
  //   layerColor *= colors1;
  //   layerColor *= 1.3;
  //   // layerColor = vec3(5.0);

  //   // CYM
  //   // layerColor = vec3(0);
  //   // layerColor += vec3(0, 1, 1) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 1.0 * dI.x + vec3(0, 0.33, 0.67))));
  //   // layerColor += vec3(1, 0, 1) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 1.2 * dI.y + vec3(0, 0.33, 0.67))));
  //   // layerColor += vec3(1, 1, 0) * 1.0 * (0.5 + 0.5 * cos(TWO_PI * (angle2C + 0.8 * dI.z + vec3(0, 0.33, 0.67))));
  //   // layerColor *= 0.65;
  //   // layerColor *= vec3(1.0, 0.6, 0.60);

  //   // layerColor *= 0.6;

  //   // Add black layer as first layer
  //   // layerColor *= step(0.5, fI);

  //   // layerColor = pow(layerColor, vec3(4 + slices));

  //   const float maxDelayLength = 0.05;
  //   float layerT = norT
  //     + maxDelayLength * (1.00 + 0.0 * sin(cosT + length(uv))) * fI / float(slices);
  //   float mask = two_dimensional(uv, layerT).x;
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

  //   // // Multiply
  //   // color *= layerColor;

  //   // // Pseudo Multiply
  //   // color = mix(color, color * layerColor, mask);
  // }

  // color = pow(color, vec3(1.50));
  // color /= float(slices);

  // // // Final layer
  // // color.rgb += 0.3 * two_dimensional(uv, 0.);

  // // // Color manipulation
  // // color.rgb = 1. - color.rgb;

  // return vec4(color, 1.);

  float time = norT;
  vec4 t = march(ro, rd, time);
  vec4 layer = shade(ro, rd, t, uv, time);
  return layer;
}

void main() {
  // Update timing
  modT = mod(time, totalT);
  norT = modT / totalT;
  cosT = TWO_PI / totalT * modT;

    const float orthoZoom = 0.5;

    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;

    float gRAngle = -TWO_PI * mod(time, totalT) / totalT;
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
            ssRo.xy += 0.0015 * vec2(
                cnoise2(238. * uv + 123. + 2384. * vec2(x, y)),
                cnoise2(323. * uv + 2034.123 * vec2(x, y)));
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
    ro.xy += 0.0015 * vec2(
        cnoise2(238. * uv + 123.),
        cnoise2(323. * uv + 20034.123));
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
