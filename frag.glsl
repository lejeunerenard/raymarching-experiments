#extension GL_OES_standard_derivatives : enable
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)
#define saturate(x) clamp(x, 0.0, 1.0)

// #define debugMapCalls
// #define debugMapMaxed
#define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform vec3 amberColor;
uniform bool BLOOM;
uniform vec3 cOffset;
uniform vec3 cameraRo;
uniform mat4 cameraMatrix;
uniform mat4 orientation;
uniform mat4 projectionMatrix;
uniform sampler2D tMatCap;
uniform sampler2D audioTexture;

uniform vec3 objectPos;
uniform float objectR;

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 1024
#define maxDistance 50.0

#define slowTime time * .2

vec3 lightPos = normalize(vec3(1., .75, 0.));
vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 gRd = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: pnoise3 = require(glsl-noise/periodic/3d)
#pragma glslify: vmax = require(./hg_sdf/vmax)
#pragma glslify: analyse = require(gl-audio-analyser)

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
  f += 0.062500*noise( p ); p = p*2.025;

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
        vfbm6(p + vec3(0.0, 0.0, 0.0)),
        vfbm6(p + vec3(3.2, 34.5, .234)),
        vfbm6(p + vec3(7.0, 2.9, -2.42)));

  s = vec3(
        vfbm6(p + scale * q + vec3(23.9, 234.0, -193.0)),
        vfbm6(p + scale * q + vec3(3.2, 852.0, 23.42)),
        vfbm6(p + scale * q + vec3(7.0, -232.0, -2.42)));

  return vfbm6(p + scale * s);
}
float vfbmWarp (vec3 p, out vec3 q) {
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

#pragma glslify: triprism = require(./model/tri-prism)

float triPrismGuide ( in vec3 p, in vec2 h ) {
  float outD = 1000.0;

  float axis = sdCapsule(p, vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, 1.0), 0.0125);
  outD = min(outD, axis);

  float point1 = sdCapsule(p, vec3(0.0, 1.0, -1.0), vec3(0.0, 1.0, 1.0), 0.0125);
  outD = min(outD, point1);

  float p2Angle = TWO_PI * 0.3333;
  float point2 = sdCapsule(p, vec3(sin(p2Angle), cos(p2Angle), -1.0), vec3(sin(p2Angle), cos(p2Angle), 1.0), 0.0125);
  outD = min(outD, point2);

  float p3Angle = -TWO_PI * 0.3333;
  float point3 = sdCapsule(p, vec3(sin(p3Angle), cos(p3Angle), -1.0), vec3(sin(p3Angle), cos(p3Angle), 1.0), 0.0125);
  outD = min(outD, point3);

  return outD;
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
#pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
#pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM, Iterations=14)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
// #pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

#define octaPreFold 2
mat4 octaM = mat4(
scale, 0., 0., 0.,
0., scale, 0., 0.,
0., 0., scale, 0.,
1., 1., 1., 1.) * mat4(
1., 0., 0.1, .0,
0., 1., 0., 0.,
0., 0.2, 1., 0.,
0., 0., 0., 1.);
#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=octaPreFold, kifsM=kifsM, trapCalc=trapCalc)
// 
// #pragma glslify: fold = require(./folds)
#pragma glslify: foldNd = require(./foldNd)
// #pragma glslify: twist = require(./twist)

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

float gRAngle = TWO_PI * 0.05 * time + PI * 0.75;
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
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
// #pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
// #pragma glslify: ease = require(glsl-easings/bounce-in)
#pragma glslify: voronoi = require(./voronoi)
#pragma glslify: band = require(./band-filter)
// #pragma glslify: tetrahedron = require(./model/tetrahedron)

// Logistic function
float sigmoid ( in float x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

// Return value is (distance, material, orbit trap)
vec3 map (in vec3 p) {
  vec3 d;

  const float period = 20.0;
  const float transitionTime = 2.0;
  float modTime = mod(time, period);
  float mixT = saturate((modTime - transitionTime) / (period - transitionTime));

  p *= globalRot;
  vec4 q = vec4(p, 1.0);

  vec4 q11 = 0.2500 * cos(3.05 * q.yzwx + noise(q.xyz) + modTime);
  vec4 q12 = 0.2500 * cos(3.05 * q.yzwx + noise(q.xyz) + modTime - period);
  q += mix(q11, q12, mixT);

  q += 0.1250 * cos(9.1 * q.yzwx);
  q += 0.0625 * cos(27.3 * q.yzwx);

  q.z *= 2.0;

  float minD = 0.0;
  // q.xyz = octahedronFold(q.xyz, minD);

  d = vec3(sdBox(q.xyz, vec3(1)), 1.0, 0.0);
  d.x *= 0.1;

  vec3 sqrP = p;
  sqrP.xyz = sqrP.xzy;
  sqrP *= rotationMatrix(vec3(0, 0, 1), PI * cos(PI * 0.5 * slowTime));
  sqrP *= 0.75;
  sqrP.y *= 2.0;
  vec3 square = vec3(sdTorus88(sqrP, vec2(1.5, 0.0625)), 2.0, 0.0);
  square.x *= 0.25;

  d = dMin(square, d);

  return d;
}

vec4 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t = 0.00001;
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
  vec2 e = vec2(1.,0.) * .015 * eps;
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
  color = mix(#ef78FF, color, 1. - l * .0625);
  // Dark
  color = mix(#043210, color, clamp(exp(l) * .325, 0., 1.));
}

#pragma glslify: checker = require(glsl-checker)
#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: rgb2hsv = require(./rgb2hsv.glsl)
#pragma glslify: debugColor = require(./debug-color-clip)

const float n1 = 1.0;
const float n2 = 1.31;
const float amount = 0.1;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  // float n = iqFBM(8.0 * rd + 2305.0);
  // float v = smoothstep(0.0, 1.0, n);

  float n = sinoise3(5.0 * rd);
  float v = smoothstep(-1.0, 1.0, n);

  // float n = dot(cos(rd), vec3(0.125));
  // float v = saturate(n);

  color = vec3(v);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd, in float ior) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

vec3 seaGradient (in float t) {
  t = saturate(t);
  vec3 color = mix(#2A3B82, #6FFFF9, t);
  return color;
}
vec3 amberGradient (in float t) {
  t = saturate(t);
  vec3 color = mix(#3D1C0B, amberColor, smoothstep(0.3, 0.75, t));
  return mix(color, pow(#ED4F2C, vec3(2.2)), smoothstep(0.75, 1.0, t));
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=amount)

float dispersionMarch (in vec3 rayDirection) {
  vec3 rayOrigin = gPos + -gNor * 0.01;
  rayDirection = normalize(rayDirection);

  float t = 0.0001;

  for (int i = 0; i < 20; i++) {
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

  vec3 reflectionPoint = gPos - gNor * 0.1 + rd * d;
  vec3 reflectionPointNor = getNormal2(reflectionPoint, 0.001);
  dNor = reflectionPointNor;
  reflectionPointNor = normalize(reflectionPointNor);

  // float sss = saturate(1.0 - pow(length(reflectionPoint - gPos) * 0.48, 10.0));
  vec3 disp = scene(refract(rd, reflectionPointNor, ior), ior);
  disp *= min(1.5, 1.0 / d);

  return disp;
}

#pragma glslify: dispersionStep1 = require(./glsl-dispersion, scene=secondRefraction, amount=amount)

#pragma glslify: gradient = require(./gradient)

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  vec3 color = background;
  color = vec3(0.1);
  color = mix(color, vec3(0), isMaterialSmooth(m, 2.0));
  return color;
}

#pragma glslify: reflection = require(./reflection, getNormal=getNormal, diffuseColor=baseColor, map=map, maxDistance=maxDistance, epsilon=epsilon, maxSteps=maxSteps)

const vec3 glowColor = pow(#ED4F2C, vec3(2.2));

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)
#pragma glslify: matCap = require(./matCap, texture=tMatCap)

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    gPos = pos;
    if (t.x>0.) {
      vec3 color = vec3(0.0);

      vec3 nor = getNormal2(pos, 0.08 * t.x);
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
      const int NUM_OF_LIGHTS = 1;
      const float repNUM_OF_LIGHTS = 1.0; // 0.33333;
      light lights[NUM_OF_LIGHTS];
      lights[0] = light(normalize(vec3(1., .75, 1.)), #ffffff, 1.0);
      // lights[1] = light(normalize(vec3(-1., .75, 0.5)), #ffffff, 1.0);
      // lights[2] = light(normalize(vec3(-0.75, -1.0, 1.0)), #ffffff, 1.0);

      float occ = 1.0; // calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      const float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 1.00;
      float specCo = 1.00;
      float disperCo = 0.5;

      float specAll = 0.0;
      for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
        vec3 lightPos = lights[i].position;
        float dif = diffuse(nor, lightPos);
        float spec = pow(clamp( dot(ref, (lightPos)), 0., 1. ), 16.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        // dif *= saturate(0.0 + softshadow(pos, lightPos, 0.02, 1.75));
        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        lin += specCo * spec * (1. - fre);
        specAll += specCo * spec * (1. - fre);

        // Ambient
        lin += 0.2 * amb * isMaterialSmooth(t.y, 2.0);

        const float conserve = 1.0; // TODO figure out how to do this w/o grey highlights
        color +=
          saturate((conserve * occ * dif * lights[i].intensity) * lights[i].color * diffuseColor)
          + saturate(lights[i].intensity * lin * mix(diffuseColor, #ffffff, 0.4));
      }

      color *= 1.0 / float(NUM_OF_LIGHTS);
      color += 0.75 * vec3(pow(specAll, 8.0));

      color += 0.5 * dispersionStep1(nor, rayDirection, n2);
      // color += 0.90 * dispersion(nor, rayDirection, n2);

      // Fog
      // color = mix(background, color, clamp(1.0 * (maxDistance-t.x) / maxDistance, 0., 1.));
      // color *= exp(-t.x * 0.005);

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
      // color *= mix(vec4(1.), vec4(background, 1), length(uv) / 2.);

      // Glow
      // vec3 glowColor = pow(#ffffff, vec3(2.2));
      // color = mix(vec4(glowColor, 1.0), color, 1. - .99 * clamp(t.z / (1.5 * float(maxSteps)), 0., 1.));

      return color;
    }
}

vec4 oil (in vec3 ro, in vec3 rd, in vec2 uv) {
  vec2 oilUV = uv;

  const float period = 20.0;
  const float transitionTime = 2.0;
  float modTime = mod(time, period);
  float mixT = saturate((modTime - transitionTime) / (period - transitionTime));

  oilUV *= 0.8;

  // Vertical stretch
  oilUV.y *= 0.30 + 0.05 * sin(PI * slowTime);

  // Sway horizontally offset by height
  float shiftX1 = ncnoise2(oilUV + 0.2 * modTime);
  float shiftX2 = ncnoise2(oilUV + 0.2 * (modTime - period));
  oilUV.x += 0.25 * mix(shiftX1, shiftX2, mixT);

  float noiseYOff;
  float noiseYOff1 = ncnoise2(5.0 * oilUV + 0.1 * (modTime));
  float noiseYOff2 = ncnoise2(5.0 * oilUV + 0.1 * (modTime - period));
  oilUV.y += dot(oilUV, vec2(2)) + 1.5 * mix(noiseYOff1, noiseYOff2, mixT);

  vec2 uvWarp = oilUV; //  + 0.5 * vec2(vfbm4(oilUV), vfbm4(oilUV + vec2(123.0, 23423.0)));

  float n;
  float n1 = vfbm6(4.0 * uvWarp + .2 * modTime + ncnoise2(40.0 * uvWarp + vec2(2342.34, 234.534)));
  float n2 = vfbm6(4.0 * uvWarp + .2 * (modTime - period) + ncnoise2(40.0 * uvWarp + vec2(2342.34, 234.534)));
  n = mix(n1, n2, mixT);

  vec3 nor = normalize( vec3( dFdx(n)*resolution.x, 1.0, dFdy(n)*resolution.y  )  );
  const vec3 light = vec3(1, 0, 1);
  vec3 ref = reflect(nor, rd);

  float v;
  // v = dot(nor, -rd);
  v = n;

  // V modifiers
  v *= 1.6;

  vec3 color;

  // B&W
  // color = vec3(v);

  // Spectrum
  color = 0.5 + vec3(0.6, 0.5, 0.5) * cos(TWO_PI * (v + vec3(0.15, 0.275, 0.67)));
  // color = 0.5 + 0.5 * cos(TWO_PI * (v + vec3(0.0, 0.33, 0.67)));

  // Dark base on red component of cosine palette
  color *= 0.80 + 0.20 * cos(TWO_PI * (v + 0.15));

  // Lighten based on green
  color += 0.20 * cos(TWO_PI * (v + 0.22));

  // Red only colors
  color.r *= 0.3 + min(1.0, length(color.gb));

  // Specular
  // color += 0.5 * pow(saturate(dot(ref, (lightPos))), 64.0);

  // float l = length(color);
  // if (l >= 1.732) {
  //   color = #ff00ff;
  // } else if (l <= 0.002) {
  //   color = #00FF00;
  // }

  return vec4(color, 1.0);
}

float noiseRain (in vec2 uv, in float offset) {
  uv *= 10.0 - offset;
  uv.y *= 0.125 - 0.1 * offset;

  float noiseScale = 10.0 + offset;
  uv += 0.5000 * cos( 4.0 * uv + ncnoise2(noiseScale * uv + time));
  uv += 0.2500 * cos( 8.0 * uv + ncnoise2(noiseScale * uv + time));
  uv += 0.1250 * cos(16.0 * uv + ncnoise2(noiseScale * uv + time));
  uv += 0.0625 * cos(32.0 * uv);

  float v = ncnoise2(uv);
  v *= pow(v, 2.0);
  v *= 0.5 + 0.5 * ncnoise2(2.0 * uv);

  return v;
}

float gridMask (in vec2 uv, in float size) {
  vec2 c = floor(uv / size) * size;

  float mask = ncnoise2(0.3 / size * c);
  mask *= mask; mask *= mask;

  return step(4.0 * size, mask);
}

#pragma glslify: hueToIOR = require(./dispersion-ray-direction)

vec3 sineTexture (in vec3 rd) {
  // return 0.5 + 0.5 * sin(rd);
  return hsv(vec3(0.75, 2, 2) * vec3(rd));
}

vec3 dispersionColor (in float hue, in vec3 nor , in vec3 rd, in float n2) {
    float ior = hueToIOR(hue, n2, 1.0, amount);
    float ior2 = hueToIOR(hue, 1.0, n2, amount);
    vec3 refracted1 = refract(rd, nor, ior);
    vec3 refracted2 = refract(refracted1, vec3(0, 0, 1), ior2);

    return sineTexture(refracted1);
}

vec3 intDispersion (in vec3 nor, in vec3 rd, in float n2) {
  vec3 color = vec3(0);

  const int stepAmount = 20;
  for (int i = 0; i < 360; i += stepAmount) {
    color += float(stepAmount) * dispersionColor(float(i), nor, rd, n2);
  }

  return color / 360.0;
}

vec4 noiseTexture (in vec3 ro, in vec3 rd, in vec2 uv) {
  const float scale = 1.55;

  const float period = 4.0;
  const float transitionTime = 1.5;
  float modTime = mod(slowTime, period);
  float mixT = saturate((modTime - transitionTime) / (period - transitionTime));

  vec2 UV = uv;

  vec3 q = vec3(0);
  vec4 nUV1 = vec4(0.5 * uv, 0.33333 * modTime, ncnoise2(uv));
  nUV1 += 0.250 * cos(2.0 * nUV1.yzwx);
  nUV1 += 0.125 * cos(4.0 * nUV1.yzwx);

  vec4 nUV2 = vec4(0.5 * uv, 0.33333 * (modTime - period), ncnoise2(uv));
  nUV2 += 0.250 * cos(2.0 * nUV2.yzwx);
  nUV2 += 0.125 * cos(4.0 * nUV2.yzwx);

  vec4 nUV = mix(nUV1, nUV2, mixT);

  float n = vfbmWarp(nUV.xyz, q);
  // n *= 0.85;

  vec3 nor = normalize(vec3(dFdx(n)*resolution.x, 1.0, dFdy(n)*resolution.y));

  vec3 color;
  color = vec3(n);

  const vec3 light = normalize(vec3(1, 1, 1));
  color *= saturate(0.25 + 0.75 * dot(nor, light));
  color += 0.5 * intDispersion(nor, rd, n2);

  color *= 0.75 + 0.25 * cos(TWO_PI * (vec3(0.5 * UV, 0.0) + vec3(0.0, 0.33, 0.67)));

  return vec4(color, 1.0);
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
            color += sample(ro, rd, uv);
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(uv, projectionMatrix);
    rd = (vec4(rd, 1.) * cameraMatrix).xyz;
    rd = normalize(rd);
    gl_FragColor = sample(ro, rd, uv);
    #endif

    // gamma
    gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(0.454545));

    // 'Film' Noise
    gl_FragColor.rgb += 0.01 * (cnoise2((500. + 60.1 * time) * uv + sin(uv + time)) + cnoise2((500. + 300.0 * time) * uv + 253.5));
}
