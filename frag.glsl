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
uniform bool BLOOM;
uniform vec3 cOffset;
uniform vec3 cameraRo;
uniform mat4 cameraMatrix;
uniform mat4 orientation;
uniform mat4 projectionMatrix;

uniform vec3 objectPos;
uniform float objectR;

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 128
#define maxDistance 200.0
#pragma glslify: import(./background)

#define slowTime time * .05

vec3 lightPos = normalize(vec3(1., .75, 0.));
vec3 gPos = vec3(0.0);
vec3 gNor = vec3(0.0);
vec3 dNor = vec3(0.0);

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: vmax = require(./hg_sdf/vmax)

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

float iqFBM (vec3 p) {
  float f = 0.0;

  f += 0.500000*noise( p ); p = p*2.02;
  f += 0.250000*noise( p ); p = p*2.03;
  f += 0.125000*noise( p ); p = p*2.01;
  // f += 0.062500*noise( p ); p = p*2.025;

  return f * 1.066667;
}

float fbmWarp (vec3 p, out vec3 q) {
  const float scale = 4.0;

  q = vec3(
        iqFBM(p + vec3(0.0, 0.0, 0.0)),
        iqFBM(p + vec3(3.2, 34.5, .234)),
        iqFBM(p + vec3(7.0, 2.9, -2.42)));

  vec3 s = vec3(
        iqFBM(p + scale * q + vec3(23.9, 234.0, -193.0)),
        iqFBM(p + scale * q + vec3(3.2, 852.0, 23.42)),
        iqFBM(p + scale * q + vec3(7.0, -232.0, -2.42)));

  // vec3 x = vec3(
  //       iqFBM(p + scale * s + vec3(-1.0, 73.234, 0.0)),
  //       iqFBM(p + scale * s + vec3(3.2, 34.5, 2.664)),
  //       iqFBM(p + scale * s + vec3(8.0, 2.9, 222.42)));

  return iqFBM(p + scale * s);
}

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

// #pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
// #pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM)

// #pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
// #pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

#define octaPreFold 6
mat4 octaM = mat4(
scale, 0., 0., 0.,
0., scale, 0., 0.,
0., 0., scale, 0.,
1., 1., 1., 1.) * mat4(
1., 0., 0.1, .0,
0., 1., 0., 0.,
0., 0.2, 1., 0.,
0., 0., 0., 1.);
#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=octaPreFold, kifsM=octaM, trapCalc=trapCalc)
// 
// #pragma glslify: fold = require(./folds)
// #pragma glslify: foldNd = require(./foldNd)
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

float gRAngle = TWO_PI * 0.05 * time;
float gRc = cos(gRAngle);
float gRs = sin(gRAngle);
mat3 globalRot = mat3(
  gRc, 0.0, -gRs,
  0.0, 1.0,  0.0,
  gRs, 0.0,  gRc);

#pragma glslify: rotationMatrix = require(./rotation-matrix3)
#pragma glslify: rotMat2 = require(./rotation-matrix2)

// IQ
float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

// p as usual, e exponent (p in the paper), r radius or something like that
#pragma glslify: octahedral = require(./model/octahedral)
// #pragma glslify: dodecahedral = require(./model/dodecahedral)
// #pragma glslify: icosahedral = require(./model/icosahedral)

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
#pragma glslify: ease = require(glsl-easings/bounce-in)
#pragma glslify: voronoi = require(./voronoi)
#pragma glslify: band = require(./band-filter)
#pragma glslify: tetrahedron = require(./model/tetrahedron)
#
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
  vec3 outD = vec3(10000., 0., 0.);

  // p *= globalRot;
  vec3 q = p;

  // Space Warp
  float blockID = pModInterval1(q.z, 0.1, -30.0, 5.0);
  // float blockID = pMod1(q.z, 0.5);

  const int sides = 7;
  const float angle = TWO_PI / float(sides);
  const float width = 0.05;
  q *= rotationMatrix(vec3(0, 0, 1), PI * 0.25 * time * (1.0 + 0.05 * blockID));

  for (int i = 0; i < sides; i++) {
    vec3 rP = q * rotationMatrix(vec3(0, 0, 1), angle * float(i));
    vec3 q_3 = rP;

    q_3.x += 1.00 + 0.5 * sin(2.0 * slowTime * (2.0 + 0.33 * blockID));

    vec3 b2 = vec3(sdBox(q_3, vec3(width, 20.0, width * 0.5)), 1.0, blockID);
    outD = dMin(outD, b2);
  }

  return outD;
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
#pragma glslify: debugColor = require(./debug-color-clip)

const float n1 = 1.0;
const float n2 = 1.55;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float v = cnoise3(3.5 * rd + 2305.0);
  v = smoothstep(-0.1, 0.5, v);

  color = vec3(v);

  float align1 = saturate(dot(-rd, vec3(1.0, 0.0, 0.0)));
  color *= mix(vec3(1.0), #00cc42, 0.4 * align1);

  float align2 = saturate(dot(-rd, vec3(0.0, 1.0, 0.0)));
  color *= mix(vec3(1.0), #cc0000, align2);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=0.055)

float dispersionMarch (in vec3 rayDirection) {
  vec3 rayOrigin = gPos + -gNor * 0.01;
  rayDirection = normalize(rayDirection);

  float t = 0.0001;

  for (int i = 0; i < 40; i++) {
    float d = map(rayOrigin + rayDirection * t).x;
    if (d >= 0.0) break;
    d = min(d, -0.0001);

    t += abs(d);
  }
  return t;
}

vec3 secondRefraction (in vec3 rd) {
  return scene(rd);
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
  vec3 sss = vec3(1.0 - pow(length(reflectionPoint - gPos) * 0.25, 0.125));
  vec3 disp = min(1.5, 1.0 / d) * dispersion(reflectionPointNor, rd, n2, n1);
  // vec3 disp = dispersion(reflectionPointNor, rd, n2, n1);

  return disp + sss;
}

#pragma glslify: dispersionStep1 = require(./glsl-dispersion, scene=secondRefraction, amount=0.055)

#pragma glslify: gradient = require(./gradient)
#pragma glslify: hsv = require(glsl-hsv2rgb)

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  vec3 color = vec3(0.5);

  // Is a fractional number? No
  // color = vec3(mod(trap, 1.0));

  // Short mod of blockID
  // color = hsv(vec3(trap * 0.03125, 0.8, 0.8));

  color = vec3(0.5, 0.5, 0.65) + vec3(0.5, 0.5, 0.35) * cos(TWO_PI * ((trap * 0.03125 + slowTime) + vec3(0.0, 0.33, 0.67)));

  // color = hsv(vec3(pos.z * 0.06 + sin(0.3 * slowTime) + 0.6, 1.0, 1.0));
  return color;
}

// #pragma glslify: reflection = require(./reflection, getNormal=getNormal, diffuseColor=baseColor, map=map, maxDistance=maxDistance, epsilon=epsilon, maxSteps=maxSteps)

const vec3 glowColor = #39CCA1;

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    // gPos = pos;
    if (t.x>0.) {
      vec3 color = vec3(0.0);

      vec3 nor = getNormal2(pos, 0.008 * t.x);
      gNor = nor;

      vec3 ref = reflect(rayDirection, nor);
      ref = normalize(ref);

      // Basic Diffusion
      vec3 diffuseColor = baseColor(pos, nor, rayDirection, t.y, t.w);

      // Declare lights
      struct light {
        vec3 position;
        vec3 color;
        float intensity;
      };
      const int NUM_OF_LIGHTS = 3;
      const float repNUM_OF_LIGHTS = 0.5; // 0.3333;
      light lights[NUM_OF_LIGHTS];
      lights[0] = light(normalize(vec3(1., .75, 1.)), #ffffff, 0.9);
      lights[1] = light(normalize(vec3(-1., .75, 0.5)), #ffffff, 0.9);
      lights[2] = light(normalize(vec3(-0.75, -1.0, 1.0)), #ffccaa, 0.9);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      const float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 0.7;
      float specCo = 0.4;
      float disperCo = 0.5;

      for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
        vec3 lightPos = lights[i].position;
        float dif = diffuse(nor, lightPos);
        float spec = pow(clamp( dot(ref, (lightPos)), 0., 1. ), 4.0);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        dif *= saturate(0.2 + softshadow(pos, lightPos, 0.02, 1.5));
        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * dif * occ;
        lin += fre;
        lin += specCo * spec * (1. - fre);

        // Ambient
        lin += 0.020 * amb;

        const float conserve = 1.0; // TODO figure out how to do this w/o grey highlights
        color +=
          saturate((conserve * dif * lights[i].intensity) * lights[i].color * diffuseColor)
          + saturate(lights[i].intensity * lin * mix(diffuseColor, #ffffff, 0.4));

        // color += disperCo * repNUM_OF_LIGHTS * lights[i].intensity * dispersionStep1(nor, rayDirection, n2, lights[i].color);
      }
      color *= 1.0 / float(NUM_OF_LIGHTS);

      color = 0.8 * diffuseColor;

      // color += 0.05 * reflection(pos, ref);
      // color += 0.7 * dispersionStep1(nor, rayDirection, n2);
      color += 0.5 * dispersion(nor, rayDirection, n2);

      // Fog
      color = mix(background, color, clamp((maxDistance-t.x) / maxDistance, 0., 1.));
      color *= exp(-t.x * .125);

      // Inner Glow
      // color += innerGlow(length(pos));

      // Post process
      // colorMap(color);

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
      // color.xyz *= mix(vec3(1.), background, length(uv) / 2.);

      // Glow
      // color = mix(vec4(#E8F6FF, 1.0), color, 1. - .99 * clamp(t.z / (50.0 * float(maxSteps)), 0., 1.));

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
    // gl_FragColor.rgb += .02 * (cnoise2((500. + 60.1 * time) * uv + sin(uv + time)) + cnoise2((500. + 300.0 * time) * uv + 253.5));
}
