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

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 512
#define maxDistance 50.0
#pragma glslify: import(./background)

#define slowTime time * .1

vec3 lightPos = normalize(vec3(1., .75, 0.));

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
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
  f += 0.062500*noise( p ); p = p*2.025;

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

// IQ
float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

// p as usual, e exponent (p in the paper), r radius or something like that
// #pragma glslify: octahedral = require(./model/octahedral)
// #pragma glslify: dodecahedral = require(./model/dodecahedral)
// #pragma glslify: icosahedral = require(./model/icosahedral)

float timeIndex (in float x) {
  return 100.0 * x + 3.0; //  + 10.0 * time;
}

float bumpStripes (in float x) {
  float v = 0.75 + 0.25 * sin(x + 0.25 * noise(vec3(20.0 * x, x, 0.0)));
  v = 0.45 + 0.5 * smoothstep(0.55, 0.6, v);
  return v;
}

mat3 stripesRot = mat3(0.);

#define goldenRatio 2.399827721492203
// In radians 👆

float polarDistance (in vec2 p1, in vec2 p2) {
  // p1.x is radius
  // p1.y is angle

  return sqrt( p1.x * p1.x + p2.x * p2.x - 2.0 * p1.x * p2.x * cos(p2.y - p2.y) );
}

float petal (in vec3 p, in float radius) {
  float groundedY = p.y + radius;

  p.x = abs(p.x);

  p.x /= groundedY * sqrt((2.0 - groundedY) / 3.0) + 0.01;
  p.x = max(p.x, 0.0); // clip negative x

  // Z Trimming
  float d = sdBox(p, vec3(1.0, 1.0, 0.05));
  d *= 0.25;

  return d;
}

float phyllotaxis (in vec3 p) {
  const float c = 0.5;
  const float radius = 0.5;

  float r0 = length(p.xy);

  float phi = 0.0;
  float r = 0.0;
  float minD = 10000.0;

  for (float n = 1.0; n < 100.0; n++) {
    float phiN = n * goldenRatio;
    float rN = c * sqrt(n);

    if (n >= 2.0 * time) { break; }

    vec3 z = p;
    if ( rN >= r0 + 3.0 * radius ) { break; }

      vec3 center = vec3(
          rN * cos(phiN),
          rN * sin(phiN),
          0.0);

      z -= center;

      // Rotate space to align w/ phi
      // z *= rotationMatrix(vec3(0.0, 0.0, -1.0), phiN);

      // float d = petal(z - center, radius);
      // float d = sdBox(z, vec3(0.5));
      float d = length(z - center) - radius;
      minD = min(minD, d);
  }

  return minD;
}

float parametric (in vec3 p) {
  float angle = atan(p.y, p.x);
  // angle += TWO_PI * step(0.0, -angle);
  float k = 6.0 + 4.0 * sin(0.1 * slowTime);
  // float k = 4.0;

  float minD = 1000.0;

  for (int i = 0; i < 20; i++) {
    float formula = cos(k * angle);
    vec3 center = vec3(
        formula * cos(angle),
        formula * sin(angle),
        0.0);

    minD = min(minD, length(p - center) - 0.05);
    angle += PI;
  }

  return minD;
}

#pragma glslify: ease = require(glsl-easings/bounce-in)
// Return value is (distance, material, orbit trap)
vec3 map (in vec3 p) {
  vec3 outD = vec3(10000., 0., 0.);

  p *= globalRot;
  vec3 q = p;

  // Timing
  float oscTime = mod(time, 5.0);

  // p += cos(p);

  // float dr = 1.0;
  const float scale = 4.0;
  const float repScale = 0.25;
  float a = scale + 0.5 * noise(p);
  float b = repScale;

  for (int i = 0; i < 4; i++) {
    // dr = length(-a * sin(a * p)) * dr + 1.0;
    p += b * cos(a * p.yzx + time);
    // p *= rotationMatrix(vec3(0.0, 1.0, 0.0), PI * 0.33);
    // p += noise(p);

    a *= scale;
    b *= repScale;
  }

  // Sphere
  vec3 s = vec3(length(p) - 1.0, 2.0, 0.0);
  outD = dMin(outD, s);

  // Safety coefficient
  outD.x *= 0.3;

  // vec3 floor = vec3(sdPlane(q, vec4(0.,1.,0.,0.)), 1.0, 0.0);
  // outD = dMin(outD, floor);

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
  return dot(lightPos, nor);
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

#pragma glslify: hsv = require(glsl-hsv2rgb)
#pragma glslify: checker = require(glsl-checker)
#pragma glslify: debugColor = require(./debug-color-clip)

const float n1 = 1.0;
const float n2 = 2.50;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  float v = noise(20.0 * rd);
  // v = smoothstep(-1.0, 1.0, v);

  color = vec3(v);
  // color = mix(#FF1F99, #FF7114, v);
  // color = .5 + vec3(.5, .3, .6) * cos(TWO_PI * (v + vec3(0.0, 0.33, 0.67)));
  color += #61FF77 * (0.5 * abs(rd.y));

  color *= 1.1 * cos(rd);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=0.75)
bool isMaterial( float m, float goal ) {
  return m < goal + 1. && m > goal - .1;
}
float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

#pragma glslify: gradient = require(./gradient)

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  vec3 color = vec3(0.0);

  color = mix(color, #ffffff, isMaterialSmooth(m, 2.));

  return color;
}

vec4 marchRef (in vec3 rayOrigin, in vec3 rayDirection) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  for (int i = 0; i < maxSteps / 3; i++) {
    vec3 d = map(rayOrigin + rayDirection * t);
    if (d.x < epsilon) return vec4(t + d.x, d.y, float(i), d.z);
    t += d.x;
    maxI = float(i);
    trap = d.z;
    if (t > maxDistance) break;
  }
  return vec4(-1., 0., maxI, trap);
}
vec3 reflection (in vec3 ro, in vec3 rd) {
  rd = normalize(rd);
  vec4 t = marchRef(ro + rd * .09, rd);
  vec3 pos = ro + rd * t.x;
  vec3 color = vec3(0.);
  if (t.x > 0.) {
    vec3 nor = getNormal(pos, .0001);
    color = baseColor(pos, nor, rd, t.y, t.w);
  }

  return color;
}

const vec3 glowColor = #39CCA1;

#pragma glslify: innerGlow = require(./inner-glow, glowColor=glowColor)

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 color = vec3(0.0);
      float isFloor = isMaterialSmooth(t.y, 2.0);
      float isWall = isMaterialSmooth(t.y, 3.0);

      vec3 nor = getNormal2(pos, 0.08 * t.x);
      // vec3 nNorP = pos * 50.0;
      // nor = normalize(nor +
      //   0.09 * vec3(noise(nNorP), noise(nNorP + 203.0), noise(nNorP - 523.0))
      //     * (isFloor + 0.5 * isWall));

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
      const float repNUM_OF_LIGHTS = 0.3333;
      light lights[NUM_OF_LIGHTS];
      lights[0] = light(normalize(vec3(1., .75, 0.)), #ffffff, 0.8);
      lights[1] = light(normalize(vec3(-1., -.5, 0.5)), #ffffff, 0.6);
      lights[2] = light(normalize(vec3(-1., 1.0, -0.5)), #ffffff, 0.2);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      const float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);

      float freCo = 1.0;
      float specCo = 1.0;
      float disperCo = 0.5;

      for (int i = 0; i < NUM_OF_LIGHTS; i++ ) {
        vec3 lightPos = lights[i].position;
        float dif = diffuse(nor, lightPos);
        float spec = pow(clamp( dot(ref, (lightPos)), 0., 1. ), 4.);
        float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

        // dif *= min(0.1 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
        vec3 lin = vec3(0.);

        // Specular Lighting
        fre *= freCo * occ;
        lin += fre;
        lin += specCo * spec * (1. - fre);

        // Ambient
        // lin += 0.01 * amb * occ * #ffcccc;

        const float conserve = 1.0; // TODO figure out how to do this w/o grey highlights
        color +=
          clamp((conserve * dif * lights[i].intensity) * lights[i].color * diffuseColor, 0.0, 1.0)
          + clamp(lights[i].intensity * lin, 0., 1.);

        color += disperCo * repNUM_OF_LIGHTS * lights[i].intensity * dispersion(nor, rayDirection, n2, lights[i].color);
      }
      color *= 2.0 / float(NUM_OF_LIGHTS);

      // color += 0.0125 * reflection(pos, ref);
      // color += 0.50 * dispersion(nor, rayDirection, n2);

      // Fog
      color = mix(background, color, clamp(1.1 * ((maxDistance-t.x) / maxDistance), 0., 1.));
      color *= exp(-t.x * .05);

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

      return vec4(color, 1.);
    } else {
      vec4 color = vec4(background, 0.);
      if (!BLOOM) {
        color.a = 1.0;
      }

      // Radial Gradient
      // color.xyz *= mix(vec3(1.), background, length(uv) / 2.);

      // Glow
      // color = mix(vec4(#E8F6FF, 1.0), color, 1. - .99 * clamp(t.z / (2.0 * float(maxSteps)), 0., 1.));

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
    stripesRot = rotationMatrix(vec3(0.0, 1.0, 0.0), PI * 0.4);

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
    gl_FragColor.rgb += .03 * (cnoise2((500. + 1.1 * time) * uv + sin(uv + time)) + cnoise2((500. + time) * uv + 253.5));
}
