#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)

// #define debugMapCalls
// #define debugMapMaxed
// #define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
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
#define maxSteps 1024
#define maxDistance 50.0
#pragma glslify: import(./background)

#define slowTime time * .01
#define Iterations 8

vec3 lightPos = normalize(vec3(1., 0., 0.));

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: voronoi = require(./voronoi)

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

  f += 0.5000*noise( p ); p = p*2.02;
  f += 0.2500*noise( p ); p = p*2.03;
  f += 0.1250*noise( p ); p = p*2.01;
  f += 0.0625*noise( p );

  return f * 1.066667;
}

// Orbit Trap
float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
  vec2 u = max(vec2(r - a,r - b), vec2(0));
  return max(r, min (a, b)) - length(u);
}

// IQ's capsule
float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}
float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}
float sdHexPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

#pragma glslify: mandelbox = require(./mandelbox, trap=Iterations, maxDistance=maxDistance, foldLimit=1., s=scale, minRadius=0.5, rotM=kifsM)
// #pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM)

#pragma glslify: dodecahedron = require(./dodecahedron, Iterations=Iterations, scale=scale, kifsM=kifsM)
#pragma glslify: mengersphere = require(./menger-sphere, intrad=1., scale=scale, kifsM=kifsM)

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

#pragma glslify: fold = require(./folds)
#pragma glslify: foldNd = require(./foldNd)
#pragma glslify: twist = require(./twist)

vec3 equation (in vec3 p, in float dt) {
  const float a = 2.4;
  const float b = -3.78;
  const float c = 14.0;
  const float d = -11.0;
  const float e = 4.0;
  const float f = 5.58;
  const float r = 1.0;

  float dxdt = dot(p.yx, vec2(a, b)) + c * p.y * p.z;
  float dydt = dot(p.yz, vec2(d, -1)) + e * p.x * p.z;
  float dzdt = f * p.z + r * p.x * p.y;

  return p + dt * vec3(dxdt, dydt, dzdt);
}
vec3 chen (in vec3 p, in float dt) {
  const float a = 40.0;
  const float b = 28.0;
  const float c = 3.0;

  float dxdt = dot(p.yx, vec2(a, -a));
  float dydt = (c - a) * p.x - p.x * p.z + c * p.y;
  float dzdt = p.x * p.y - b * p.z;

  return p + dt * vec3(dxdt, dydt, dzdt);
}

vec3 lorenz (in vec3 p, in float dt) {
  const float r = 28.0;
  const float s = 10.0;
  const float b = 2.666667; // 8/3

  float dxdt = s * dot(p.yx, vec2(1.0, -1.0));
  float dydt = p.x * ( r - p.z ) - p.y;
  float dzdt = p.x * p.y - b * p.z;

  return p + dt * vec3(dxdt, dydt, dzdt);
}

vec3 beep (in vec3 p, in float dt) {
  const float a = 2.1;
  const float b = 5.1;
  const float c = 0.5;

  float dxdt = b * cos( p.x + p.z * cos(p.y) );
  float dydt = c * cos( p.y + p.x * cos(p.x) );
  float dzdt = a * sin(dot(p, vec3(1.)) + exp(p.x));

  return p + dt * vec3(dxdt, dydt, dzdt);
}

#define MAX_IT 100

vec3 de (in vec3 x, in float endt, in float dt) {
  float t = 0.0;

  for (int i = 0; i < MAX_IT; i ++) {
    if (t >= endt) break;

    // x = noiseField(x, dt);

    t += dt;
  }

  return x;
}

vec2 dMin (vec2 d1, vec2 d2) {
  return (d1.x < d2.x) ? d1 : d2;
}
vec3 dMin (vec3 d1, vec3 d2) {
  return (d1.x < d2.x) ? d1 : d2;
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

// --------------------------------------------------------
// http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
// --------------------------------------------------------

mat3 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat3(
        oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
        oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
        oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c
    );
}

#pragma glslify: dodeca = require(./dodec)

vec3 quartz (in vec3 p, in float angle) {
  vec3 quartzBody = vec3(sdHexPrism(p.xzy - vec3(0., -0.30, 0.0), vec2(0.5, 3.)), 1., 0.);

  const float dodecaScale = 1.00;
  vec3 dp = p;
  dp.y *= 0.416;
  dp *= rotationMatrix(vec3(0.000,PHI,1.), PI / 5. + angle) / dodecaScale;
  vec3 d = vec3(dodeca(dp, 1.) * dodecaScale, 1., 0.);

  return dMax(quartzBody, d);
  // return d;
}

vec3 cluster (in vec3 x) {
  const float scale = 0.4;
  x /= scale;

  const float angleSpan = PI * 0.85;

  vec3 res = vec3(1000.0, 1., 0.);
  for (int k=-1; k<=1; k++) {
    for (int j=-1; j<=1; j++) {
      for (int i=-1; i<=1; i++) {
        vec3 b = 1.5 * vec3(i, j, k);

        vec3 p = x + b;
        p *= rotationMatrix(vec3(
          noise(b),
          noise(1.1 * b - 20.3),
          noise(0.9 * b + 414.9)), angleSpan * noise(20.0 * b));

        vec3 d = quartz(p, PI / 4.0 * noise(13.0 * b));
        d.x *= scale;
        res = dMin( res, d );
      }
    }
  }
  return res;

  return res;
}

// Return value is (distance, material, orbit trap)
vec3 map (in vec3 p) {
  p *= globalRot;

  vec3 outD = vec3(10000., 0., 0.);

  vec3 cl = cluster(p);
  outD = dMin(outD, cl);

  vec3 s = vec3(length(p) - 0.65, 2., 0.);
  s.x -= 0.7 * iqFBM(p + iqFBM(p + iqFBM(p + 20.)));
  // s.x -= 0.9 * voronoi(p - voronoi(2.0 * p));
  s.x *= 0.8;

  float dPre = outD.x;
  outD = dMin(outD, s);
  outD.x = fOpUnionRound(dPre, s.x, 0.0625);

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
  vec2 e = vec2(1.,-1.) * .015 * eps;
  return normalize(vec3(
    map(p + e.xyy).x - map(p - e.xyy).x,
    map(p + e.yxy).x - map(p - e.yxy).x,
    map(p + e.yyx).x - map(p - e.yyx).x));
}

// Material Functions
float diffuse (in vec3 nor, in vec3 lightPos) {
  return clamp(dot(lightPos, nor), 0., 1.);
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
const float n2 = 2.57;

vec3 textures (in vec3 rd) {
  vec3 color = vec3(0.);

  // float v = 2.1 * noise(rd);
  float v = cnoise3(5. * rd);
  v = smoothstep(-1.0, 1.0, v);

  // vec3 maxRd = abs(rd);
  // float v = max(maxRd.x, maxRd.y);

  color = vec3(v);
  // color = mix(#FF1F99, #FF7114, v);
  // color = .5 + vec3(.5, .3, .6) * cos(TWO_PI * (v + vec3(0.0, 0.33, 0.67)));
  color += #61FF77 * (0.5 * abs(rd.y));

  // color *= 1.1 * cos(rd);

  return clamp(color, 0., 1.);
}

vec3 scene (in vec3 rd) {
  vec3 color = vec3(0.);

  rd = normalize(rd);
  color = textures(rd);

  return color;
}

#pragma glslify: dispersion = require(./glsl-dispersion, scene=scene, amount=0.5)
bool isMaterial( float m, float goal ) {
  return m < goal + 1. && m > goal - .1;
}
float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

vec3 baseColor(in vec3 pos, in vec3 nor, in vec3 rd, in float m, in float trap) {
  // return #ffffff;
  return dispersion(nor, rd, n2) * isMaterialSmooth(m, 1.) + #7D8091 * isMaterialSmooth(m, 2.);
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
      vec3 color = background;

      vec3 nor = getNormal(pos, .0001);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = baseColor(pos, nor, rayDirection, t.y, t.w);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);
      float spec = pow(clamp( dot(ref, (lightPos)), 0., 1. ), 4.);
      const float ReflectionFresnel = pow((n1 - n2) / (n1 + n2), 2.);
      float fre = ReflectionFresnel + pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 5.) * (1. - ReflectionFresnel);

      dif *= min(0.1 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      vec3 lin = vec3(0.);

      // Specular Lighting
      lin += spec * (1. - fre) * dif * color.g;
      lin += fre * occ;

      // Ambient
      lin += 0.04 * amb * occ * #ccccff;

      float conserve = (1. - (dot(lin, vec3(1.)) * .3333));
      lin += conserve * dif;

      color *= lin;

      color += .09 * reflection(pos, ref) * isMaterialSmooth(t.y, 1.);

      // Fog
      color = mix(background, color, clamp(1.1 * ((maxDistance-t.x) / maxDistance), 0., 1.));
      color *= exp(-t.x * .1);

      // Inner Glow
      color += innerGlow(t.w);

      // Post process
      // colorMap(color);

      // Debugging
      #ifdef debugMapMaxed
      if (t.z / float(maxSteps) > 0.9) {
        color = vec3(1., 0., 1.);
      }
      #endif

      #ifdef debugMapCalls
      color = vec3(t.z / float(maxSteps));
      #endif

      return vec4(color, 1.);
    } else {
      vec4 color = vec4(background, 0.);
      // Radial Gradient
      // color.xyz *= mix(vec3(1.), background, length(uv) / 2.);

      // Glow
      // color = mix(vec4(1.), color, 1. - .95 * clamp(t.z / float(maxSteps), 0., 1.));
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

    lightPos = normalize(vec3(1., .75, 0.));

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
    // gl_FragColor.rgb += .03 * (cnoise2((500. + 1.1 * time) * uv + sin(uv + time)) + cnoise2((500. + time) * uv + 253.5));
}
