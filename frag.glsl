#define PI 3.1415926536

// #define debugMapCalls
// #define debugMapMaxed
#define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform vec3 cOffset;
uniform vec3 cameraRo;
uniform mat4 cameraMatrix;
uniform mat4 orientation;
uniform mat4 projectionMatrix;
uniform sampler2D texture;

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 1024
#define maxDistance 20.
#define background #dddddd

const vec3 lightPos = vec3(-.6, 0., .6);

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
#pragma glslify: rot4 = require(./rotation-matrix4.glsl)

// Folds
#pragma glslify: fold = require(./folds)
#pragma glslify: foldInv = require(./foldInv)
#pragma glslify: sphereFold = require(./sphere-fold)
#pragma glslify: twist = require(./twist)
void foldNd (inout vec3 z, vec3 n1) {
  z-=2.0 * min(0.0, dot(z, n1)) * n1;
}

#define Iterations 25
#pragma glslify: mandelbox = require(./mandelbox, trap=19, maxDistance=maxDistance, foldLimit=1.25, s=scale, minRadius=0.1, rotM=kifsM)
#pragma glslify: octahedron = require(./octahedron, scale=scale, kifsM=kifsM)
#pragma glslify: dodecahedron = require(./dodecahedron, scale=scale, kifsM=kifsM)

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

vec3 map (in vec3 p) {
  vec4 pp = vec4(p, 1);
  vec3 q = p;

  // Sphere
  // return vec3(length(q) - 1., 1., 0.);

  // Square
  // return vec3(sdBox(q, vec3(.5)), 1., 0.);

  vec2 fractal = dodecahedron(q);
  return vec3(fractal.x, 1., fractal.y);
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
#
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
#pragma glslify: matCap = require(./matCap, texture=texture)

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#90B8FF, color, 1. - l * .0625);
  // Dark
  color = mix(#7782E8, color, clamp(exp(l) * .325, 0., 1.));
}

float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

vec3 stripsGeneral (in float t) {
  const int num = 2;
  vec3 colors[num];
  colors[0] = #FF0000;
  colors[1] = #1485CC;

  const float period = 2.;
  t = mod(t + snoise2(vec2(t, t)), period);

  // Gradient
  const float fraction = period / float(num);
  const float delta = fraction;

  vec3 c = vec3(1.);
  for (int i = 0; i < num; i++) {
    c = mix(c, colors[i], smoothstep(float(i) * fraction - delta * .5, float(i) * fraction + delta * .5, t));
  }

  return c;
}

#define THICKNESS_SCALE 32.0     // film thickness scaling factor
vec3 attenuation(float filmThickness, vec3 wavelengths, vec3 normal, vec3 rd) {
  return 0.5 + 0.5 * cos(((THICKNESS_SCALE * filmThickness)/(wavelengths + 1.0)) * dot(normal, rd));    
}

#pragma glslify: hsv = require(glsl-hsv2rgb)

vec3 baseColor (in vec3 p, in vec3 nor, in vec3 rd, float m) {
  vec3 color = vec3(1.);

  color = #34FFE0;
  color = mix(color, #FF4DAC, (1. + dot(rd, nor)) / 2.);
  color += .45 * hsv(vec3(length(p) + dot(rd, nor), .75, 1.));
  color *= .8;

  return clamp(color, 0., 1.);
}

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv ) {
    vec3 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, .001 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = baseColor(pos, nor, rayDirection, t.y);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);
      float spec = pow(clamp( dot(ref, lightPos), 0., 1. ), 32.);
      const float ReflectionFresnel = 0.99;
      float fre = ReflectionFresnel * pow(clamp( 1. + dot(nor, rayDirection), 0., 1. ), 2.) + (1. - ReflectionFresnel);

      dif *= min(0.1 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      vec3 lin = vec3(0.);
      lin += 1. * vec3(dif);
      lin += 0.4 * amb * occ * #ffbb66;
      lin += .25 * fre * occ * dif;
      lin += 2. * spec * dif * color.g;
      color *= lin;

      // Fog
      color = mix(background, color, clamp(1.1 * ((maxDistance-t.x) / maxDistance), 0., 1.));

      // Inner Glow
      vec3 glowColor = #FF3356 * 5.0;
      float fGlow = clamp(t.w * 0.1, 0.0, 1.0);
      fGlow = pow(fGlow, 3.5);
      color += glowColor * 3.5 * fGlow;

      color *= exp(-t.x * .1);

      // colorMap(color);

      #ifdef debugMapMaxed
      if (t.z / float(maxSteps) > 0.9) {
        color = vec3(1., 0., 1.);
      }
      #endif

      #ifdef debugMapCalls
      color = vec3(t.z / float(maxSteps));
      #endif
    } else {
      // Radial Gradient
      color *= mix(vec3(1.), background, length(uv) / 2.);

      // Glow
      // color = mix(vec3(1.), color, 1. - .95 * clamp(t.z / float(maxSteps), 0., 1.));
    }

    return vec4(color, 1.);
}

void main() {
    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec4 t = vec4(0.);
    vec2 R = resolution * 2.;

    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 rd = getRayDirection(vec2(
                  float(x) / R.y + uv.x,
                  float(y) / R.y + uv.y),
                  projectionMatrix);
            rd = (vec4(rd, 1.) * cameraMatrix).xyz;
            t = march(ro, rd);
            color += shade(ro, rd, t, uv);
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(uv, projectionMatrix);
    rd = (vec4(rd, 1.) * cameraMatrix).xyz;
    vec4 t = march(ro, rd);
    gl_FragColor = shade(ro, rd, t, uv);
    #endif

    gl_FragColor += .0125 * snoise2(uv.xy * resolution * .1 + 100000. * time);
}
