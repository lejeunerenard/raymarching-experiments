#define PI 3.1415926536

// #define debugMapCalls
// #define debugMapMaxed
#define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform vec3 cOffset;
uniform mat4 orientation;
uniform mat4 projectionMatrix;
uniform sampler2D texture;
uniform float d;

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

// Greatest precision = 0.000001;
float epsilon = .001;
const int maxSteps = 512;
float maxDistance = 20.;
const vec3 background = #222222;

const vec3 lightPos = vec3(3., 1., 5.);

mat4 rotationMatrix(vec3 axis, float angle) {
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

#pragma glslify: fold = require(./folds)
#pragma glslify: foldInv = require(./foldInv)

void bfold (inout vec2 p) {
  const vec2 n = vec2(1., -1.);
  p -= min(0., dot(p, n)) * n;
}
void sphereFold (inout vec3 p, inout float dz) {
  const float fixedRadius2 = 1.;
  const float minRadius2 = .25;

  float r2 = dot(p, p);
  if (r2 < minRadius2) {
    float temp = (fixedRadius2 / minRadius2);
    p *= temp;
    dz *= temp;
  } else if (r2 < fixedRadius2) {
    float temp = (fixedRadius2 / r2);
    p *= temp;
    dz *= temp;
  }
} 

// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolarC(inout vec2 p, float repetitions) {
  float angle = 2.*PI/repetitions;
  float a = atan(p.y, p.x) + angle/2.;
  float r = length(p);
  float c = floor(a/angle);
  a = mod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
  // For an odd number of repetitions, fix cell index of the cell in -x direction
  // (cell index would be e.g. -5 and 5 in the two halves of the cell):
  if (abs(c) >= (repetitions/2.)) c = abs(c);
  return c;
}
void pModPolar(inout vec2 p, float repetitions) {
  float angle = 2.*PI/repetitions;
  float a = atan(p.y, p.x) + angle/2.;
  float r = length(p);
  a = mod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
}
vec3 opTwist( vec3 p, float angle ) {
    float c = cos(angle);
    float s = sin(angle);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}

vec2 kifs( inout vec3 p ) {
  int maxI = 0;

  const int Iterations = 20;

  float trap = maxDistance;

  float dz = 1.;

  p = opTwist(p, p.y * 2. * sin(.4 * PI * time));
  // pModPolar(p.yz, 4.);

  for (int i = 0; i < Iterations; i++) {
    p = abs(p);

    // Folding
    fold(p.xy);
    fold(p.xz);
    fold(p.yz);

    p.z -= 1.;
    p = abs(p);
    p.z += 1.;

    // Rot2 & Stretch
    p.xyz = (vec4(p, 1.) * kifsM).xyz;

    if(p.z > .5 * offset.z * (scale -1.)) {
      p.z -= offset.z * (scale - 1.);
    }

    trap = min(trap, length(p));
  }

  return vec2(.25 * (length(p) - .2) * pow(scale, - float(Iterations)) / abs(dz), trap);
}

vec3 map (in vec3 p) {
  vec4 pp = vec4(p, 1);
  vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), 0.2 * PI * time) * pp).xyz;
  // vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), PI / 4.) * pp).xyz;

  vec2 fractal = kifs(q);
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

// Material Functions
float diffuse (in vec3 nor, in vec3 lightPos) {
  return clamp(dot(nor, lightPos) / length(lightPos), 0., 1.);
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

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t ) {
    vec3 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, .01 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = mix(#F80596, #F86805, clamp(length(pos) / 3., 0., 1.));

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      dif *= min(0.1 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      color *= vec3(dif) + (0.8*amb*occ)*vec3(0.50,0.70,1.00);
      color *= 1.3;

      // Fog
      // color = mix(background, color, (maxDistance-t.x) / maxDistance);

      // Inner Glow
      vec3 glowColor = #FFFFFF * 5.0;
      float fGlow = clamp(t.w * 0.1, 0.0, 1.0);
      fGlow = pow(fGlow, 3.0);
      color += glowColor * 7.5 * fGlow;

      color *= exp(-t.x * .1);

      colorMap(color);

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
      color *= sqrt(1.75 + 1.25 * rayDirection.z);
      // Glow
      // color = mix(vec3(1.), color, 1. - .3 * clamp(t.z / float(maxSteps), 0., 1.));
    }

    return vec4(color, 1.);
}

void main() {
    // float cameraTime = time * .35 + PI / 2.;
    // float d = 2.5 + clamp(2.5 * sin(cameraTime), -2., 2.);
    // float turn = - PI / 2. * clamp(-2. + 3. * sin(cameraTime + PI), 0., 1.);
    // vec3 turnAxis = normalize(vec3(0., 1., 1.));

    vec3 turnAxis = normalize(vec3(0., 1., 0.));
    float turn = 0.;

    vec3 ro = vec3(0.,0.,d) + cOffset;

    mat4 cameraMatrix = rotationMatrix(vec3(0., 1., 0.), 0.);

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec4 t = vec4(0.);
    vec2 R = resolution * 2.;

    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 rd = getRayDirection(vec2(
                  float(x) / R.y + fragCoord.x,
                  float(y) / R.y + fragCoord.y),
                  projectionMatrix);
            rd = normalize(vec4(rd, 1.) * cameraMatrix * rotationMatrix(turnAxis, turn)).xyz;
            t = march(ro, rd);
            color += shade(ro, rd, t);
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(fragCoord, projectionMatrix);
    rd = normalize((vec4(rd, 1.) * cameraMatrix * rotationMatrix(turnAxis, turn)).xyz);
    vec4 t = march(ro, rd);
    gl_FragColor = shade(ro, rd, t);
    #endif

    gl_FragColor += .0125 * snoise2(fragCoord.xy * resolution * .1);
}
