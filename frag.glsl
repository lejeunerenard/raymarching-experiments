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

#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)

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

mat3 rotation3 (float time, float tOff) {
  float rtime1 = PI + 1. * PI * .23 * sin((time + tOff) * .12);
  float rtime2 = PI + 1. * PI * .13 * sin((time + tOff + 1.) * .93);
  float rtime3 = PI + 1. * PI * .23 * sin((time + tOff) * .71);
  return
    mat3(cos(rtime1),0,sin(rtime1),0,1,0,-sin(rtime1),0,cos(rtime1)) *
    mat3(cos(rtime2),sin(rtime2),.0,-sin(rtime2),cos(rtime2),.0,0,0,1) *
    mat3(1,0,0,0,cos(rtime3),sin(rtime3),0,-sin(rtime3),cos(rtime3));
}

void bfold (inout vec2 p) {
  if (p.x - p.y < 0.) {
    float x1 = p.y;
    p.y = p.x;
    p.x = x1;
  }
}

vec2 kifs( inout vec3 p ) {
  const float scale = 2.;

  int maxI = 0;
  const vec3 Offset = vec3(2., 3.1, .0);

  const int Iterations = 14;

  float trap = maxDistance;
  mat4 rot = rotationMatrix(vec3(0., 1., 0.), PI / 4.);

  float r = dot(p, p);
  for (int i = 0; i < Iterations; i++) {
    vec4 pp = vec4(p, 1.) * rot;
    p = pp.xyz;

    p = abs(p);
    // Folding
    bfold(p.xy);
    bfold(p.xz);
    bfold(p.yz);

    // Stretch
    p *= scale;
    p.xy -= Offset.xy * (scale - 1.);

    if(p.z > .5 * Offset.z * (scale -1.)) {
      p.z -= Offset.z * (scale - 1.);
    }

    trap = min(trap, length(p));
  }

  return vec2((length(p) - .2) * pow(scale, - float(Iterations)), trap);
}

vec3 map (in vec3 p) {
  vec4 pp = vec4(p, 1);
  // vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), 2. * PI * sin(20. * time) / 4.) * pp).xyz;
  vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), PI / 4.) * pp).xyz;

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
      vec3 nor = getNormal(pos, 0.01 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = vec3(0.15, .9, .4);

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      dif *= min(0.4 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      vec3 lin = vec3(0.);
      lin += dif;
      lin += 0.4*amb*vec3(0.50,0.70,1.00)*occ;
      color *= lin;

      // Fog
      // color = mix(background, color, (maxDistance-t.x) / maxDistance);

      // Inner Glow
      vec3 glowColor = vec3(0.1, .4, .7) * 5.0;
      float fGlow = clamp(t.w * 0.1, 0.0, 1.0);
      fGlow = pow(fGlow, 3.0);
      color += glowColor * 15.0 * fGlow;

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
      color *= sqrt(1.75 - 1.25 * dot(rayDirection, vec3(0., 0., -1.)));
      // Glow
      color = mix(vec3(1.), color, 1. - .3 * clamp(t.z / float(maxSteps), 0., 1.));
    }

    return vec4(color, 1.);
}

void main() {
    float cameraTime = time * .35 + PI / 2.;
    float d = 2.5 + clamp(2.5 * sin(cameraTime), -2., 2.);
    float turn = - PI / 2. * clamp(-2. + 3. * sin(cameraTime + PI), 0., 1.);
    vec3 turnAxis = normalize(vec3(0., 1., 1.));
    // float turn = PI / 2.;
    // float d = .5;

    vec3 ro = vec3(0.,0.,d) + cOffset;

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec4 t = vec4(0.);
    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 rd = getRayDirection(vec2(
                  float(x) / resolution.y + fragCoord.x,
                  float(y) / resolution.y + fragCoord.y),
                  projectionMatrix);
            rd = (vec4(rd, 1.) * rotationMatrix(turnAxis, turn)).xyz;
            t = march(ro, rd);
            color += shade(ro, rd, t);
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(fragCoord, projectionMatrix);
    rd = (vec4(rd, 1.) * rotationMatrix(turnAxis, turn)).xyz;
    vec4 t = march(ro, rd);
    gl_FragColor = shade(ro, rd, t);
    #endif
}
