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

// KIFS
uniform mat4 kifsM;
uniform float scale;
uniform vec3 offset;

// Greatest precision = 0.000001;
uniform float epsilon;
#define maxSteps 1024
#define maxDistance 50.
#pragma glslify: import(./background)

#define slowTime time * .01
#define Iterations 9

vec3 lightPos = normalize(vec3(1., 0., 0.));

const vec3 un = vec3(1., -1., 0.);

// Utils
#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)
#pragma glslify: cnoise3 = require(glsl-noise/classic/3d)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: rot4 = require(./rotation-matrix4.glsl)
#pragma glslify: bounce = require(glsl-easings/bounce-out)
#pragma glslify: circ = require(glsl-easings/circular-in-out)
#pragma glslify: rotateEase = require(glsl-easings/quintic-in-out)

// Folds
#pragma glslify: fold = require(./folds)
#pragma glslify: foldInv = require(./foldInv)
#pragma glslify: sphereFold = require(./sphere-fold)
#pragma glslify: twist = require(./twist)
void foldNd (inout vec3 z, vec3 n1) {
  z-=2.0 * min(0.0, dot(z, n1)) * n1;
}

float hash(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

// Orbit Trap
float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
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

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

#pragma glslify: fold = require(./folds)

vec3 map (in vec3 p) {
  // Sphere
  // vec3 s = vec3(length(p) - 1., 1., 0.);
  // return s;

  // Square
  float minD = maxDistance;
  p = octahedronFold(p, minD);
  float octaScale = pow(scale, -float(octaPreFold));
  // return vec3(sdBox(p, vec3(1.)) * octaScale, 1., 0.);

  // Fractal
  vec2 f = dodecahedron(p);
  vec3 fractal = vec3(f.x * octaScale, 1., min(minD, f.y));

  // vec2 f2 = mengersphere(p);
  // vec3 fractal2 = vec3(f2.x, 1., f2.y);

  return fractal;
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

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#ef78FF, color, 1. - l * .0625);
  // Dark
  color = mix(#043210, color, clamp(exp(l) * .325, 0., 1.));
}

#pragma glslify: hsv = require(glsl-hsv2rgb)

vec3 scene (in vec3 pos) {
  // return .5 + .5 * sin(2. * PI * pos);
  return vec3(length(pos.zy + vec2(slowTime, 0.)));
  // return vec3(mod(pos.x, .5) + mod(pos.y, .5));
  // return vec3(dot(vec3(0., 0., -1.), pos));
}

float fbm (in vec2 p) {
  float n = 0.5 * cnoise2(p);
  p *= 2.01;

  n += 0.25 * cnoise2(p);
  p *= 2.03;

  n += 0.125 * cnoise2(p);
  p *= 2.02;

  // n += 0.0625 * cnoise2(p);
  // p *= 2.04;

  // n += 0.03125 * cnoise2(p);
  // p *= 2.01;

  // n += 0.015625 * cnoise2(p);
  // p *= 2.02;

  return n;
}

vec3 scene (in vec2 pos) {
  const float d = .15;
  const float r = .9;
  // return vec3(smoothstep(r, r + d, length(pos)));
  return vec3(smoothstep(0., .35, cnoise2(2. * pos + 5. * slowTime)));
}

float normalNoise (in vec2 p) {
  p *= .5;
  // p.x += 18.;
  p += vec2(slowTime, smoothstep(-1., 1., cnoise2(p)) * slowTime);

  // return cos(p.x + fbm(p));
  return .5 + .5 * abs(fbm(p));
  // return .5 + .5 * abs(fbm(p)) + .00005 * cnoise2(500. * p - 134.2);
}
vec3 makeNormal (in vec2 eye) {
  const vec2 eps = vec2(.01, 0.);

  // return vec3(
  //   normalNoise(eye + eps.xy) - normalNoise(eye - eps.xy),
  //   normalNoise(eye + eps.yx) - normalNoise(eye - eps.yx),
  //   cnoise2(eye + cnoise2(eye + slowTime)));
  return vec3(
    normalNoise(eye + eps.xy) - normalNoise(eye - eps.xy),
    normalNoise(eye + eps.yx) - normalNoise(eye - eps.yx),
    0.);
  // return vec3(
  //   normalNoise(eye + eps.xy) - normalNoise(eye - eps.xy),
  //   normalNoise(eye + eps.yx) - normalNoise(eye - eps.yx),
  //   normalNoise(eye.yx + eps.yx) - normalNoise(eye.yx - eps.yx));
}

vec4 sample (in vec3 ro, in vec3 rd, in vec2 uv) {
  vec3 color = vec3(0.);
  vec2 eye = rd.xy;

  const vec2 eps = vec2(.01, 0.);
  // vec3 nor = vec3(0., 0., 1.);
  vec3 nor = makeNormal(eye);

  nor = normalize(nor);

  const float between = .08;
  const float greenIOR = 1.57;
  float redIORRatio = 1./(greenIOR - 2. * between);
  float yellowIORRatio = 1./(greenIOR - 1. * between);
  float greenIORRatio = 1./greenIOR;
  float cyanIORRatio = 1./(greenIOR + 1. * between);
  float blueIORRatio = 1./(greenIOR + 2. * between);
  float purpleIORRatio = 1./(greenIOR + 3. * between);

  vec3 eye3D = vec3(eye, 1.);
  vec3 redRefract = refract(eye3D, nor, redIORRatio);
  vec3 yellowRefract = refract(eye3D, nor, yellowIORRatio);
  vec3 greenRefract = refract(eye3D, nor, greenIORRatio);
  vec3 cyanRefract = refract(eye3D, nor, cyanIORRatio);
  vec3 blueRefract = refract(eye3D, nor, blueIORRatio);
  vec3 purpleRefract = refract(eye3D, nor, purpleIORRatio);

  float r = scene(eye + redRefract.xy).r * 0.5;
  float y = dot(scene(eye + yellowRefract.xy), vec3(2.0, 2.0, -1.0)) / 6.0;
  float g = scene(eye + greenRefract.xy).g * 0.5;
  float c = dot(scene(eye + cyanRefract.xy), vec3(-1.0, 2.0, 2.0)) / 6.0;
  float b = scene(eye + blueRefract.xy).b * 0.5;
  float p = dot(scene(eye + purpleRefract.xy), vec3(2.0, -1.0, 2.0)) / 6.0;

  float R = r + (2.0*p + 2.0*y - c)/3.0;
  float G = g + (2.0*y + 2.0*c - p)/3.0;
  float B = b + (2.0*c + 2.0*p - y)/3.0;

  vec3 refractC = vec3(R, G, B);

  const vec3 light = normalize(vec3(0., .25, -1.));

  vec3 diffuseColor = #00FF00;

  color = vec3(1.);

  vec3 diff = vec3(.3) + vec3(.7) * diffuseColor * vec3(1.) * max(dot(nor, light), 0.);

  vec3 bdrf; // ganked from iq
  // bdrf  = vec3(0.70,0.90,0.95)*(nor.y*0.5+0.5);
  // bdrf += vec3(0.15,0.10,0.05)*diff;
  bdrf  = #c0efff*(nor.y*0.5+0.5);
  bdrf += vec3(0.15,0.10,0.05)*diff;
  color *= 1.0*bdrf;

  color = mix(color, refractC, .4);
  // color = nor;

  #if 0
  if (color == vec3(1.)) {
    color = vec3(1., 0., 1.);
  }
  if (color == vec3(0.)) {
    color = vec3(0., 1., 0.);
  }
  #endif

  return vec4(color, 1.);
}

void main() {
    vec3 ro = cameraRo + cOffset;

    vec2 uv = fragCoord.xy;

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

    // delayRotate(1. * (time - 4.), uv);
    // float trans = clamp(time - 6., 0., 1.);

    float mask = 0.; // Block out everything
    float eps = .001;
    uv = abs(uv);
    float dist = max(uv.x, uv.y);

    const int circles = 2;
    float startAmount = .65;
    float endAmount = .5;
    float delta = (startAmount - endAmount) / float(circles);

    // Circles
    for(int i = 0; i < circles; i++) {
      float radius = startAmount - float(i) * delta;
      mask += (1. - 2. * mod(float(i), 2.)) * (1. - smoothstep(radius, radius + eps, dist));
    }

    mask = 1. - clamp(mask, 0., 1.);
    gl_FragColor.rgb *= mask;
    gl_FragColor.rgb += (1. - mask);

    // gamma
    gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(0.454545));

    gl_FragColor += .05 * cnoise2(500. * uv) + .05 * cnoise2(500. * uv + 253.5);
}
