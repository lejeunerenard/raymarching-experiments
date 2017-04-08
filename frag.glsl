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

float isMaterialSmooth( float m, float goal ) {
  const float eps = .1;
  return 1. - smoothstep(0., eps, abs(m - goal));
}

#pragma glslify: hsv = require(glsl-hsv2rgb)

vec3 baseColor (in vec3 p, in vec3 nor, in vec3 rd, float m, float trap) {
  vec3 color = vec3(.7);

  float n = clamp(1. + dot(rd, nor), 0., 1.) - .15;
  n = smoothstep(.2, 1., n);
  n += .75 * snoise3(p * 1.3);
  // float n = clamp(1. + dot(vec3(-1., 0., 0.), p), 0., 1.);
  // float n = .5 * p.x;
  // color = vec3(.5) + vec3(.5) * cos( 2. * PI * ( vec3(1.) * n + vec3(0., .33, .67) ) );
  // float mask = clamp(pow(smoothstep(.1, 1., 1. + dot(rd, nor)), .8), 0., 1.);
  // color = mask * hsv(vec3(1. + n, .75, 1.));
  // color = mix(color, hsv(vec3(1. + n, .9, 1.)), mask);
  color = hsv(vec3(0.69 + .5 * n, .9, 1.));

  return clamp(color, 0., 1.);
}

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
      const float n1 = 1.000277; // Air
      const float n2 = 2.42; // Diamond
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

      // Fog
      color = mix(background, color, clamp(1.1 * ((maxDistance-t.x) / maxDistance), 0., 1.));

      // Inner Glow
      // vec3 glowColor = #6699FF * 5.0;
      // float fGlow = clamp(t.w * 0.1, 0.0, 1.0);
      // fGlow = pow(fGlow, 3.5);
      // color += glowColor * 3.5 * fGlow;

      color *= exp(-t.x * .2);

      // colorMap(color);

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

float mask1 (in float time, in vec2 uv, in float minRadius, in float maxRadius) {
  time = clamp(time, 0., 1.);

  float radius = min((maxRadius - minRadius) * min(bounce(time), 1.) + minRadius, maxRadius);
  const float eps = .005;

  vec2 absUV = abs(uv);

  float metric = max(absUV.x, absUV.y);

  return 1. - smoothstep(radius - eps, radius, metric);
}
float mask1 (in float time, in vec2 uv) {
  time = clamp(time, 0., 1.);

  const float minRadius = .3;
  const float maxRadius = .65;
  return mask1(time, uv, minRadius, maxRadius);
}

float mask2 (in float time, in vec2 uv) {
  time = clamp(time, 0., 1.);

  const float minRadius = .1;
  const float maxRadius = .5;
  float radius = min((maxRadius - minRadius) * min(bounce(time), 1.) + minRadius, maxRadius);
  const float eps = .005;

  vec2 absUV = abs(uv);

  float metric = length(absUV);

  return 1. - smoothstep(radius - eps, radius, metric);
}

void delayRotate(in float t, inout vec2 uv) {
  t = clamp(t, 0., 1.);

  float a = .5 * PI * circ(t);
  float c = cos(a);
  float s = sin(a);
  mat2 rot = mat2(
     c, s,
    -s, c);

  uv *= rot;
}

vec3 scene (in vec3 pos) {
  // return .5 + .5 * sin(2. * PI * pos);
  return vec3(length(pos.zy));
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

  return n;
}

vec3 scene (in vec2 pos) {
  const float d = .15;
  const float r = .9;
  // return vec3(smoothstep(r, r + d, length(pos)));
  return vec3(smoothstep(0., .35, cnoise2(2. * pos)));
}

float iorRatio (in float frequency) {
  const float green = 588.;
  const float greenIOR = 1.57;
  const float airIOR = 1.;
  const float dispersionFactor = .03;

  float ior = greenIOR + dispersionFactor * sin(.5 * PI * (frequency - green) / 125. );
  return airIOR / ior;
}

float normalNoise (in vec2 p) {
  p *= 0.03125;
  p += vec2(slowTime, smoothstep(-1., 1., cnoise2(p)) * slowTime);

  vec2 q = vec2(
    fbm(p + vec2(0., slowTime)),
    fbm(p + vec2(1., 3.4 * sin(.1 * slowTime))));

  vec2 r = vec2(
    fbm(p + 4.0 * q + vec2(6.)),
    fbm(p + 4.0 * q + vec2(103., 245.)));

  vec2 s = vec2(
    fbm(p + 4.0 * r + vec2(1.34)),
    fbm(p + 4.0 * r + vec2(13., 275.)));

  return fbm(p + 4.0 * s);
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
    cnoise2(1. * (eye + 20. * slowTime + cnoise2(eye + slowTime))));
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
  // vec3 nor = vec3(
  //   cnoise3(rd + eps.xyy) - cnoise3(rd - eps.xyy),
  //   cnoise3(rd + eps.yxy) - cnoise3(rd - eps.yxy),
  //   cnoise3(rd + eps.yyx) - cnoise3(rd - eps.yyx));

  nor = normalize(nor);

  const float between = .08;
  const float greenIOR = 1.57;
  float redIORRatio = 1./(greenIOR - 2. * between); // iorRatio(461.);
  float yellowIORRatio = 1./(greenIOR - 1. * between); // iorRatio(508.);
  float greenIORRatio = 1./greenIOR; // iorRatio(588.);
  float cyanIORRatio = 1./(greenIOR + 1. * between); // iorRatio(609.);
  float blueIORRatio = 1./(greenIOR + 2. * between); // iorRatio(631.);
  float purpleIORRatio = 1./(greenIOR + 3. * between); // iorRatio(674.);

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

  color = vec3(R, G, B);

  const vec3 light = normalize(vec3(0.5, .95, -1.));
  vec3 ref = reflect(eye3D, nor);
  float spec = pow(clamp( dot(ref, light), 0., 1. ), 4.);
  color = mix(color, color * vec3(spec), .85);

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

    // float mask = 0.; // Block out everything
    // float eps = .01;
    // float dist = length(uv);

    // const int circles = 6;
    // float startAmount = .9;
    // float endAmount = .001;
    // float delta = (startAmount - endAmount) / float(circles);

    // // Circles
    // for(int i = 0; i < circles; i++) {
    //   float radius = startAmount - float(i) * delta;
    //   mask += (1. - 2. * mod(float(i), 2.)) * (1. - smoothstep(radius, radius + eps, dist));
    // }

    // gl_FragColor *= clamp(mask, 0., 1.);

    // gamma
    gl_FragColor.rgb = pow(gl_FragColor.rgb, vec3(0.454545));

    // gl_FragColor *= smoothstep(.25, .3, abs(cnoise2(8. * uv)));
}
