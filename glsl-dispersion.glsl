#ifndef PI
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#endif

// #define RGBCMY 1
// #define REFR_INTEGRAL 1
#define HUE 1
#define HUE_NUM 20
// #define COS_HUE 1
#pragma glslify: hsv = require(glsl-hsv2rgb)

#pragma glslify: hue2IOR = require(./dispersion-ray-direction)
// #pragma glslify: hue2IOR = require(./dispersion/hue-to-ior-exponential)
// #pragma glslify: hue2IOR = require(./dispersion/hue-to-ior-sigmoid)
// #pragma glslify: hue2IOR = require(./dispersion/hue-to-ior-polynomial)

vec3 intRefract (in vec3 l, in vec3 n, in float r1, in float r2) {
  float sqrt1minusC = sqrt(1.0 - dot(l, -n));
  return 0.5 * (
      (r2 * refract(l, n, r2)
        - n * asin(sqrt1minusC * r2) / sqrt1minusC)
    -
      (r1 * refract(l, n, r1)
        - n * asin(sqrt1minusC * r1) / sqrt1minusC));
}

vec3 refractColors (in vec3 nor, in vec3 eye, in float n2, in float n1, in vec3 lightColor) {
  float between = amount;
  float greenIOR = n2;

  #ifdef RGBCMY
  float redIORRatio = hue2IOR(0.0, greenIOR, n1, between);
  float yellowIORRatio = hue2IOR(60.0, greenIOR, n1, between);
  float greenIORRatio = hue2IOR(120.0, greenIOR, n1, between);
  float cyanIORRatio = hue2IOR(180.0, greenIOR, n1, between);
  float blueIORRatio = hue2IOR(240.0, greenIOR, n1, between);
  float purpleIORRatio = hue2IOR(270.0, greenIOR, n1, between);

  vec3 redRefract = refract(eye, nor, redIORRatio);
  vec3 yellowRefract = refract(eye, nor, yellowIORRatio);
  vec3 greenRefract = refract(eye, nor, greenIORRatio);
  vec3 cyanRefract = refract(eye, nor, cyanIORRatio);
  vec3 blueRefract = refract(eye, nor, blueIORRatio);
  vec3 purpleRefract = refract(eye, nor, purpleIORRatio);

  float r = scene(redRefract, redIORRatio).r * 0.5;
  float y = dot(scene(yellowRefract, yellowIORRatio), vec3(2.0, 2.0, -1.0)) / 6.0;
  float g = scene(greenRefract, greenIORRatio).g * 0.5;
  float c = dot(scene(cyanRefract, cyanIORRatio), vec3(-1.0, 2.0, 2.0)) / 6.0;
  float b = scene(blueRefract, blueIORRatio).b * 0.5;
  float p = dot(scene(purpleRefract, purpleIORRatio), vec3(2.0, -1.0, 2.0)) / 6.0;

  float R = r + (2.0*p + 2.0*y - c)/3.0;
  float G = g + (2.0*y + 2.0*c - p)/3.0;
  float B = b + (2.0*c + 2.0*p - y)/3.0;

  #else

  #ifdef HUE
  const float hueStep = 1.0 / float(HUE_NUM);
  vec3 color = vec3(0.);

  for (int i = 0; i < HUE_NUM; i++) {
    float hue = float(i) * hueStep;
    float ior = hue2IOR(360.0 * hue, greenIOR, n1, between);
    vec3 iorRefract = refract(eye, nor, ior);

    #ifdef COS_HUE
    color += (0.5 + 0.5 * cos(TWO_PI * (hue + vec3(0, 0.33, 0.67)))) * scene(iorRefract, ior);

    // Red / Teal
    // color += (vec3(0.8, 0.5, 0.4) + vec3(0.2, 0.4, 0.2) * cos(TWO_PI * (vec3(2, 1, 0) * hue + vec3(0, 0.25, 0.25)))) * scene(iorRefract, ior);

    // Simple
    // color += (0.5 + 0.5 * cos(TWO_PI * (hue + vec3(0, 0.1, 0.2)))) * scene(iorRefract, ior);

    // Simple 2
    // color += (0.5 + 0.5 * cos(TWO_PI * (vec3(1.5) * hue + vec3(0, 0.1, 0.2)))) * scene(iorRefract, ior);

    // Neon
    // color += (0.5 + 0.5 * cos(TWO_PI * (vec3(2, 1, 0) * hue + vec3(0.5, 0.2, 0.25)))) * scene(iorRefract, ior);

    // Something 1
    // color += (vec3(0.5, 0.3, 0.5) + vec3(0.25, 0.5, 0.7) * cos(TWO_PI * (vec3(2, 1, 0.5) * hue + vec3(0.6, 0.1, 0.25)))) * scene(iorRefract, ior);

    // RED CYAN
    // color += mix(#FF0000, #00FFFF, hue) * scene(iorRefract, ior);

    #else
    // color += hsv(vec3(hue, 1.0, 1.0)) * scene(iorRefract, ior);
    const vec3 magenta = pow(#FF1799, vec3(2.2));
    const vec3 purple = pow(#9636FF, vec3(2.2));

    float dI = dot(nor, -eye);

    vec3 sceneResult = scene(iorRefract, ior);
    vec3 mixI = clamp(0.5 + 0.5 * sin(4.0 * dI + sin(nor)), 0.0, 1.0);

    vec3 thisColor = vec3(0);
    thisColor += mix(#FF0000, #00FFFF, mixI)
      * sceneResult;

    color += thisColor * (0.75 + 0.25 * cos(TWO_PI * (dI + vec3(0., 0.33, 0.67))));
    #endif
  }

  color *= pow(hueStep, 0.8);

  float R = color.r;
  float G = color.g;
  float B = color.b;

  #else

  #ifdef REFR_INTEGRAL
  float r1 = hue2IOR(0.0, greenIOR, n1, between);
  float r2 = hue2IOR(240.0, greenIOR, n1, between);
  // float r1 = (n2 - amount) / n1;
  // float r2 = (n2 + amount) / n1;
  vec3 intRefracted = intRefract(eye, nor, r1, r2);
  // intRefracted = 0.5 + 0.5 * cos(TWO_PI * (intRefracted + vec3(0, 0.33, 0.67)));
  // intRefracted = 0.5 + 0.5 * cos(intRefracted);
  vec3 color = scene(intRefracted, r1);
  // vec3 color = intRefracted;

  float R = color.r;
  float G = color.g;
  float B = color.b;

  #else

  float redIORRatio = hue2IOR(0.0, greenIOR, n1, between);
  float greenIORRatio = hue2IOR(120.0, greenIOR, n1, between);
  float blueIORRatio = hue2IOR(240.0, greenIOR, n1, between);

  vec3 redRefract = refract(eye, nor, redIORRatio);
  vec3 greenRefract = refract(eye, nor, greenIORRatio);
  vec3 blueRefract = refract(eye, nor, blueIORRatio);

  float r = scene(redRefract, redIORRatio).r;
  float g = scene(greenRefract, greenIORRatio).g;
  float b = scene(blueRefract, blueIORRatio).b;

  float R = r;
  float G = g;
  float B = b;
  #endif
  #endif
  #endif

  return vec3(R, G, B) * lightColor;
}
vec3 refractColors (in vec3 nor, in vec3 eye, in float n2) {
  return refractColors(nor, eye, n2, 1., vec3(1.0));
}
vec3 refractColors (in vec3 nor, in vec3 eye, in float n2, in float n1) {
  return refractColors(nor, eye, n2, n1, vec3(1.0));
}
vec3 refractColors (in vec3 nor, in vec3 eye, in float n2, in vec3 lightColor) {
  return refractColors(nor, eye, n2, 1., lightColor);
}

#pragma glslify: export(refractColors)
