#ifndef PI
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#endif

#define RGBCMY 1
// #define REFR_INTEGRAL 1
// #define HUE 1
// #pragma glslify: hsv = require(glsl-hsv2rgb)

// #pragma glslify: hue2IOR = require(./dispersion-ray-direction)
// #pragma glslify: hue2IOR = require(./dispersion/hue-to-ior-exponential)
#pragma glslify: hue2IOR = require(./dispersion/hue-to-ior-sigmoid)
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
  const float between = amount;
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
  float ior1 = hue2IOR(0.0, greenIOR, n1, between);
  float ior2 = hue2IOR(90.0, greenIOR, n1, between);
  float ior3 = hue2IOR(180.0, greenIOR, n1, between);
  float ior4 = hue2IOR(240.0, greenIOR, n1, between);

  vec3 ior1Refract = refract(eye, nor, ior1);
  vec3 ior2Refract = refract(eye, nor, ior2);
  vec3 ior3Refract = refract(eye, nor, ior3);
  vec3 ior4Refract = refract(eye, nor, ior4);

  vec3 color = vec3(0.);
  color += hsv(vec3(0.0, 1.0, 1.0)) * scene(ior1Refract, ior1);
  color += hsv(vec3(90.0, 1.0, 1.0)) * scene(ior2Refract, ior2);
  color += hsv(vec3(180.0, 1.0, 1.0)) * scene(ior3Refract, ior3);
  color += hsv(vec3(240.0, 1.0, 1.0)) * scene(ior4Refract, ior4);

  color *= 0.25;

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
