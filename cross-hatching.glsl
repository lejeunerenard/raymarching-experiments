// source:
// http://learningwebgl.com/blog/?p=2858
float crossHatching (in vec3 color) {
  float v = 1.0;
  float l = color.x;

  float period = 0.021;
  float limit = period * 0.2;
  if (mod(dot(fragCoord.xy, vec2(1)), period) <= limit && l > 1.0) {
    v = 0.0;
  }
  if (mod(dot(fragCoord.xy, vec2(1, -1)), period) <= limit && l > 0.75) {
    v = 0.0;
  }
  if (mod(dot(fragCoord.xy, vec2(1)) - 5.0, period) <= limit && l > 0.5) {
    v = 0.0;
  }
  if (mod(dot(fragCoord.xy, vec2(1, -1)) - 5.0, period) <= limit && l > 0.3465) {
    v = 0.0;
  }

  return 1.0 - v;
}

#pragma glslify: export(crossHatching)
