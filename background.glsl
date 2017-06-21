vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  // DARK Grey
  return mix(vec3(0.0, 0.0, 0.0), vec3(0.02, 0.01, 0.01), coord);

  return pow(#D1CC9B, vec3(2.2));

  // return pow(mix(#FFBABE, #D5A9CA, coord), vec3(2.2));
  // return pow(mix(#ff0000, #0000ff, coord), vec3(2.2));
  vec3 color =  mix(#ff7777, #aa77ff, coord);

  if (length(color) >= 1.732051) {
    color = #00ff00;
  } else if (length(color) <= 0.005) {
    color = #ff0000;
  }

  return color;
}
vec3 background = vec3(0.);
