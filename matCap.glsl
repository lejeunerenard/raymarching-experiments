vec3 matCap (vec3 ref) {
  float m = 2. * sqrt( pow( ref.x, 2. ) + pow( ref.y, 2. ) + pow( ref.z + 1., 2. ) );
  vec2 vN = ref.xy / m + .5;
  return texture2D(texture, vN).rgb;
}
#pragma glslify: export(matCap)
