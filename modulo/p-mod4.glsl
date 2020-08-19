// Repeat in 4 dimensions
vec4 pMod4(inout vec4 p, vec4 size) {
  vec4 c = floor((p + size*0.5)/size);
  p = mod(p + size*0.5,size) - size*0.5;
  return c;
}
#pragma glslify: export(pMod4)
