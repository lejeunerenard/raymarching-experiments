mat2 rotationMatrix2 (in float a ) {
  float c = cos(a);
  float s = sin(a);
  return mat2(c, -s, s, c);
}

#pragma glslify: export(rotationMatrix2)
