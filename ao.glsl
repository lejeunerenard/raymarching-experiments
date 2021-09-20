float calcAO( in vec3 pos, in vec3 nor, in float t ) {
  float occ = 0.0;
  float sca = 1.0;
  for( int i=0; i<4; i++  ) {
    float hr = 0.01 + 0.12*float(i)/4.0;
    vec3 aopos =  nor * hr + pos;
    vec3 dd = map(aopos, t);
    occ += -(dd.x-hr)*sca;
    sca *= 0.95;
  }
  return clamp( 1.0 - 3.0*occ, 0.0, 1.0  );
}

float calcAO( in vec3 pos, in vec3 nor  ) {
  return calcAO(pos, nor, 0.);
}

#pragma glslify: export(calcAO)
