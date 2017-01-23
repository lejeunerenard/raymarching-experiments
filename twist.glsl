vec3 opTwist( vec3 p, float angle ) {
    float c = cos(angle);
    float s = sin(angle);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}

#pragma glslify: export(opTwist)
