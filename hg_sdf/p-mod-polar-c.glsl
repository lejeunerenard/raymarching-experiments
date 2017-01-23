// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolarC(inout vec2 p, float repetitions) {
  float angle = 2.*PI/repetitions;
  float a = atan(p.y, p.x) + angle/2.;
  float r = length(p);
  float c = floor(a/angle);
  a = mod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
  // For an odd number of repetitions, fix cell index of the cell in -x direction
  // (cell index would be e.g. -5 and 5 in the two halves of the cell):
  if (abs(c) >= (repetitions/2.)) c = abs(c);
  return c;
}
void pModPolar(inout vec2 p, float repetitions) {
  float angle = 2.*PI/repetitions;
  float a = atan(p.y, p.x) + angle/2.;
  float r = length(p);
  a = mod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
}
