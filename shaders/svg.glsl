// Source:
// https://www.shadertoy.com/view/ldXyRn

// improved version of https://www.shadertoy.com/view/MlVSWc
// === SVG Player ====      short spec: https://www.w3.org/TR/2008/REC-SVGTiny12-20081222/paths.html

float bezier(vec2,vec2,vec2,vec2,vec2);
float line(vec2,vec2,vec2);
void  draw(float,inout vec4);
const float FILL=1., CONTOUR=0., DIR=1., INDIR=-1.;
vec4 COL = vec4(0); float fill=FILL, S=DIR;  // style state

float dh2=1e38, dh1;

// smin: min of |dist| keeping signed distance.
// smin0() is for easy cases (i.e. along splines)
// smin() deals with ambiguous distance sign in the exterior fan for angles > 90°
#define smin0(a,b)  abs(_b=b) < abs(a) ? _b : a  
#define smin(a,b)   abs(_b=b) < abs(a)+1e-2 ? abs(abs(_b)-abs(a))<1e-2 ? dh1<dh2 ? _b:a: _b : a

// --- spline interpolation ( inspired from revers https://www.shadertoy.com/view/MlGSz3 )
vec2 interpolateSVG (vec2 G1, vec2 G2, vec2 G3, vec2 G4, float t) {
  vec2 A = G4-G1 + 3.*(G2-G3),
       B = 3.*(G1-2.*G2+G3),
       C = 3.*(G2-G1),
       D = G1;
  return t * (t * (t * A + B) + C) + D;
}

float line (vec2 p, vec2 a, vec2 b) {
  vec2 pa = p - a, ba = b - a;
  float l = dot(ba, ba),
        h0 = dot(pa, ba) / l,                    // parameterization of proj on line
        h = clamp(h0 , 0., 1.);                  // parameterization of proj on segment
  if (l<1e-6) return 1e38;                       // degenerated line
  dh2 = dh1, dh1 = abs(h-h0)*length(ba);         // how much beyond (for sign desambiguation)
  vec2 d = pa - ba * h;                          // distance vector to segment
  if  ( (a.y>p.y) != (b.y>p.y) &&
      pa.x < ba.x * pa.y / ba.y ) S = -S;        // track interior vs exterior
  return dot(d,d);                              // optimization by deferring sqrt
}

float bezier ( vec2 uv, vec2 A, vec2 B, vec2 C, vec2 D) {
  //float svgD = 1e5;                               // for correct sign desambiguation with prev
  vec2 p = A;
  for (float t = 1.; t <= N; t++) {
    vec2 q = interpolateSVG(A, B, C, D, t/N);
    float l = line(uv, p, q);
    svgD = t==1. ? smin (svgD, l ) : smin0(svgD, l );   // t==1 :  same thing
    p = q;
  }
  return svgD;
}

void draw(float d, inout float O) 
{
  d = sqrt(d);  // optimization by deferring sqrt here

  d = fill>0. ? S * d : d; // Apply sidedness
  O = min(O, d);
}

// === SVG drawing ===============================================================
float SVG (vec2 uv) {
  float O = 1e10;
  float _x, _y, x0, y0;

  // Copy-paste your SVG pathes here.  Slight adaptations :
  //  - add () around command params and  comma between points,
  //  - split polylines and polybéziers into sets of 1 vs 3 pairs of coordinates
  //  - path( style( FILL/CONTOUR, color(hexa) )
  //          commands 
  //        )

#define SVG_PATHS 1

  return O;
}
// #pragma glslify: export(mainImage)
