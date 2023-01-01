#extension GL_OES_standard_derivatives : enable
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define PHI (1.618033988749895)
#define saturate(x) clamp(x, 0.0, 1.0)

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float scale;
uniform vec3 cameraRo;

// --- SVG ---
float svgD = 1e38;
float _b;

#define N 20.  // splines discretization. Lower it on slow GPUs
// absolute main SVG commands
#define M(x,y)             x0 = _x = x;   y0 = _y = y;
#define L(x,y)             svgD = smin(svgD, line(uv, vec2(_x,_y), vec2(x,y)) ); _x=x,_y=y;
#define C(x1,y1,x2,y2,x,y) svgD = smin(svgD, bezier(uv, vec2(_x,_y), vec2(x1,y1),vec2(x2,y2), vec2(x,y)) ); _x=x,_y=y; 
#define H(x)               svgD = smin(svgD, line(uv, vec2(_x,_y), vec2(x,_y)) ); _x=x;
#define V(y)               svgD = smin(svgD, line(uv, vec2(_x,_y), vec2(_x,y)) ); _y=y;
#define Z                  svgD = smin(svgD, line(uv, vec2(_x,_y), vec2(x0,y0)) ); _x=x0,_y=y0;
// relative main SVG commands
#define m(x,y)             M(_x+x,_y+y)
#define l(x,y)             L(_x+x,_y+y)
#define c(x1,y1,x2,y2,x,y) C(_x+x1,_y+y1,_x+x2,_y+y2,_x+x,_y+y)
#define h(x)               H(_x+x)
#define v(y)               V(_y+y)
#define svgz               Z

#define style(f,c,d)       fill=f; S=d; COL= mod(vec4((c)/65536,(c)/256,c,1),256.)/255.;
#define path(cmd)          svgD = 1e38; cmd; if (fill>0.) svgz; draw(svgD,O);

#pragma glslify: import(./svg.glsl)

vec3 two_dimensional (in vec2 uv) {
  vec3 color = vec3(0);

  vec2 q = uv; // UV is centered at (0, 0)

  vec2 wQ = q;

  wQ.xy += cameraRo.xy;
  wQ *= cameraRo.z;

  wQ.y *= -1.; // SVG Coords are upsidedown

  vec2 totalScale = vec2(1.);

  // Adjust into pixel space : [-resolution, resolution]
  // Scale is the ratio of the FBO resolution divided by SVG resolution
  totalScale *= resolution / scale;
  // This Account for it being [-1, 1] & resolution being [0, 1]. So makes the
  // space [-resolution/2, resolution/2]
  totalScale *= 0.5;

  wQ *= totalScale;

  q = wQ;

  float logo = SVG(q);
  logo /= min(totalScale.x, totalScale.y);

  color = vec3(logo);
  // color = vec3(step(0.5, logo));

  // // Debug space
  // color.rg = q;
  // color.rg = color.gr; // Look at y

  // Convert to [0, 1] so it fits in a texture pixel value
  color = 0.5 * (color + 1.0);

  // // Show max clip
  // vec3 clipColor = color;
  // if (color.x >= 1.) {
  //   clipColor = vec3(1, 0, 1);
  // }
  // // "Zero" before [0, 1] scale
  // if (color.x < 0.5) {
  //   clipColor = vec3(0, 0, 1);
  // }
  // // Sub zero
  // if (color.x < 0.) {
  //   clipColor = vec3(0, 1, 0);
  // }
  // color = clipColor;

  return color.rgb;
}

void main() {
    vec2 uv = fragCoord.xy;
    gl_FragColor = saturate(vec4(two_dimensional(uv), 1));
}
