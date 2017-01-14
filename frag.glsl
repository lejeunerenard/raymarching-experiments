#define PI 3.1415926536

// #define debugMapCalls
// #define debugMapMaxed
#define SS 2

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform vec3 cOffset;
uniform mat4 orientation;
uniform mat4 projectionMatrix;
uniform sampler2D texture;

#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)

float epsilon = 0.01;
const int maxSteps = 128;
float maxDistance = 19.;
vec3 background = #222222;

vec3 lightPos = vec3(0, 0, 5.);

const vec2 un = vec2(1., 0.);

vec2 dMin (vec2 d1, vec2 d2) {
  return (d1.x < d2.x) ? d1 : d2;
}
mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

void fold (inout vec2 p) {
  if (p.x + p.y < 0.) {
    float x1 = -p.y;
    p.y = -p.x;
    p.x = x1;
  }
}
void foldInv (inout vec2 p) {
  if (p.x - p.y < 0.) {
    float x1 = p.y;
    p.y = p.x;
    p.x = x1;
  }
}

mat3 rotation3 (float time, float tOff) {
  float rtime1 = PI + 1. * PI * .23 * sin((time + tOff) * .12);
  float rtime2 = PI + 1. * PI * .13 * sin((time + tOff + 1.) * .93);
  float rtime3 = PI + 1. * PI * .23 * sin((time + tOff) * .71);
  return
    mat3(cos(rtime1),0,sin(rtime1),0,1,0,-sin(rtime1),0,cos(rtime1)) *
    mat3(cos(rtime2),sin(rtime2),.0,-sin(rtime2),cos(rtime2),.0,0,0,1) *
    mat3(1,0,0,0,cos(rtime3),sin(rtime3),0,-sin(rtime3),cos(rtime3));
}

vec2 kifs( inout vec3 p ) {
  float r = dot(p, p);
  const float scale = 2.1;

  int maxI = 0;
  const vec3 Offset = vec3(2., -1., -.2);

  const int Iterations = 10;

  //mat4 rot = rotationMatrix(vec3(0.5,-0.05,-0.5), 10. * sin(time));
  mat3 rot = rotation3(10. * sin(time), 0.);

  float t = 1.5+sin(0.*.5)*.5;
  p.xy *= mat2(cos(t), sin(t), -sin(t), cos(t));

  float trap = maxDistance;

  for (int i = 0; i < Iterations; i++) {
    p.yz = abs(p.yz);

    // Folding
    fold(p.xy);
    fold(p.xz);
    foldInv(p.xy);
    foldInv(p.xz);

    //p = (vec4(p, 1.) * rot).xyz;
    p *= rot;

    // Stretch
    p *= scale;
    p -= Offset * (scale - 1.);

    trap = min(trap, length(p));
  }

  return vec2((length(p) - .1) * pow(scale, - float(Iterations)), trap);
}

vec3 map (in vec3 p) {
  vec4 pp = vec4(p, 1);
  vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), 2. * PI * sin(20. * time) / 4.) * pp).xyz;
  // vec3 q = vec3(orientation * pp).xyz;

  vec2 fractal = kifs(q);
  return vec3(fractal.x, 1., fractal.y);
}

vec4 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  for (int i = 0; i < maxSteps; i++) {
    vec3 d = map(rayOrigin + rayDirection * t);
    if (d.x < epsilon) return vec4(t + d.x, d.y, float(i), d.z);
    t += d.x;
    maxI = float(i);
    trap = d.z;
    if (t > maxDistance) break;
  }
  return vec4(-1., 0., maxI, trap);
}

vec3 getNormal( in vec3 p, in float eps ) {
    vec2 e = vec2(1.0,-1.0)*.05773*eps;
    return normalize(
      e.xyy * map( p + e.xyy ).x + 
      e.yyx * map( p + e.yyx ).x + 
      e.yxy * map( p + e.yxy ).x + 
      e.xxx * map( p + e.xxx ).x );
}

// Material Functions
float diffuse (in vec3 nor, in vec3 lightPos) {
  return clamp(dot(nor, lightPos) / length(lightPos), 0., 1.);
}

// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax ) {
  float res = 1.0;
    float t = mint;
    for( int i=0; i<16; i++ ) {
      vec3 h = map(ro + rd*t);
      res = min( res, 4.0*h.x/t );
      t += clamp( h.x, 0.02, 0.10 );
      if( h.x<0.001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

float calcAO( in vec3 pos, in vec3 nor  ) {
  float occ = 0.0;
  float sca = 1.0;
  for( int i=0; i<4; i++  ) {
    float hr = 0.01 + 0.12*float(i)/4.0;
    vec3 aopos =  nor * hr + pos;
    vec3 dd = map(aopos);
    occ += -(dd.x-hr)*sca;
    sca *= 0.95;
  }
  return clamp( 1.0 - 3.0*occ, 0.0, 1.0  );    
}

vec3 matCap (vec3 ref) {
  float m = 2. * sqrt( pow( ref.x, 2. ) + pow( ref.y, 2. ) + pow( ref.z + 1., 2. ) );
  vec2 vN = ref.xy / m + .5;
  return texture2D(texture, vN).rgb;
}

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#90B8FF, color, 1. - l * .0625);
  // Dark
  color = mix(#7782E8, color, clamp(exp(l) * .325, 0., 1.));
}

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t ) {
    vec3 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, 0.001 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = length(pos) * #FF8055;

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      dif *= min(0.3 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      vec3 lin = vec3(0.);
      lin += dif;
      lin += 0.2*amb*vec3(0.50,0.70,1.00)*occ;
      color *= lin;

      // Fog
      color = mix(background, color, (maxDistance-t.x) / maxDistance);

      // Inner Glow
      vec3 glowColor = vec3(1.0, 0.075, 0.01) * 5.0;
      float fGlow = clamp(t.w * 0.1, 0.0, 1.0);
      fGlow = pow(fGlow, 3.0);
      color += glowColor * 15.0 * fGlow;

      color *= exp(-t.x * .4);

      colorMap(color);

      #ifdef debugMapMaxed
      if (t.z / float(maxSteps) > 0.9) {
        color = vec3(1., 0., 1.);
      }
      #endif

      #ifdef debugMapCalls
      color = vec3(t.z / float(maxSteps));
      #endif
    } else {
      // Radial Gradient
      color *= sqrt(1.75 - 1.25 * dot(rayDirection, vec3(0., 0., -1.)));
      // Glow
      color = mix(vec3(1.), color, 1. - .3 * clamp(t.z / float(maxSteps), 0., 1.));
    }

    return vec4(color, 1.);
}

// Original

void main() {
    const float d = 3.;

    vec3 ro = vec3(0.,0.,d) + cOffset;

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec4 t = vec4(0.);
    for (int x = - SS / 2; x < SS / 2; x++) {
        for (int y = - SS / 2; y < SS / 2; y++) {
            vec3 rd = getRayDirection(vec2(
                  float(x) / resolution.y + fragCoord.x,
                  float(y) / resolution.y + fragCoord.y),
                  projectionMatrix);
            t = march(ro, rd);
            color += shade(ro, rd, t);
        }
    }
    gl_FragColor = color / float(SS * SS);

    #else
    vec3 rd = getRayDirection(fragCoord, projectionMatrix);
    vec4 t = march(ro, rd);
    gl_FragColor = shade(ro, rd, t);
    #endif
}
