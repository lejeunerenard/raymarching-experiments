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

#pragma glslify: snoise3 = require(glsl-noise/simplex/3d)

float lumPeriod (vec3 p, int i) {
  return .2 * (snoise3(p * vec3(1., .3, 1.)) + 2. * float(i) + 1.);
}

#pragma glslify: iridescant = require(./iridescent, lumPeriod=lumPeriod)

const float epsilon = 0.001;
const int maxSteps = 40;
const float maxDistance = 15.;
const vec3 background = vec3(.15);

const vec3 lightPos = vec3(0, 0, 5.);

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

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

const float tOff = 168.;
float rtime1 = (time + tOff) * .12;
float rtime2 = (time + tOff) * .27;
float rtime3 = (time + tOff) * .13;
mat3 rot =
  mat3(cos(rtime1),0,sin(rtime1),0,1,0,-sin(rtime1),0,cos(rtime1)) *
  mat3(cos(rtime2),sin(rtime2),.0,-sin(rtime2),cos(rtime2),.0,0,0,1) *
  mat3(1,0,0,0,cos(rtime3),sin(rtime3),0,-sin(rtime3),cos(rtime3));

void rotFold (inout vec3 p) {
  float e=.4;
  for (int i=0; i < 9; i++) {
    p = abs(p*rot) - e;
    p.y -= p.x*.1;
    p.x -= p.z*.1;
    e = e * .8 + e *e * .1;
  }
  p = abs(p * rot) - e;
  p = abs(p * rot) - e;
}

vec2 map (in vec3 p) {
  vec4 pp = vec4(p, 1);
  vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), 2. * PI * sin(time) / 2.) * pp).xyz;

  // Space Transforms
  rotFold(q);

  vec2 d = vec2(sdBox(q, vec3(.05)), 1.);

  return d;
}

vec2 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t  = 1.5;
  float maxI = 0.;

  for (int i = 0; i < maxSteps; i++) {
    vec2 d = map(rayOrigin + rayDirection * t);
    if (d.x < epsilon) return vec2(t + d.x, d.y);
    t += d.x;
    maxI = float(i);
    if (t > maxDistance) break;
  }
  return vec2(-1., maxI);
}

vec3 getNormal( in vec3 p, in float eps ) {
    vec2 e = vec2(1.0,-1.0)*0.05773*eps;
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
    for( int i=0; i<16; i++ )
    {
    vec2 h = map(ro + rd*t);
        res = min( res, 4.0*h.x/t );
        t += clamp( h.x, 0.02, 0.10 );
        if( h.x<0.001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

float calcAO( in vec3 pos, in vec3 nor  ) {
  float occ = 0.0;
  float sca = 1.0;
  for( int i=0; i<5; i++  ) {
    float hr = 0.01 + 0.12*float(i)/4.0;
    vec3 aopos =  nor * hr + pos;
    vec2 dd = map(aopos);
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

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec2 t ) {
    vec3 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, 0.001 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = vec3(1.);
      // color.r *= dot(nor, ref);
      // color.g *= max(nor.x, t.y) / 6.;

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      //dif *= min(0.3 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
      vec3 lin = vec3(0.);
      lin += 1.1*dif*vec3(1.);
      lin += 0.20*amb*vec3(0.50,0.70,1.00)*occ;
      color *= lin;

      // Fog
      color = mix(background, color, (maxDistance-t.x) / maxDistance);

      color *= (dot(rayDirection,nor)*.5+.5) * 3.;
      color += matCap(ref);

      // Color Map
      color = mix(vec3(0., .7, 1.), color, clamp(exp(length(vec4(color, 1.))) * .325, 0., 1.));
      color = mix(vec3(1., .3, 0.), color, 1. - length(vec4(color, 1.)) *.05);

      #ifdef debugMapMaxed
      if (t.y / float(maxSteps) > 0.9) {
        color = vec3(1., 0., 1.);
      }
      #endif

      #ifdef debugMapCalls
      color = vec3(t.y / float(maxSteps));
      #endif
    } else {
      // Radial Gradient
      color *= sqrt(1.75 - 1.25 * dot(rayDirection, vec3(0., 0., -1.)));
      // Glow
      // color = mix(vec3(pos * .7), color, 1. - clamp(t.y / float(maxSteps), 0., 1.));
    }

    return vec4(color, 1.);
}

// Original

void main() {
    const float d = 7.;

    vec3 ro = vec3(0.,0.,d) + cOffset;

    #ifdef SS
    // Antialias by averaging all adjacent values
    vec4 color = vec4(0.);
    vec2 t = vec2(0.);
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
    vec2 t = march(ro, rd);
    gl_FragColor = shade(ro, rd, t);
    #endif
}
