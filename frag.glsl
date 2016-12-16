#define PI 3.1415926536

// #define debugMapCalls
// #define debugMapMaxed
// #define SS 2

precision highp float;

#pragma glslify: getRayDirection = require(./ray-apply-proj-matrix)

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform vec3 cOffset;
uniform mat4 orientation;
uniform mat4 projectionMatrix;

const float epsilon = 0.003;
const int maxSteps = 512;
const float maxDistance = 20.;
const vec4 background = vec4(vec3(.6), 1.);

const vec3 lightPos = vec3(0, 0, 5.);

float sdSphere (in vec3 p, float radius) {
  return length(p) - radius;
}

void opReflect (inout vec3 p, in vec4 plane) {
  vec3 normal = normalize(plane.xyz);
  vec3 v = p - plane.w * normal;
  float d = dot(v, normal);
  vec3 reflected = p - 2. * d * normal;

  // Distance to plane
  float t = (plane.w - dot(p, normal)) / dot(normal, normal);

  float which = step(0., t);
  p = which * p + (1. - which) * reflected;
}
vec2 rotate(vec2 p, float ang) {
    float c = cos(ang), s = sin(ang);
    return vec2(p.x*c - p.y*s, p.x*s + p.y*c);
}

vec3 repeatAngS(vec2 p, float n) {
    float ang = 2.0*PI/n;
    float sector = floor(atan(p.x, p.y)/ang + 0.5);
    p = rotate(p, sector*ang);
    return vec3(p.x, p.y, mod(sector, n));
}
vec2 repeatAng(vec2 p, float n) {
    float ang = 2.0*PI/n;
    float sector = floor(atan(p.x, p.y)/ang + 0.5);
    p = rotate(p, sector*ang);
    return vec2(p.x, p.y);
}

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

const vec2 un = vec2(1., 0.);

float obj1 (in vec3 p) {
  vec3 bSize = vec3(.25);
  float box1 = sdBox(p, bSize);
  float sph1 = sdSphere(p, .35);
  opReflect(p, un.xyyy);
  opReflect(p, un.xxyy);

  p.x += bSize.x;
  float sph2 = sdSphere(p, .1);

  return min(
    max(box1, sph1),
    sph2
  );
}

float map (in vec3 p) {
  float t = time;

  vec4 pp = vec4(p, 1);
  vec3 q = vec3(orientation * pp).xyz;

  q.xy = repeatAng(q.xy, 5.0);
  opReflect(q, un.xyyy);

  q.zx = repeatAng(q.zx, 5.0);
  opReflect(q, un.yyxy);

  float org = sdBox(q, vec3(.15));

  vec3 scale = sin(vec3(1., 2., 3.) * time);
  vec3 obj1P = scale * vec3(.5, .5, 0.) + q;
  float d = min(sdSphere(obj1P, .25), org);

  const vec3 bSize = vec3(.25);

  for (int i = 0; i < 8; i++) {
    d = min(sdBox(q + sin(time + float(i)*vec3(1., 2., 3.)), bSize), d);
  }

  for (int i = 0; i < 8; i++) {
    d = min(sdSphere(q + sin(time + float(i)*vec3(5., 7., 11.)), float(i) / 16.), d);
  }
  return d;
}

vec2 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t  = 1.5;
  float maxI = 0.;

  for (int i = 0; i < maxSteps; i++) {
    float d = map(rayOrigin + rayDirection * t);
    if (d < epsilon) return vec2(t + d, i);
    t += d;
    if (t > maxDistance) break;
  }
  return vec2(-1., maxI);
}

vec3 getNormal( in vec3 p, in float eps ) {
    vec2 e = vec2(1.0,-1.0)*0.05773*eps;
    return normalize( e.xyy*map( p + e.xyy ) + 
            e.yyx*map( p + e.yyx ) + 
            e.yxy*map( p + e.yxy ) + 
            e.xxx*map( p + e.xxx ) );
}

float diffuse (in vec3 nor, in vec3 lightPos) {
  return clamp(dot(nor, lightPos) / length(lightPos), 0., 1.);
}
// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax ) {
  float res = 1.0;
    float t = mint;
    for( int i=0; i<16; i++ )
    {
    float h = map(ro + rd*t);
        res = min( res, 4.0*h/t );
        t += clamp( h, 0.02, 0.10 );
        if( h<0.001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}
float calcAO( in vec3 pos, in vec3 nor  ) {
  float occ = 0.0;
  float sca = 1.0;
  for( int i=0; i<5; i++  ) {
    float hr = 0.01 + 0.12*float(i)/4.0;
    vec3 aopos =  nor * hr + pos;
    float dd = map(aopos);
    occ += -(dd-hr)*sca;
    sca *= 0.95;
  }
  return clamp( 1.0 - 3.0*occ, 0.0, 1.0  );    
}

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec2 t ) {
    vec4 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, 0.001 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      color = vec4(1.);
      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      dif *= softshadow(pos, lightPos, 0.02, 1.5);
      vec3 lin = vec3(0.);
      lin += 1.1*dif*vec3(1.);
      lin += 0.20*amb*vec3(0.50,0.70,1.00)*occ;
      lin += 0.20*(t.y / float(maxSteps))*vec3(0.30,0.90,1.00);
      color.xyz *= lin;
      color.a = 1.;

      // Fog
      color = mix(background, color, (maxDistance-t.x) / maxDistance);

      // Color Map
      color  = mix(vec4(0., .7, 1., 1.), color, clamp(exp(length(color)) * .325, 0., 1.));
      color  = mix(vec4(1., .3, 0., 1.), color, 1. - length(color) *.05);

      #ifdef debugMapMaxed
      if (t.y / float(maxSteps) > 0.9) {
        color = vec4(1., 0., 1., 1.);
      }
      #endif

      #ifdef debugMapCalls
      color = vec4(vec3(t.y / float(maxSteps)), 1.);
      #endif
    } else {
      color *= sqrt(1.45 - .5 * dot(rayDirection, vec3(0., 0., -1.)) * 2.5);
      color = mix(vec4(1.), color, 1. - clamp(t.y / float(maxSteps), 0., 1.));
      color.a = 1.;
    }
    return color;
}

// Original

void main() {
    const float d = 3.;

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
