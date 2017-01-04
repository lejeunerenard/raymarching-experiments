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
uniform sampler2D texture;

const float epsilon = 0.003;
const int maxSteps = 512;
const float maxDistance = 20.;
const vec4 background = vec4(vec3(1.), 1.);

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

vec2 map (in vec3 p) {
  float t = time;

  vec4 pp = vec4(p, 1);
  vec3 q = vec3(orientation * rotationMatrix(vec3(0., 1. ,0.), 2. * PI * sin(time) / 2.) * pp).xyz;

  q.xy = repeatAng(q.xy, 5.0);
  opReflect(q, un.xyyy);

  q.zx = repeatAng(q.zx, 5.0);
  opReflect(q, un.yyxy);

  vec3 scale = sin(vec3(1., 2., 3.) * time);
  vec3 obj1P = scale * vec3(.5, .5, 0.) + q;
  vec2 d = vec2(min(sdSphere(obj1P, .25), 1000.), 0.);

  const vec3 bSize = vec3(.25);

  d = dMin(vec2(sdBox(q - vec3(.5) - .5 * sin(time + vec3(1., 2., 3.)), bSize), 1.), d);
  d = dMin(vec2(sdBox(q - vec3(.5) - .5 * sin(time + 2. * vec3(1., 2., 3.)), bSize), 2.), d);
  d = dMin(vec2(sdBox(q - vec3(.5) - .5 * sin(time + 3. * vec3(3., 5., 7.)), bSize), 3.), d);
  d = dMin(vec2(sdBox(q - vec3(.5) - .5 * sin(time + 4. * vec3(11., 13., 17.)), bSize), 4.), d);
  d = dMin(vec2(sdBox(q - vec3(.5) - .5 * sin(time + 5. * vec3(3., 5., 3.)), bSize), 5.), d);

  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + vec3(5., 7., 11.)), 4. / 16.), 6.), d);
  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + 2. * vec3(5., 7., 11.)), 2. / 16.), 7.), d);
  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + 3. * vec3(5., 7., 11.)), 3. / 16.), 8.), d);
  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + 4. * vec3(5., 7., 11.)), 4. / 16.), 9.), d);
  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + 5. * vec3(5., 7., 11.)), 5. / 16.), 10.), d);
  d = dMin(vec2(sdSphere(q + vec3(.2) + sin(time + 6. * vec3(5., 7., 11.)), 6. / 16.), 11.), d);

  return d;
}

vec2 march (in vec3 rayOrigin, in vec3 rayDirection) {
  float t  = 1.5;
  float maxI = 0.;

  for (int i = 0; i < maxSteps; i++) {
    vec2 d = map(rayOrigin + rayDirection * t);
    if (d.x < epsilon) return vec2(t + d.x, d.y);
    t += d.x;
    if (t > maxDistance) break;
  }
  return vec2(-1., maxI);
}

vec3 getNormal( in vec3 p, in float eps ) {
    vec2 e = vec2(1.0,-1.0)*0.05773*eps;
    return normalize( e.xyy*map( p + e.xyy ).x + 
            e.yyx*map( p + e.yyx ).x + 
            e.yxy*map( p + e.yxy ).x + 
            e.xxx*map( p + e.xxx ).x );
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

vec4 shade( in vec3 rayOrigin, in vec3 rayDirection, in vec2 t ) {
    vec4 color = background;
    vec3 pos = rayOrigin + rayDirection * t.x;
    if (t.x>0.) {
      vec3 nor = getNormal(pos, 0.001 * t.x);
      vec3 ref = reflect(rayDirection, nor);

      // Basic Diffusion
      vec3 a = vec3(.25, .75, .75);
      vec3 b = vec3(.12, .25, .25);
      vec3 c = vec3(1., .9, .5);
      vec3 d = vec3(0., .33, .67);
      if (t.y > 6.){
        color = texture2D(texture, nor.xy);
      } else {
        color = vec4(1.);
      }
      color.rgb = max(color.rgb, vec3(.4));
      color.rgb *= a + b * cos(2. * PI * ( c * sin(time / 100.) * t.y / 11. + d ));
      color.a = 1.;

      float occ = calcAO(pos, nor);
      float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0  );
      float dif = diffuse(nor, lightPos);

      dif *= min(0.3 + softshadow(pos, lightPos, 0.02, 1.5), 1.);
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
      color *= sqrt(1.75 - .5 * dot(rayDirection, vec3(0., 0., -1.)) * 2.5);
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
