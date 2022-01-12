precision highp float;

uniform vec2 resolution;
uniform sampler2D buffer;
uniform float minBright;

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);

  vec4 color = texture2D(buffer, uv);
  vec3 colorLinear = pow(color.rgb, vec3(2.2));


  // TODO Finish debugging bloom filter intensity

  // Check whether fragment output is higher than threshold, if so output as brightness color
  // float brightnessStandard = dot(colorLinear.rgb, vec3(0.2126, 0.7152, 0.0722));
  // ITU BT.601
  // source: http://stackoverflow.com/a/596243/630490
  float brightnessPerceived = dot(colorLinear.rgb, vec3(0.299, 0.587, 0.114));

  // float selectiveBrightness = 1.414214 - length(colorLinear.rgb - #FF0000);

  float brightness = brightnessPerceived;

  if (brightness >= minBright) {
    gl_FragColor = vec4(color.rgb, 1.0);
  } else {
    gl_FragColor = vec4(vec3(0.), 1.);
  }
}
