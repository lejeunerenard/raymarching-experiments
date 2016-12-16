precision highp float;

attribute vec2 position;
uniform vec2 resolution;
varying vec2 fragCoord;

void main() {
    // Adjust coordinates to (-1, -1) -> (1, 1)
    fragCoord = position.xy * resolution / max(resolution.y, resolution.x);

    gl_Position = vec4(position.xy, 0.0, 1.0);
}
