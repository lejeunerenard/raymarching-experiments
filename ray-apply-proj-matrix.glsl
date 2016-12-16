#pragma glslify: inverse = require(glsl-inverse)

// getRayDirection source:
// https://bitbucket.org/jimbo00000/opengl-with-luajit/src/33a9abcd04a468d9974d7b35f869d3b0e8e696f4/scene/shadertoy_scene2.lua?at=master&fileviewer=file-view-default#shadertoy_scene2.lua-215
// Construct the usual eye ray frustum oriented down the negative z axis.
// http://antongerdelan.net/opengl/raycasting.html
vec3 getRayDirection(vec2 uv, in mat4 projectionMatrix) {
    vec4 ray_clip = vec4(uv.x, uv.y, -1., 1.);
    vec4 ray_eye = inverse(projectionMatrix) * ray_clip;
    return normalize(vec3(ray_eye.x, ray_eye.y, -1.));
}

#pragma glslify: export(getRayDirection)
