const glslify = require('glslify')

module.exports = require('shader-reload')({
  vertex: glslify('../vert.glsl'),
  fragment: glslify('../bright.glsl')
})
