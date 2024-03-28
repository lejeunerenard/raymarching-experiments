const glslify = require('glslify')

// Base GLSL
const SVG_BASE = glslify('./base-svg-sdf.glsl')

function createSVGShader (svgPaths) {
  return SVG_BASE.replace(/\s*#define SVG_PATHS 1/, '\n' + svgPaths)
}

module.exports = function createSVGSDF (svgPaths) {
  return {
    vertex: glslify('../vert.glsl'),
    fragment: createSVGShader(svgPaths)
  }
}
