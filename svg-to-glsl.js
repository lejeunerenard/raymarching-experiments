const { parse } = require('svg-parser')
const parseSVGPath = require('svg-path-parser')

function float (number) {
  let str = number + ''
  if (!str.match(/\./)) {
    str += '.'
  }

  return str
}

function walk (node, metadata = {}) {
  if (node.type === 'element' && node.tagName === 'path') {
    if ('d' in node.properties) {
      const pathCommands = parseSVGPath(node.properties.d)
      console.log('pathCommands', pathCommands)
      const offset = metadata.offset ? metadata.offset : { x: 0, y: 0 }

      let bounds = { x: { min: Infinity, max: -Infinity }, y: { min: Infinity, max: -Infinity } }
      let pathGLSL = 'path(\n'

      // Style
      pathGLSL += '  style(FILL, 0x000000, DIR)\n'

      // Commands
      let hasMidPathClose = false
      let prevCommand = null
      let prevControlPoint = {}
      let currentPosition = {}
      for (let i = 0; i < pathCommands.length; i++) {
        const command = pathCommands[i]

        // Translate coordinate system
        if (!command.relative) {
          // Xs
          if ('x' in command) {
            command.x += offset.x
          }
          if ('x1' in command) {
            command.x1 += offset.x
          }
          if ('x2' in command) {
            command.x2 += offset.x
          }

          // Ys
          if ('y' in command) {
            command.y += offset.y
          }
          if ('y1' in command) {
            command.y1 += offset.y
          }
          if ('y2' in command) {
            command.y2 += offset.y
          }
        }

        switch (command.code) {
          case 'M':
            pathGLSL += `  M(${float(command.x)},${float(command.y)})\n`

            // Movement
            currentPosition.x = command.x
            currentPosition.y = command.y
            break
          case 'm':
            pathGLSL += `  m(${float(command.x)},${float(command.y)})\n`

            // Movement
            currentPosition.x += command.x
            currentPosition.y += command.y
            break
          case 'l':
            pathGLSL += `  l(${float(command.x)},${float(command.y)})\n`

            // Movement
            currentPosition.x += command.x
            currentPosition.y += command.y
            break
          case 'L':
            pathGLSL += `  L(${float(command.x)},${float(command.y)})\n`

            // Movement
            currentPosition.x = command.x
            currentPosition.y = command.y
            break
          case 'c':
            pathGLSL += `  c(${float(command.x1)},${float(command.y1)},${float(command.x2)},${float(command.y2)},${float(command.x)},${float(command.y)})\n`
            prevControlPoint.x = currentPosition.x + command.x2
            prevControlPoint.y = currentPosition.y + command.y2

            // Movement
            currentPosition.x += command.x
            currentPosition.y += command.y
            break
          case 'C':
            pathGLSL += `  C(${float(command.x1)},${float(command.y1)},${float(command.x2)},${float(command.y2)},${float(command.x)},${float(command.y)})\n`
            prevControlPoint.x = command.x2
            prevControlPoint.y = command.y2

            // Movement
            currentPosition.x = command.x
            currentPosition.y = command.y
            break
          case 's':
            const controlPoint1 = {}
            if (prevCommand.code.toLowerCase() === 'c' ||
              prevCommand.code.toLowerCase() === 's') {
              // Compute control point as mirrored from previous control point and as relative
              // delta = prevControl - current
              // mirroredDelta = -delta
              // newControl = current + mirroredDelta
              // newControl = current - prevControl + current
              // newRelativeControl + current = newControl
              // newRelativeControl = newControl - current
              // newRelativeControl = current - prevControl + current - current
              // newRelativeControl = current - prevControl
              controlPoint1.x = currentPosition.x - prevControlPoint.x
              controlPoint1.y = currentPosition.y - prevControlPoint.y
            }

            pathGLSL += `  c(${float(controlPoint1.x)},${float(controlPoint1.y)},${float(command.x2)},${float(command.y2)},${float(command.x)},${float(command.y)})\n`

            prevControlPoint.x = currentPosition.x + command.x2
            prevControlPoint.y = currentPosition.y + command.y2

            // Movement
            currentPosition.x += command.x
            currentPosition.y += command.y
            break
          case 'h':
            pathGLSL += `  h(${float(command.x)})\n`

            // Movement
            currentPosition.x += command.x
            break
          case 'H':
            pathGLSL += `  H(${float(command.x)})\n`

            // Movement
            currentPosition.x = command.x
            break
          case 'v':
            pathGLSL += `  v(${float(command.y)})\n`

            // Movement
            currentPosition.y += command.y
            break
          case 'V':
            pathGLSL += `  V(${float(command.y)})\n`

            // Movement
            currentPosition.y = command.y
            break
          case 'Z':
            if (i !== pathCommands.length - 1) hasMidPathClose = true

            // Disabling the `Z` fixes a rendering issue that only shows when a
            // path without `Z`s mid way through it ends on a `Z`. This maybe
            // is caused by the last command ending on the starting coordinates
            // instead of relying on the `Z` command to complete the path?
            const disable = i !== pathCommands.length - 1 || hasMidPathClose ? '' : '// '
            pathGLSL += `  ${disable}Z\n`
            break
          case 'z':
            pathGLSL += `  svgz\n`
            break
          default:
            console.log('skipped command', command)
            break
        }

        prevCommand = command

        bounds.x.min = Math.min(bounds.x.min, currentPosition.x)
        bounds.y.min = Math.min(bounds.y.min, currentPosition.y)
        bounds.x.max = Math.max(bounds.x.max, currentPosition.x)
        bounds.y.max = Math.max(bounds.y.max, currentPosition.y)
      }

      pathGLSL += ')\n'

      console.log('bounds', bounds)
      return { glsl: pathGLSL, metadata }
    }
  } else if (node.type === 'element' && node.tagName === 'style') {
    return { glsl: '', metadata }
  } else {
    let output = ''

    if (node.type === 'element' && node.tagName === 'svg' && 'viewBox' in node.properties) {
      const coordinates = node.properties.viewBox.split(' ')
      metadata.viewBox = {
        x1: coordinates[0],
        y1: coordinates[1],
        x2: coordinates[2],
        y2: coordinates[3]
      }
      metadata.offset = metadata.offset || { x: 0, y: 0 }
      metadata.offset.x = -(metadata.viewBox.x2 - metadata.viewBox.x1) / 2
      metadata.offset.y = -(metadata.viewBox.y2 - metadata.viewBox.y1) / 2
      console.log('metadata.viewBox', metadata.viewBox)
    }

    for (let child of node.children) {
      const { glsl, metadata: childMetadata } = walk(child, metadata)
      metadata = childMetadata
      output += glsl
    }

    return { glsl: output, metadata }
  }
}

module.exports = function convert (svgText) {
  const parsed = parse(svgText)
  console.log('parsed', parsed)

  return walk(parsed)
}
