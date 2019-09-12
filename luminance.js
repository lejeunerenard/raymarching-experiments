const color = [0, 0, 0]

const redLuminanceCoefficient = 0.2126
const greenLuminanceCoefficient = 0.7152
const blueLuminanceCoefficient = 0.0722

function getLuminance (color) {
  return redLuminanceCoefficient * color[0] +
    greenLuminanceCoefficient * color[1] +
    blueLuminanceCoefficient * color[2]
}

function getColorWFixedLuminance (red, green, blue, targetLuminance) {
  let luminanceCoefficient = 0

  if (!(!!red || !!green || !!blue)) {
    throw Error('One color must be omitted to calculate the new color')
  }

  if (!red) {
    luminanceCoefficient = redLuminanceCoefficient
    console.log('select red')
  }

  if (!green) {
    luminanceCoefficient = greenLuminanceCoefficient
    console.log('select green')
  }

  if (!blue) {
    luminanceCoefficient = blueLuminanceCoefficient
    console.log('select blue')
  }

  let tempLuminence = getLuminance([red || 0, green || 0, blue || 0])

  if (tempLuminence > targetLuminance) throw Error('Target luminance is greater than provided color')

  return (targetLuminance - tempLuminence) / luminanceCoefficient
}

// console.log(getColorWFixedLuminance(0.5, null, 0.8, 42 / 100)) // 42
// console.log(getColorWFixedLuminance(null, 0.8, 1.0, 80 / 100)) // 80

console.log(getColorWFixedLuminance(1, null, 0.1, 42 / 100)) // 42
console.log(getColorWFixedLuminance(null, 0.9, 0.4, 80 / 100)) // 80
console.log(getLuminance([1, 0x6A / 255, 0x95 / 255]))
