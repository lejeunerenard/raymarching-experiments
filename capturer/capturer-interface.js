const assert = require('assert')

class Capturer {
  constructor (opt = {}) {
    assert.strictEqual(typeof opt, 'object', 'Options must be an object')
    assert.ok('filename' in opt, 'Options must define a "filename"')
    assert.ok('framerate' in opt, 'Options must define a "framerate"')
    assert.ok('quality' in opt, 'Options must define a "quality"')
    assert.ok('startTime' in opt, 'Options must define a "startTime"')
    assert.ok('timeLength' in opt, 'Options must define a "timeLength"')

    Object.assign(this, opt)
  }

  start () {
    assert.fail('start() was not defined for renderer')
  }

  async capture (canvas) {
    assert.fail('capture() was not defined for renderer')
  }

  async save () {
    assert.fail('save() was not defined for renderer')
  }
}

module.exports = Capturer
