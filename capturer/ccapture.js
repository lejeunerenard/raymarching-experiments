const Capturer = require('./capturer-interface')

class CCaptureCapturer extends Capturer {
  constructor (opt = {}) {
    super(opt)

    this.capturer = new CCapture({
      format: opt.format || 'jpg',
      framerate: this.framerate,
      name: this.filename,
      autoSaveTime: opt.autoSaveTime || 5,
      quality: this.quality,
      startTime: this.startTime,
      timeLimit: this.timeLength,
      verbose: true
    })
  }

  start () {
    this.capturer.start()
  }

  async capture (canvas) {
    this.capturer.capture(canvas)
  }

  async save () {
    console.log('ccapture.js save')
  }
}

module.exports = CCaptureCapturer
