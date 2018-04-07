import Octavian from 'octavian'

function createSwell (audioCtx, note, start, length) {
  start += audioCtx.currentTime

  const baseNote = new Octavian.Note(note)
  const offNote = baseNote.perfectFourth()

  const baseGain = audioCtx.createGain()
  baseGain.gain.setValueAtTime(1, audioCtx.currentTime)

  const osc1 = audioCtx.createOscillator()
  osc1.type = 'sawtooth'
  osc1.frequency.setValueAtTime(baseNote.frequency, audioCtx.currentTime)

  const filter1 = audioCtx.createBiquadFilter()
  filter1.type = 'lowpass'
  filter1.frequency.setValueAtTime(0, start)
  filter1.frequency.linearRampToValueAtTime(1000, start + length / 2)
  filter1.frequency.linearRampToValueAtTime(0, start + length)
  filter1.detune.setValueAtTime(0, start)
  filter1.detune.linearRampToValueAtTime(400, start + length / 2)
  filter1.detune.linearRampToValueAtTime(0, start + length)

  osc1.connect(filter1)
  filter1.connect(baseGain)

  const osc2 = audioCtx.createOscillator()
  osc2.type = 'sawtooth'
  osc2.frequency.setValueAtTime(offNote.frequency - 1.5, audioCtx.currentTime)

  const filter2 = audioCtx.createBiquadFilter()
  filter2.type = 'lowpass'
  filter2.frequency.setValueAtTime(0, start)
  filter2.frequency.linearRampToValueAtTime(880, start + length / 2)
  filter2.frequency.linearRampToValueAtTime(0, start + length)

  osc2.connect(filter2)
  filter2.connect(baseGain)

  osc1.start(start)
  osc2.start(start)

  return baseGain
}

module.exports = createSwell
