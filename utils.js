export function isIOS () {
  return !!navigator.platform && /iPad|iPhone|iPod/.test(navigator.platform)
}
export function isAndroid () {
  let ua = navigator.userAgent.toLowerCase()
  return ua.indexOf("android") > -1
}
