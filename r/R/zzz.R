.onUnload = function(libpath) {
  message("Unloading tumopp")
  library.dynam.unload("tumopp", libpath)
}
