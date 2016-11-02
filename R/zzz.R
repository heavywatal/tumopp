.onUnload = function (libpath) {
    message('Unloading tumorr')
    library.dynam.unload('tumorr', libpath)
}
