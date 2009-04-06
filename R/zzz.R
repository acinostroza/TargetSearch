.First.lib <- function(lib, pkg) {
  library.dynam("TargetSearch", pkg, lib)
}