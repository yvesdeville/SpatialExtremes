".First.lib" <-
function(lib, pkg)
{
  library.dynam("SpatialExtremes", package = pkg, lib.loc = lib)
  return(invisible(0))
}

