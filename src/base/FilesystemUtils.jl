module FilesystemUtils

# ---------------------

export toext

# ---------------------

# change the path to have the given extension
toext(path::String, ext::String) = splitext(path)[1] * _cleanext(ext)

# append the dot for an extension if need be
_cleanext(ext::String) = ext[1] == '.' ? ext : "." * ext

# ---------------------

end # end module