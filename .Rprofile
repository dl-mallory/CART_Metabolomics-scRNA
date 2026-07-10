# Put the pinned library ahead of the system library BEFORE any library() call.
# msigdbr must resolve to 24.1.0 (MSigDB v2024.1); see env/install.R.
local({
    lib <- file.path(getwd(), "env", "rlib")
    if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
})
