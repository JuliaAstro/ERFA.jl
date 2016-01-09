using BinDeps
@BinDeps.setup

version = "1.2.0"
url = "https://github.com/liberfa/erfa/releases/download/v$version/erfa-$version.tar.gz"

# This function returns true if a library satisfies our
# requirements. BinDeps passes to this function the library path and
# the dlopen'ed handle. In this case, we require the library have the
# symbol `eraAb`, which is not present in liberfa 1.0.0.
validate(l, h) = Libdl.dlsym_e(h, "eraAb") != C_NULL

erfa = library_dependency("liberfa"; validate = validate)
provides(Sources, URI(url), erfa)
provides(BuildProcess, Autotools(libtarget="src/liberfa.la"), erfa)

@BinDeps.install Dict(:liberfa => :liberfa)
