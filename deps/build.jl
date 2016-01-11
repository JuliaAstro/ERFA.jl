using BinDeps
using Compat
@BinDeps.setup

version = "1.2.0"
url = "https://github.com/liberfa/erfa/releases/download/v$version/erfa-$version.tar.gz"

# This function returns true if a library satisfies our
# requirements. BinDeps passes to this function the library path and
# the dlopen'ed handle. We test is certain symbols we need exist.
validate(l, h) = (Libdl.dlsym_e(h, "eraAb") != C_NULL &&
                  Libdl.dlsym_e(h, "eraG2icrs") != C_NULL &&
                  Libdl.dlsym_e(h, "eraIcrs2g") != C_NULL)

erfa = library_dependency("liberfa"; validate = validate)
provides(Sources, URI(url), erfa)
provides(BuildProcess, Autotools(libtarget="src/liberfa.la"), erfa)

@BinDeps.install @compat Dict(:liberfa => :liberfa)
