using BinDeps
@BinDeps.setup

version = "1.3.0"
url = "https://github.com/liberfa/erfa/releases/download/v$version/erfa-$version.tar.gz"

# This function returns true if a library satisfies our
# requirements. BinDeps passes to this function the library path and
# the dlopen'ed handle. We test is certain symbols we need exist.
validate(l, h) = (Libdl.dlsym_e(h, "eraAb") != C_NULL &&
                  Libdl.dlsym_e(h, "eraG2icrs") != C_NULL &&
                  Libdl.dlsym_e(h, "eraIcrs2g") != C_NULL &&
                  Libdl.dlsym_e(h, "eraEceq06") != C_NULL &&
                  Libdl.dlsym_e(h, "eraLtpb") != C_NULL)

erfa = library_dependency("liberfa"; validate = validate)

provides(Sources, URI(url), erfa)
provides(BuildProcess, Autotools(libtarget="src/liberfa.la"), erfa)

# Windows binaries cross-compiled in erfa with
#
# mkdir -p /usr/lib
# i686-w64-mingw32-gcc -Isrc -o usr/lib/liberfa.dll -O3 -shared \
#     -static-libgcc src/*.c
# 7za a erfa-win32.7z usr
# x86_64-w64-mingw32-gcc -Isrc -o usr/lib/liberfa.dll -O3 -shared \
#     -static-libgcc src/*.c
# 7za a erfa-win64.7z usr

provides(Binaries, URI("https://dl.bintray.com/kbarbary/generic/erfa-$(version)-win$(Sys.WORD_SIZE).7z"), erfa, os = :Windows)

@BinDeps.install Dict(:liberfa => :liberfa)
