using BinDeps
@BinDeps.setup

version = "1.2.0"
url = "https://github.com/liberfa/erfa/releases/download/v$version/erfa-$version.tar.gz"

erfa = library_dependency("liberfa")
provides(Sources, URI(url), erfa)
provides(BuildProcess, Autotools(libtarget="src/liberfa.la"), erfa)

@BinDeps.install [:liberfa => :liberfa]
