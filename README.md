# GSW-Haskell
Haskell interface to [Gibbs-SeaWater (GSW) Oceanographic Toolbox](https://www.teos-10.org/).
Previous versions were interface to FORTRAN toolbox, but from version 2.0.0.0,
this Haskell version uses the C toolbox. The latest version of the C toolbox
at the time of writing is 3.05. If a different version is to be used, it might
be necessary to tweak Oceanogr/GSWtools.chs.

## Install
Download the C version of the GSW Toolbox from [here](https://www.teos-10.org/software.htm). Install in your favourite directory $(LIB).

    % cd $(LIB)
    % unzip $(SRC)/gsw_c_v3.05.zip

It's a good idea to test the C version.

    % cd gsw_c_v3.05
    % make
    % ./gsw_check

A patch is recommended to deal with a [problem](https://github.com/TEOS-10/GSW-C/pull/2) related to land value of saar.

    % patch < $(SOMEWHERE)/patch.gsw

Build the shared library

    % make library

The interfaces (i.e. `foreign import`s) are automatically generated from Oceanogr/GSWtools.chs by [c2hs](https://wiki.haskell.org/C2hs). The attached version is for [GSW C version 3.05](http://www.teos-10.org/software/gsw_C_v3_05.zip). It might be necessary to modify Oceanogr/GSWtools.chs if different version is to be used.

Before compiling, modify `include-dirs` and `extra-lib-dirs` in GSW.cabal. Note that the paths must be absolute ([issue](https://github.com/haskell/cabal/issues/2641) with cabal). If the C library in $(LIB) cannot be found at runtime, it is possible to add $(LIB) to the runtime library search path by

    ghc-options: -optl-Wl,-rpath,$(LIB)/gsw_c_v3.05/

entry in the GSW.cabal file.
Then,

    % stack install
    
will install the library.

### Changes
0.2.0.2 Bug fixes

0.2.0.1 With gsw_c_v3.05_1.zip

0.2.0.0 With gsw_c_v3.05.zip
