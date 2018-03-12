# GSW-Haskell
Haskell interface to [Gibbs-SeaWater (GSW) Oceanographic Toolbox](https://www.teos-10.org/).
Previous versions were interface to FORTRAN toolbox, but from version 2.0.0.0,
this Haskell version uses the C toolbox. The latest version of the C toolbox
at the time of writing is 3.05.0-4. If a different version is to be used, it might
be necessary to tweak Oceanogr/GSWtools.chs.

## Install
Download the C version of the GSW Toolbox. Usually the latest from [Github repo](https://github.com/TEOS-10/GSW-C.git) is newer than [official repo](https://www.teos-10.org/software.htm). Install in your favourite directory $(LIB).

    % cd $(LIB)
    % git clone https://github.com/TEOS-10/GSW-C.git
    % mv GSW-C gsw_c_v3.05_0-4

If you do not have a Python environment, apply the following patch to avoid loosing `gsw_chec_data.c`.

   --- GNUmakefile.dist    2017-04-27 12:28:18.273726163 +0900
   +++ GNUmakefile 2017-04-27 12:28:36.329540133 +0900
   @@ -83,9 +83,9 @@
           rm -f $(Library)
           ln -s $(Library).$(LibVersion) $(Library)
    
   -gsw_check_data.c:      $(GSW_3_DATA)
   -                       rm -f $@; \
   -                       ./make_check_data.py
   +#gsw_check_data.c:     $(GSW_3_DATA)
   +#                      rm -f $@; \
   +#                      ./make_check_data.py
    
    gsw_saar_data.c:       $(GSW_3_DATA)
                           rm -f $@; \

It's a good idea to test the C version.

    % cd gsw_c_v3.05_0-4
    % ./gsw_check

The interfaces (i.e. `foreign import`s) are automatically generated from Oceanogr/GSWtools.chs by [c2hs](https://wiki.haskell.org/C2hs). The attached version is for [GSW C version 3.05](http://www.teos-10.org/software/gsw_C_v3_05.zip). It might be necessary to modify Oceanogr/GSWtools.chs if different version is to be used.

Before compiling, modify `include-dirs` and `extra-lib-dirs` in GSW.cabal. Note that the paths must be absolute ([issue](https://github.com/haskell/cabal/issues/2641) with cabal). If the C library in $(LIB) cannot be found at runtime, it is possible to add $(LIB) to the runtime library search path by

    ghc-options: -optl-Wl,-rpath,$(LIB)/gsw_c_v3.05/

entry in the GSW.cabal file.
Then,

    % stack install
    
will install the library.

### Changes
0.2.0.3 Update for gsw_C_V3_05.0-4 and stack LTS-9.10

0.2.0.2 Bug fixes

0.2.0.1 With gsw_c_v3.05_1.zip

0.2.0.0 With gsw_c_v3.05.zip
