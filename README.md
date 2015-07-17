# GSW-Haskell
Haskell interface to Gibbs-SeaWater (GSW) Oceanographic Toolbox in Fortran.
Although [the development version](https://github.com/TEOS-10/GSW-Fortran)
appears far ahead of [the published version](http://www.teos-10.org/software.htm),
the latter is used in this implementation.

## Install

The interfaces (i.e. `foreign import`) are in Oceanogr/GSWtools.hs.
The attached version is built for [GSW Fortran version 3.03](http://www.teos-10.org/software/gsw_fortran_v3_03.zip). If this is the version to be used

    % cd GSW-Haskell
    % mkdir gsw_fortran
    % (cd gsw_fortran && unzip $(DOWNLOAD)/gsw_fortran_v3_03.zip)
    % edit the location of gsw_data_v3.0.dat if necessary (see below)
    % make obj
    % stack build
    % stack test

will build the library. Building and testing rely on [stack](https://github.com/commercialhaskell/stack). See files under Test/ for usage.


If a different version from 3.03 is to be used;

- Download the TEOS-10 Fortran version and place it under gsw_fortran.

    ```
    % cd GSW-Haskell
    % mkdir gsw_fortran
    % cd gsw_fortran
    % unzip $(DOWNLOAD)/gsw_fortran_vX_Y.zip
    ```

- The Fortran toolbox needs to locate the data file *gsw_data_v3_0.dat*.

        ```
        --- gsw_oceanographic_toolbox.f90.dist  2015-07-14 10:38:42.000000000 +0900
        +++ gsw_oceanographic_toolbox.f90       2015-07-14 10:39:33.000000000 +0900
        @@ -3835,7 +3835,7 @@
        if(icalled.eq.0d0) then
        icalled = 1
        -   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
        +   open(10,file='/opt/lib/GSW/gsw_fortran/gsw_data_v3_0.dat',status='old',err=1)
        flag_saar = 1
        read(10,*) (longs_ref(i), i=1,nx)
        read(10,*) (lats_ref(i), i=1,ny)
        ```

- `% make gentool` will yield generate GSWtools.hs. Edit this.

- When ready, `% mv GSWtools.hs Oceangr/`.
