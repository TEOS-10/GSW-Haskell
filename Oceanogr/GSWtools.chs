{-# LANGUAGE ForeignFunctionInterface #-}
--
-- | Raw binding for Gibbs-Seawater Oceanographic Tools C version (gsw_c_v3.05)
--
-- The `pure` keyword in `#fun` hook seems to be implemented with `unsafePerformIO`
-- (maybe in response to [Issue 130](https://github.com/haskell/c2hs/issues/130)?)
-- I really do not understand how `unsafePerformIO` works (see e.g. [禁断の機能...](http://itpro.nikkeibp.co.jp/article/COLUMN/20090512/329783/)).
-- So I did not make these functions `pure` (instead `unsafe`).
-- If performance really matters, consider using `pure` for better optimisation.
--

module Oceanogr.GSWtools where

import Foreign
import Foreign.C.Types
-- import Foreign.Ptr
-- import Foreign.Storable
-- import Foreign.Marshal.Alloc
-- import Foreign.Marshal.Array

#include "gswteos-10.h"
--
-- Error values
--
gswInvalidValue, gswErrorLimit :: Double
gswInvalidValue = 9e15
gswErrorLimit   = 1e10

--
-- Marshalling
--

-- The following function is copied (and modified) from
-- C2HS.hs which 'is DEPRECATED "The C2HS module should no longer be used."' (version 0.28.1)

-- peekFloatConv :: (Storable a, RealFloat a, RealFloat b) 
--              => Ptr a -> IO b
peekFloatConv :: Ptr CDouble -> IO Double
peekFloatConv = fmap realToFrac . peek
--
--
withArrayConv :: [Double] -> (Ptr CDouble -> IO a) -> IO a
withArrayConv = withArray . map realToFrac

withArrayConvI :: [Int] -> (Ptr CInt -> IO a) -> IO a
withArrayConvI = withArray . map fromIntegral

allocaArray4 :: (Ptr CDouble -> IO a) -> IO a
allocaArray4 = allocaArray 4

peekArray4 :: Ptr CDouble -> IO [Double]
peekArray4 d = map realToFrac `fmap` peekArray 4 d

-- void
-- gsw_add_barrier(double *input_data, double lon, double lat,
--		double long_grid, double lat_grid, double dlong_grid,
--		double dlat_grid, double *output_data)
{#fun unsafe gsw_add_barrier
        { withArrayConv* `[Double]',
         `Double', `Double',
         `Double', `Double',
         `Double', `Double',
          allocaArray4- `[Double]' peekArray4*} -> `()' #}

-- void
-- gsw_add_mean(double *data_in, double *data_out)
{#fun unsafe gsw_add_mean
        { withArrayConv* `[Double]',
          allocaArray4- `[Double]' peekArray4*} -> `()' #}

-- double
-- gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p)
{#fun unsafe gsw_adiabatic_lapse_rate_from_ct
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_adiabatic_lapse_rate_ice(double t, double p)
{#fun unsafe gsw_adiabatic_lapse_rate_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_alpha(double sa, double ct, double p)
{#fun unsafe gsw_alpha
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_alpha_on_beta(double sa, double ct, double p)
{#fun unsafe gsw_alpha_on_beta
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_alpha_wrt_t_exact(double sa, double t, double p)
{#fun unsafe gsw_alpha_wrt_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_alpha_wrt_t_ice(double t, double p)
{#fun unsafe gsw_alpha_wrt_t_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_beta(double sa, double ct, double p)
{#fun unsafe gsw_beta
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_beta_const_t_exact(double sa, double t, double p)
{#fun unsafe gsw_beta_const_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_c_from_sp(double sp, double t, double p)
{#fun unsafe gsw_c_from_sp
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_cabbeling(double sa, double ct, double p)
{#fun unsafe gsw_cabbeling
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_chem_potential_water_ice(double t, double p)
{#fun unsafe gsw_chem_potential_water_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_chem_potential_water_t_exact(double sa, double t, double p)
{#fun unsafe gsw_chem_potential_water_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_cp_ice(double t, double p)
{#fun unsafe gsw_cp_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_cp_t_exact(double sa, double t, double p)
{#fun unsafe gsw_cp_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_ct_first_derivatives(double sa, double pt, double *ct_sa, double *ct_pt)
{#fun unsafe gsw_ct_first_derivatives
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p, 
-- 	double *ct_sa_wrt_t, double *ct_t_wrt_t, double *ct_p_wrt_t)
{#fun unsafe gsw_ct_first_derivatives_wrt_t_exact
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_ct_freezing(double sa, double p, double saturation_fraction)
{#fun unsafe gsw_ct_freezing
        {`Double', `Double', `Double'} -> `Double' #}

-- removed (https://github.com/TEOS-10/GSW-C/pull/12)
-- double
-- gsw_ct_freezing_exact(double sa, double p, double saturation_fraction)
-- {#fun unsafe gsw_ct_freezing_exact
--        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_ct_freezing_first_derivatives(double sa, double p,
--       double saturation_fraction, double *ctfreezing_sa, double *ctfreezing_p)
{#fun unsafe gsw_ct_freezing_first_derivatives
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_ct_freezing_first_derivatives_poly(double sa, double p,
--		double saturation_fraction, double *ctfreezing_sa,
--		double *ctfreezing_p)
{#fun unsafe gsw_ct_freezing_first_derivatives_poly
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_ct_freezing_poly(double sa, double p, double saturation_fraction)
{#fun unsafe gsw_ct_freezing_poly
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_ct_from_enthalpy(double sa, double h, double p)
{#fun unsafe gsw_ct_from_enthalpy
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_ct_from_enthalpy_exact(double sa, double h, double p)
{#fun unsafe gsw_ct_from_enthalpy_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_ct_from_entropy(double sa, double entropy)
{#fun unsafe gsw_ct_from_entropy
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_ct_from_pt(double sa, double pt)
{#fun unsafe gsw_ct_from_pt
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_ct_from_rho(double rho, double sa, double p, double *ct,
--		double *ct_multiple)
{#fun unsafe gsw_ct_from_rho
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_ct_from_t(double sa, double t, double p)
{#fun unsafe gsw_ct_from_t
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_ct_maxdensity(double sa, double p)
{#fun unsafe gsw_ct_maxdensity
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_ct_second_derivatives(double sa, double pt, double *ct_sa_sa,
-- 	double *ct_sa_pt, double *ct_pt_pt)
{#fun unsafe gsw_ct_second_derivatives
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_deltasa_from_sp(double sp, double p, double lon, double lat)
{#fun unsafe gsw_deltasa_from_sp
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_dilution_coefficient_t_exact(double sa, double t, double p)
{#fun unsafe gsw_dilution_coefficient_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_dynamic_enthalpy(double sa, double ct, double p)
{#fun unsafe gsw_dynamic_enthalpy
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_enthalpy(double sa, double ct, double p)
{#fun unsafe gsw_enthalpy
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_enthalpy_ct_exact(double sa, double ct, double p)
{#fun unsafe gsw_enthalpy_ct_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_enthalpy_diff(double sa, double ct, double p_shallow, double p_deep)
{#fun unsafe gsw_enthalpy_diff
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_enthalpy_first_derivatives(double sa, double ct, double p, double *h_sa,
--				double *h_ct)
{#fun unsafe gsw_enthalpy_first_derivatives
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct, double p,
--		double *h_sa, double *h_ct)
{#fun unsafe gsw_enthalpy_first_derivatives_ct_exact
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_enthalpy_ice(double t, double p)
{#fun unsafe gsw_enthalpy_ice
        {`Double', `Double'} -> `()' #}

-- void
-- gsw_enthalpy_second_derivatives(double sa, double ct, double p,
--		double *h_sa_sa, double *h_sa_ct, double *h_ct_ct)
{#fun unsafe gsw_enthalpy_second_derivatives
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_enthalpy_second_derivatives_ct_exact(double sa, double ct, double p,
--		double *h_sa_sa, double *h_sa_ct, double *h_ct_ct)
{#fun unsafe gsw_enthalpy_second_derivatives_ct_exact
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_enthalpy_sso_0(double p)
{#fun unsafe gsw_enthalpy_sso_0
        {`Double'} -> `Double' #}

-- double
-- gsw_enthalpy_t_exact(double sa, double t, double p)
{#fun unsafe gsw_enthalpy_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_entropy_first_derivatives(double sa, double ct, double *eta_sa,
-- 	double *eta_ct)
{#fun unsafe gsw_entropy_first_derivatives
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_entropy_from_ct(double sa, double ct)
{#fun unsafe gsw_entropy_from_ct
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_entropy_from_pt(double sa, double pt)
{#fun unsafe gsw_entropy_from_pt
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_entropy_from_t(double sa, double t, double p)
{#fun unsafe gsw_entropy_from_t
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_entropy_ice(double t, double p)
{#fun unsafe gsw_entropy_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_entropy_part(double sa, double t, double p)
{#fun unsafe gsw_entropy_part
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_entropy_part_zerop(double sa, double pt0)
{#fun unsafe gsw_entropy_part_zerop
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_entropy_second_derivatives(double sa, double ct,
-- 	double *eta_sa_sa, double *eta_sa_ct, double *eta_ct_ct)
{#fun unsafe gsw_entropy_second_derivatives
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_fdelta(double p, double lon, double lat)
{#fun unsafe gsw_fdelta
        {`Double', `Double', `Double'} -> `Double' #}


-- void
-- gsw_frazil_properties(double sa_bulk, double h_bulk, double p,
--	double *sa_final, double *ct_final, double *w_ih_final)
{#fun unsafe gsw_frazil_properties
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_frazil_properties_potential(double sa_bulk, double h_pot_bulk, double p,
-- 	double *sa_final, double *ct_final, double *w_ih_final)
{#fun unsafe gsw_frazil_properties_potential
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_frazil_properties_potential_poly(double sa_bulk, double h_pot_bulk,
-- 	double p, double *sa_final, double *ct_final, double *w_ih_final)
{#fun unsafe gsw_frazil_properties_potential_poly
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_frazil_ratios_adiabatic (double sa, double p, double w_ih,
-- 	double *dsa_dct_frazil, double *dsa_dp_frazil, double *dct_dp_frazil)
{#fun unsafe gsw_frazil_ratios_adiabatic
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih,
-- 	double *dsa_dct_frazil, double *dsa_dp_frazil, double *dct_dp_frazil)
{#fun unsafe gsw_frazil_ratios_adiabatic_poly
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double	*	/* Returns NULL on error, dyn_height if okay */
-- gsw_geo_strf_dyn_height(double *sa, double *ct, double *p, double p_ref,
--	int n_levels, double *dyn_height)
gsw_geo_strf_dyn_height :: [Double] -> [Double] -> [Double] -> Double -> Int
                              -> IO (Maybe [Double])
gsw_geo_strf_dyn_height sa ct p pRef nLevels =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
            do dynHeight' <- mallocArray nLevels :: IO (Ptr CDouble)
               ret <- {#call unsafe gsw_geo_strf_dyn_height as ^ #}
                    sa' ct' p' (realToFrac pRef) (fromIntegral nLevels) dynHeight'
               if ret == nullPtr
                 then free dynHeight' >> return Nothing
                 else do dynHeight <- map realToFrac `fmap` peekArray nLevels dynHeight'
                         free dynHeight'
                         return $ Just dynHeight

-- double *
-- gsw_geo_strf_dyn_height_pc(double *sa, double *ct, double *delta_p, int n_levels,
--	double *geo_strf_dyn_height_pc, double *p_mid)
gsw_geo_strf_dyn_height_pc :: [Double] -> [Double] -> [Double] -> Int
                              -> IO (Maybe ([Double], [Double]))
gsw_geo_strf_dyn_height_pc sa ct deltaP nLevels =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv deltaP $ \deltaP' ->
            do geoStrfDynHeightPc' <- mallocArray nLevels :: IO (Ptr CDouble)
               pMid'               <- mallocArray nLevels :: IO (Ptr CDouble)
               ret <- {#call unsafe gsw_geo_strf_dyn_height_pc as ^ #}
                    sa' ct' deltaP' (fromIntegral nLevels) geoStrfDynHeightPc' pMid'
               if ret == nullPtr
                 then do free geoStrfDynHeightPc'
                         free pMid'
                         return Nothing
                 else do geoStrDynHeightPc <-
                                 map realToFrac `fmap` peekArray nLevels geoStrfDynHeightPc'
                         pMid <- map realToFrac `fmap` peekArray nLevels pMid'
                         free geoStrfDynHeightPc'
                         free pMid'
                         return $ Just (geoStrDynHeightPc, pMid)

-- double
-- gsw_gibbs(int ns, int nt, int np, double sa, double t, double p)
{#fun unsafe gsw_gibbs
        {`Int', `Int', `Int',
         `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_grav(double lat, double p)
{#fun unsafe gsw_grav
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_helmholtz_energy_ice(double t, double p)
{#fun unsafe gsw_helmholtz_energy_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_hill_ratio_at_sp2(double t)
{#fun unsafe gsw_hill_ratio_at_sp2
        {`Double'} -> `Double' #}

-- void
-- gsw_ice_fraction_to_freeze_seawater(double sa, double ct, double p, double t_ih,
-- 	double *sa_freeze, double *ct_freeze, double *w_ih)
{#fun unsafe gsw_ice_fraction_to_freeze_seawater
        {`Double', `Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_internal_energy(double sa, double ct, double p)
{#fun unsafe gsw_internal_energy
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_internal_energy_ice(double t, double p)
{#fun unsafe gsw_internal_energy_ice
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p, double p_ref,
--	int nz, double *ipv_vs_fnsquared_ratio, double *p_mid)
gsw_ipv_vs_fnsquared_ratio :: [Double] -> [Double] -> [Double] -> Double -> Int
                                -> IO ([Double], [Double])
gsw_ipv_vs_fnsquared_ratio sa ct p pRef nz =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
          -- IPV_vs_fNsquared_ratio
          --          : The ratio of the vertical gradient of potential density
          --            referenced to p_ref, to the vertical gradient of locally-
          --           referenced potential density.  It is ouput on the same
          --           vertical (M-1)xN grid as p_mid. 
          --           IPV_vs_fNsquared_ratio is dimensionless.          [ unitless ]
         --  p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
          do ipvVsFnsquaredRatio' <- mallocArray (nz-1) :: IO (Ptr CDouble)
             pMid'                <- mallocArray (nz-1) :: IO (Ptr CDouble)
             {#call unsafe gsw_ipv_vs_fnsquared_ratio as ^ #}
                sa' ct' p' (realToFrac pRef) (fromIntegral nz) ipvVsFnsquaredRatio' pMid'
             ipvVsFnsquaredRatio  <- map realToFrac `fmap` peekArray (nz-1) ipvVsFnsquaredRatio'
             pMid                 <- map realToFrac `fmap` peekArray (nz-1) pMid'
             free ipvVsFnsquaredRatio'
             free pMid'
             return (ipvVsFnsquaredRatio, pMid)

-- double
-- gsw_kappa(double sa, double ct, double p)
{#fun unsafe gsw_kappa
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_kappa_const_t_ice(double t, double p)
{#fun unsafe gsw_kappa_const_t_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_kappa_ice(double t, double p)
{#fun unsafe gsw_kappa_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_kappa_t_exact(double sa, double t, double p)
{#fun unsafe gsw_kappa_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_latentheat_evap_ct(double sa, double ct)
{#fun unsafe gsw_latentheat_evap_ct
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_latentheat_evap_t(double sa, double t)
{#fun unsafe gsw_latentheat_evap_t
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_latentheat_melting(double sa, double p)
{#fun unsafe gsw_latentheat_melting
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_linear_interp_sa_ct(double *sa, double *ct, double *p, int np,
--	double *p_i, int npi, double *sa_i, double *ct_i)
gsw_linear_interp_sa_ct :: [Double] -> [Double] -> [Double] -> Int -> [Double] -> Int
                        -> IO ([Double], [Double])
gsw_linear_interp_sa_ct sa ct p np pI npI =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
          withArrayConv pI $ \pI' ->
            do saI' <- mallocArray npI :: IO (Ptr CDouble)
               ctI' <- mallocArray npI :: IO (Ptr CDouble)
               {#call unsafe gsw_linear_interp_sa_ct as ^ #}
                                sa' ct' p' (fromIntegral np) pI' (fromIntegral npI) saI' ctI'
               saI  <- map realToFrac `fmap` peekArray npI saI'
               ctI  <- map realToFrac `fmap` peekArray npI ctI'
               free saI'
               free ctI'
               return (saI, ctI)

-- double
-- gsw_melting_ice_equilibrium_sa_ct_ratio(double sa, double p)
{#fun unsafe gsw_melting_ice_equilibrium_sa_ct_ratio
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa, double p)
{#fun unsafe gsw_melting_ice_equilibrium_sa_ct_ratio_poly
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_melting_ice_into_seawater(double sa, double ct, double p, double w_ih,
-- 	double t_ih, double *sa_final, double *ct_final, double *w_ih_final)
{#fun unsafe gsw_melting_ice_into_seawater
        {`Double', `Double', `Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_melting_ice_sa_ct_ratio(double sa, double ct, double p, double t_ih)
{#fun unsafe gsw_melting_ice_sa_ct_ratio
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_melting_ice_sa_ct_ratio_poly(double sa, double ct, double p, double t_ih)
{#fun unsafe gsw_melting_ice_sa_ct_ratio_poly
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa, double p)
{#fun unsafe gsw_melting_seaice_equilibrium_sa_ct_ratio
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa, double p)
{#fun gsw_melting_seaice_equilibrium_sa_ct_ratio_poly
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_melting_seaice_into_seawater(double sa, double ct, double p,
-- 	double w_seaice, double sa_seaice, double t_seaice,
--	double *sa_final, double *ct_final)
{#fun unsafe gsw_melting_seaice_into_seawater
        {`Double', `Double', `Double',
         `Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_melting_seaice_sa_ct_ratio(double sa, double ct, double p,
--	double sa_seaice, double t_seaice)
{#fun unsafe gsw_melting_seaice_sa_ct_ratio
        {`Double', `Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_melting_seaice_sa_ct_ratio_poly(double sa, double ct, double p,
-- 	double sa_seaice, double t_seaice)
{#fun unsafe gsw_melting_seaice_sa_ct_ratio_poly
        {`Double', `Double', `Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_nsquared(double *sa, double *ct, double *p, double *lat, int nz,
--	double *n2, double *p_mid)
gsw_nsquared :: [Double] -> [Double] -> [Double] -> [Double] -> Int
                -> IO ([Double], [Double])
gsw_nsquared sa ct p lat nz =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
          withArrayConv lat $ \lat' ->
            -- ! n2     : Brunt-Vaisala Frequency squared  (length nz-1)        [s^-2]
            -- ! p_mid  : Mid pressure between p grid      (length nz-1)        [dbar]
            do n2'   <- mallocArray (nz-1) :: IO (Ptr CDouble)
               pMid' <- mallocArray (nz-1) :: IO (Ptr CDouble)
               {#call unsafe gsw_nsquared as ^ #}
                    sa' ct' p' lat' (fromIntegral nz) n2' pMid'
               n2    <- map realToFrac `fmap` peekArray (nz-1) n2'
               pMid  <- map realToFrac `fmap` peekArray (nz-1) pMid'
               free n2'
               free pMid'
               return (n2, pMid)

-- double
-- gsw_p_from_z(double z, double lat)
{#fun unsafe gsw_p_from_z
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_pot_enthalpy_from_pt_ice(double pt0_ice)
{#fun unsafe gsw_pot_enthalpy_from_pt_ice
        {`Double'} -> `Double' #}

-- double
-- gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice)
{#fun unsafe gsw_pot_enthalpy_from_pt_ice_poly
        {`Double'} -> `Double' #}

-- double
-- gsw_pot_enthalpy_ice_freezing(double sa, double p)
{#fun unsafe gsw_pot_enthalpy_ice_freezing
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa, double p,
--     double *pot_enthalpy_ice_freezing_sa, double *pot_enthalpy_ice_freezing_p)
{#fun unsafe gsw_pot_enthalpy_ice_freezing_first_derivatives
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa, double p,
--     double *pot_enthalpy_ice_freezing_sa, double *pot_enthalpy_ice_freezing_p)
{#fun unsafe gsw_pot_enthalpy_ice_freezing_first_derivatives_poly
        {`Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_pot_enthalpy_ice_freezing_poly(double sa, double p)
{#fun unsafe gsw_pot_enthalpy_ice_freezing_poly
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref)
{#fun unsafe gsw_pot_rho_t_exact
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_pressure_coefficient_ice(double t, double p)
{#fun unsafe gsw_pressure_coefficient_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_pressure_freezing_ct(double sa, double ct, double saturation_fraction)
{#fun unsafe gsw_pressure_freezing_ct
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_pt0_cold_ice_poly(double pot_enthalpy_ice)
{#fun unsafe gsw_pt0_cold_ice_poly
        {`Double'} -> `Double' #}

-- double
-- gsw_pt0_from_t(double sa, double t, double p)
{#fun unsafe gsw_pt0_from_t
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_pt0_from_t_ice(double t, double p)
{#fun unsafe gsw_pt0_from_t_ice
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_pt_first_derivatives (double sa, double ct, double *pt_sa, double *pt_ct)
{#fun unsafe gsw_pt_first_derivatives
        {`Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_pt_from_ct(double sa, double ct)
{#fun unsafe gsw_pt_from_ct
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_pt_from_entropy(double sa, double entropy)
{#fun unsafe gsw_pt_from_entropy
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice)
{#fun unsafe gsw_pt_from_pot_enthalpy_ice
        {`Double'} -> `Double' #}

-- double
-- gsw_pt_from_pot_enthalpy_ice_poly
{#fun unsafe gsw_pt_from_pot_enthalpy_ice_poly
        {`Double'} -> `Double' #}

-- double
-- gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice)
{#fun unsafe gsw_pt_from_pot_enthalpy_ice_poly_dh
        {`Double'} -> `Double' #}

-- double
-- gsw_pt_from_t(double sa, double t, double p, double p_ref)
{#fun unsafe gsw_pt_from_t
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_pt_from_t_ice(double t, double p, double p_ref)
{#fun unsafe gsw_pt_from_t_ice
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_pt_second_derivatives (double sa, double ct, double *pt_sa_sa,
-- 	double *pt_sa_ct, double *pt_ct_ct)
{#fun unsafe gsw_pt_second_derivatives
        {`Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_rho(double sa, double ct, double p)
{#fun unsafe gsw_rho
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_rho_alpha_beta (double sa, double ct, double p, double *rho, double *alpha,
--			double *beta)
{#fun unsafe gsw_rho_alpha_beta
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_rho_first_derivatives(double sa, double ct, double p,
--	double *drho_dsa, double *drho_dct, double *drho_dp)
{#fun unsafe gsw_rho_first_derivatives
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_rho_first_derivatives_wrt_enthalpy (double sa, double ct, double p,
--	double *rho_sa, double *rho_h)
{#fun unsafe gsw_rho_first_derivatives_wrt_enthalpy
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_rho_ice(double t, double p)
{#fun unsafe gsw_rho_ice
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_rho_second_derivatives(double sa, double ct, double p, double *rho_sa_sa,
--	double *rho_sa_ct, double *rho_ct_ct, double *rho_sa_p,
--	double *rho_ct_p)
{#fun unsafe gsw_rho_second_derivatives
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_rho_second_derivatives_wrt_enthalpy(double sa, double ct, double p,
--	double *rho_sa_sa, double *rho_sa_h, double *rho_h_h)
{#fun unsafe gsw_rho_second_derivatives_wrt_enthalpy
       {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_rho_t_exact(double sa, double t, double p)
{#fun unsafe gsw_rho_t_exact
       {`Double', `Double', `Double'} -> `()' #}

-- void
-- gsw_rr68_interp_sa_ct(double *sa, double *ct, double *p, int mp, double *p_i,
--	int mp_i, double *sa_i, double *ct_i)
gsw_rr68_interp_sa_ct :: [Double] -> [Double] -> [Double] -> Int -> [Double] -> Int
                        -> IO ([Double], [Double])
gsw_rr68_interp_sa_ct sa ct p mp pI mpI =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
          withArrayConv pI $ \pI' ->
            do saI' <- mallocArray mpI :: IO (Ptr CDouble)
               ctI' <- mallocArray mpI :: IO (Ptr CDouble)
               {#call unsafe gsw_rr68_interp_sa_ct as ^ #}
                                sa' ct' p' (fromIntegral mp) pI' (fromIntegral mpI) saI' ctI'
               saI  <- map realToFrac `fmap` peekArray mpI saI'
               ctI  <- map realToFrac `fmap` peekArray mpI ctI'
               free saI'
               free ctI'
               return (saI, ctI)
---
--- The following function has strange argument passing (double *ct, double *t --
--- although they are INTENT IN).
--- Haskell interface will not be provided.
-- !==========================================================================
-- elemental function gsw_sa_freezing_estimate (p, saturation_fraction, ct, t)
-- !==========================================================================
-- !
-- ! Form an estimate of SA from a polynomial in CT and p 
-- !
-- !--------------------------------------------------------------------------
-- */
--double
--gsw_sa_freezing_estimate(double p, double saturation_fraction, double *ct,
--	double *t)
--{
--	GSW_TEOS10_CONSTANTS;
--	double	ctx, ctsat, sa,
--		/*note that aa = 0.502500117621d0/35.16504*/
--		aa = 0.014289763856964,
--		bb = 0.057000649899720,
--
--		p0  =  2.570124672768757e-1,
--		p1  = -1.917742353032266e1,
--		p2  = -1.413382858617969e-2,
--		p3  = -5.427484830917552e-1,
--		p4  = -4.126621135193472e-4,
--		p5  = -4.176407833276121e-7,
--		p6  =  4.688217641883641e-5,
--		p7  = -3.039808885885726e-8,
--		p8  = -4.990118091261456e-11,
--		p9  = -9.733920711119464e-9,
--		p10 = -7.723324202726337e-12,
--		p11 =  7.121854166249257e-16,
--		p12 =  1.256474634100811e-12,
--		p13 =  2.105103897918125e-15,
--		p14 =  8.663811778227171e-19;
--
--	/*A very rough estimate of sa to get the saturated ct*/
--	if (ct != NULL) {
--	    sa = max(-(*ct + 9e-4*p)/0.06, 0.0);
--	    ctx = *ct;
--	} else if (t != NULL) {
--	    sa = max(-(*t + 9e-4*p)/0.06, 0.0);
--	    ctx = gsw_ct_from_t(sa,*t,p);
--	} else {
--	    return (0.0);
--	}
--	/*
--	! CTsat is the estimated value of CT if the seawater were saturated with
--	! dissolved air, recognizing that it actually has the air fraction
--	! saturation_fraction; see McDougall, Barker and Feistel, 2014).  
--	*/
--	ctsat = ctx - (1.0-saturation_fraction)*
--	        (1e-3)*(2.4-aa*sa)*(1.0+bb*(1.0-sa/gsw_sso));
--
--	return (p0 + p*(p2 + p4*ctsat + p*(p5 + ctsat*(p7 + p9*ctsat)
--	    + p*(p8  + ctsat*(p10 + p12*ctsat) + p*(p11 + p13*ctsat + p14*p))))
--	    + ctsat*(p1 + ctsat*(p3 + p6*p)));
--}

-- double
-- gsw_sa_freezing_from_ct(double ct, double p, double saturation_fraction)
{#fun unsafe gsw_sa_freezing_from_ct
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_freezing_from_ct_poly(double ct, double p, double saturation_fraction)
{#fun unsafe gsw_sa_freezing_from_ct_poly
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_freezing_from_t(double t, double p, double saturation_fraction)
{#fun unsafe gsw_sa_freezing_from_t
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_freezing_from_t_poly(double t, double p, double saturation_fraction)
{#fun unsafe gsw_sa_freezing_from_t_poly
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_from_rho(double rho, double ct, double p)
{#fun unsafe gsw_sa_from_rho
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_from_sp(double sp, double p, double lon, double lat)
{#fun unsafe gsw_sa_from_sp
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_from_sp_baltic(double sp, double lon, double lat)
{#fun unsafe gsw_sa_from_sp_baltic
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sa_from_sstar(double sstar, double p, double lon, double lat)
{#fun unsafe gsw_sa_from_sstar
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- int
-- gsw_sa_p_inrange(double sa, double p)
{#fun unsafe gsw_sa_p_inrange
        {`Double', `Double'} -> `Int' #}

-- void
-- gsw_seaice_fraction_to_freeze_seawater(double sa, double ct, double p,
--	double sa_seaice, double t_seaice, double *sa_freeze, double *ct_freeze,
--	double *w_seaice)
{#fun unsafe gsw_seaice_fraction_to_freeze_seawater
        {`Double', `Double', `Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_sigma0(double sa, double ct)
{#fun unsafe gsw_sigma0
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sigma1(double sa, double ct)
{#fun unsafe gsw_sigma1
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sigma2(double sa, double ct)
{#fun unsafe gsw_sigma2
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sigma3(double sa, double ct)
{#fun unsafe gsw_sigma3
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sigma4(double sa, double ct)
{#fun unsafe gsw_sigma4
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sound_speed(double sa, double ct, double p)
{#fun unsafe gsw_sound_speed
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sound_speed_ice(double t, double p)
{#fun gsw_sound_speed_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sound_speed_t_exact(double sa, double t, double p)
{#fun unsafe gsw_sound_speed_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sp_from_c(double c, double t, double p)
{#fun unsafe gsw_sp_from_c
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sp_from_sa(double sa, double p, double lon, double lat)
{#fun unsafe gsw_sp_from_sa
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sp_from_sa_baltic(double sa, double lon, double lat)
{#fun unsafe gsw_sp_from_sa_baltic
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sp_from_sk(double sk)
{#fun unsafe gsw_sp_from_sk
        {`Double'} -> `Double' #}


-- double
-- gsw_sp_from_sr(double sr)
{#fun unsafe gsw_sp_from_sr
        {`Double'} -> `Double' #}

-- double
-- gsw_sp_from_sstar(double sstar, double p, double lon, double lat)
{#fun unsafe gsw_sp_from_sstar
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_specvol(double sa, double ct, double p)
{#fun unsafe gsw_specvol
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_rho_alpha_beta (double sa, double ct, double p, double *rho, double *alpha,
--                       double *beta)
{#fun unsafe gsw_specvol_alpha_beta
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_specvol_anom_standard(double sa, double ct, double p)
{#fun unsafe gsw_specvol_anom_standard
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_specvol_first_derivatives (double sa, double ct, double p,
--				double *v_sa, double *v_ct, double *v_p)
{#fun unsafe gsw_specvol_first_derivatives
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_specvol_first_derivatives_wrt_enthalpy(double sa, double ct, double p,
--	double *v_sa, double *v_h)
{#fun unsafe gsw_specvol_first_derivatives_wrt_enthalpy
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_specvol_ice(double t, double p)
{#fun unsafe gsw_specvol_ice
        {`Double', `Double'} -> `Double' #}

-- void
-- gsw_specvol_second_derivatives (double sa, double ct, double p,
-- 	double *v_sa_sa, double *v_sa_ct, double *v_ct_ct, double *v_sa_p,
-- 	double *v_ct_p)
{#fun unsafe gsw_specvol_second_derivatives 
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_specvol_second_derivatives_wrt_enthalpy (double sa, double ct, double p,
-- 	double *v_sa_sa, double *v_sa_h, double *v_h_h)
{#fun unsafe gsw_specvol_second_derivatives_wrt_enthalpy
        {`Double', `Double', `Double',
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*,
        alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_specvol_sso_0(double p)
{#fun unsafe gsw_specvol_sso_0
        {`Double'} -> `Double' #}

-- double
-- gsw_specvol_t_exact(double sa, double t, double p)
{#fun unsafe gsw_specvol_t_exact
        {`Double', `Double', `Double'} -> `Double' #}


-- double
-- gsw_spiciness0(double sa, double ct)
{#fun unsafe gsw_spiciness0
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_spiciness1(double sa, double ct)
{#fun unsafe gsw_spiciness1
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_spiciness2(double sa, double ct)
{#fun unsafe gsw_spiciness2
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_sr_from_sp(double sp)
{#fun unsafe gsw_sr_from_sp
        {`Double'} -> `Double' #}

-- double
-- gsw_sstar_from_sa(double sa, double p, double lon, double lat)
{#fun unsafe gsw_sstar_from_sa
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_sstar_from_sp(double sp, double p, double lon, double lat)
{#fun unsafe gsw_sstar_from_sp
        {`Double', `Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_t_deriv_chem_potential_water_t_exact(double sa, double t, double p)
{#fun unsafe gsw_t_deriv_chem_potential_water_t_exact
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_t_freezing(double sa, double p, double saturation_fraction)
{#fun unsafe gsw_t_freezing
        {`Double', `Double', `Double'} -> `Double' #}

-- removed (https://github.com/TEOS-10/GSW-C/pull/12)
-- double
-- gsw_t_freezing_exact (double sa, double p, double saturation_fraction)
-- {#fun unsafe gsw_t_freezing_exact
--        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_t_freezing_first_derivatives(double sa, double p,
--	double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{#fun gsw_t_freezing_first_derivatives
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- void
-- gsw_t_freezing_first_derivatives_poly(double sa, double p,
--	double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{#fun unsafe gsw_t_freezing_first_derivatives_poly
        {`Double', `Double', `Double',
         alloca- `Double' peekFloatConv*,
         alloca- `Double' peekFloatConv*} -> `()' #}

-- double
-- gsw_t_freezing_poly(double sa, double p, double saturation_fraction)
-- 			int polynomial) -- removed (https://github.com/TEOS-10/GSW-C/pull/12)
{#fun unsafe gsw_t_freezing_poly
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_t_from_ct(double sa, double ct, double p)
{#fun unsafe gsw_t_from_ct
        {`Double', `Double', `Double'} -> `Double' #}

-- double
-- gsw_t_from_pt0_ice(double pt0_ice, double p)
{#fun unsafe gsw_t_from_pt0_ice
        {`Double', `Double'} -> `Double' #}

-- double
-- gsw_thermobaric(double sa, double ct, double p)
{#fun unsafe gsw_thermobaric
        {`Double', `Double', `Double'} -> `Double' #}

-- void
-- gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz,
--	double *tu, double *rsubrho, double *p_mid)
gsw_turner_rsubrho :: [Double] -> [Double] -> [Double] -> Int
                        -> IO ([Double], [Double], [Double])
gsw_turner_rsubrho sa ct p nz =
    withArrayConv sa $ \sa' ->
      withArrayConv ct $ \ct' ->
        withArrayConv p $ \p' ->
          do tu'      <- mallocArray nz :: IO (Ptr CDouble)
             rsubrho' <- mallocArray nz :: IO (Ptr CDouble)
             pMid'    <- mallocArray nz :: IO (Ptr CDouble)
             {#call unsafe gsw_turner_rsubrho as ^ #}
                    sa' ct' p' (fromIntegral nz) tu' rsubrho' pMid'
             tu       <- map realToFrac `fmap` peekArray nz tu'
             rsubrho  <- map realToFrac `fmap` peekArray nz rsubrho'
             pMid     <- map realToFrac `fmap` peekArray nz pMid'
             free tu'
             free rsubrho'
             free pMid'
             return (tu, rsubrho, pMid)

-- int
-- gsw_util_indx(double *x, int n, double z)
{#fun unsafe gsw_util_indx
        {withArrayConv* `[Double]',
                        `Int',
                        `Double'} -> `Int' #}

-- double *
-- gsw_util_interp1q_int(int nx, double *x, int *iy, int nxi, double *x_i,
-- 	double *y_i)
gsw_util_interp1q_int :: Int -> [Double] -> [Int] -> Int -> [Double]
                            -> IO (Maybe [Double])
gsw_util_interp1q_int nx x iy nxi xI =
    withArrayConv x $ \x' ->
      withArrayConvI iy $ \iy' ->
        withArrayConv xI $ \xI' ->
          do yI' <- mallocArray nxi :: IO (Ptr CDouble)
             ret <- {#call unsafe gsw_util_interp1q_int as ^ #}
                            (fromIntegral nx) x' iy' (fromIntegral nxi) xI' yI'
             if ret == nullPtr
                then free yI' >> return Nothing
                else do yI <- map realToFrac `fmap` peekArray nxi yI'
                        free yI'
                        return $ Just yI

----
---- No entry available in gswteos-10.h
----
------ double *
------ gsw_util_linear_interp(int nx, double *x, int ny, double *y, int nxi,
------	double *x_i, double *y_i)
----gsw_util_linear_interp :: Int -> [Double] -> Int -> [Double] -> Int -> [Double]
----                            -> IO (Maybe [Double])
----gsw_util_linear_interp nx x ny y nxi xI =
----    withArrayConv x $ \x' ->
----      withArrayConvI y $ \y' ->
----        withArrayConv xI $ \xI' ->
----          do yI' <- mallocArray nxi :: IO (Ptr CDouble)
----             ret <- {#call unsafe gsw_util_linear_interp as ^ #}
----                        (fromIntegral nx) x' (fromIntegral ny) y' (fromIntegral nxi) xI' yI'
----             if ret == nullPtr
----                then free yI' >> return Nothing
----                else do yI <- map realToFrac `fmap` peekArray nxi yI'
----                        free yI'
----                        return $ Just yI
--
---- void
---- gsw_util_sort_real(double *rarray, int nx, int *iarray)
--gsw_util_sort_real :: [Double] -> Int
--                            -> IO [Double]
--gsw_util_sort_real rarray nx =
--    withArrayConv rarray $ \rarray' ->
--      do iarray' <- mallocArray nx :: IO (Ptr CInt)
--         {#call unsafe gsw_util_sort_real as ^ #}
--                rarray' (fromIntegral nx) iarray'
--         rarray'' <- map realToFrac `fmap` peekArray nx rarray'
--         free iarray'
--         return rarray''
--
---- double
---- gsw_util_xinterp1(double *x, double *y, int n, double x0)
--{#fun unsafe gsw_util_xinterp1
--        {withArrayConv* `[Double]',
--         withArrayConv* `[Double]',
--                        `Int',
--                        `Double'} -> `Double' #}
--
---- double
---- gsw_z_from_p(double p, double lat)
--{#fun unsafe gsw_z_from_p
--        {`Double', `Double'} -> `Double' #}
--
---- double
---- gsw_saar(double p, double lon, double lat)
--{#fun unsafe gsw_saar
--        {`Double', `Double', `Double'} -> `Double' #}
--
---- double
---- gsw_deltasa_atlas(double p, double lon, double lat)
--{#fun unsafe gsw_deltasa_atlas
--        {`Double', `Double', `Double'} -> `Double' #}
