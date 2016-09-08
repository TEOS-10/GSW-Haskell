{-# LANGUAGE TypeFamilies, ScopedTypeVariables #-}

--
-- | Not all functions of GSW in matlab are not exported
--   to fortran, i.e., FFI is not a full solution.
--
module Oceanogr.GSW (
-- * distance calculated in both Double and Float
    LatLon, Loc(..), getLon, getLat, putLoc,
--
    gsw_distance, gsw_f,
    linearInterp1, interpRRD, interp1
           )
where

import Control.Monad (when)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Numeric.IEEE (nan)

class Floating a => LatLon a where
    data Loc a :: *
    getLon, getLat :: Loc a -> a
    putLoc  :: a -> a -> Loc a

instance LatLon Double where
    data Loc Double = LocD {getLonD  :: !Double,  getLatD :: !Double} deriving (Eq, Show)
    getLon = getLonD
    getLat = getLatD
    putLoc = LocD 

instance LatLon Float where
    data Loc Float  = LocF {getLonF :: !Float,   getLatF :: !Float} deriving (Eq, Show)
    getLon = getLonF
    getLat = getLatF
    putLoc = LocF

-- no depth correction
-- based on Matlab v3.4
gsw_distance :: forall a. (RealFloat a, LatLon a) => Loc a -> Loc a -> a
gsw_distance p0 p1 = earth_radius * angles
  where
    dlon = pi180 * (getLon p0 - getLon p1)
    dlat = pi180 * (getLat p0 - getLat p1)
    a = sin (0.5 * dlat)^(2::Integer)
      + cos (getLat p0 * pi180) * cos (getLat p1 * pi180) * sin (0.5 * dlon)^(2::Integer)
    angles = 2 * atan2 (sqrt a) (sqrt $ 1.0-a)

    pi180 = 3.14159265358979 / 180.0
    earth_radius = 6371000 -- [m]

---
gsw_f :: Floating a => a -> a
gsw_f lat = let deg2rad = pi / 180
                omega   = 7.292115e-5
             in 2 * omega * sin (lat * deg2rad)

---
-- |
--  Combination of interpRRD and interp1
--  First, apply interpRRD, second nan's are
--  filled by interp1
--
-- >>> let x = U.fromList [2020,2216,2413,2611,2878,3000] :: U.Vector Double
-- >>> let y = U.fromList [2.211,2.011,1.894,1.788,1.554,1.380] :: U.Vector Double
-- >>> let xi = U.fromList [2100,2200 .. 3000] :: U.Vector Double
-- >>> interp1 x y xi
-- >>> Right (fromList [2.1293673469387753,2.027326530612245,1.9570207317973867,1.9013255208670579,1.8455641980510156,1.7947405764986781,1.7292355861090267,1.6473991803675845,1.5226229508196722,1.38])
--
--
interp1 :: U.Vector Double -- ^ bottle pressure
        -> U.Vector Double -- ^ bottle data
        -> U.Vector Double -- ^ pressure to interpolate onto
        -> IO (Either String (U.Vector Double)) -- ^ output

interp1 x y xi = do

    y0 <- interpRRD x y xi

    case y0 of
      Left e0  -> return $ Left e0
      Right y1 -> do
        let idx = U.fromList [0 .. (U.length xi - 1)] :: U.Vector Int
            gap = U.filter (\(_,_,y10) -> isNaN y10) $ U.zip3 idx xi y1
        if U.null gap
        then return $ Right y1
        else do
            let (ig, xg, _) = U.unzip3 gap
            yg <- linearInterp1 x y xg
            case yg of
              Left e1  -> return $ Left e1
              Right y2 -> let toUpdate = U.zip ig y2
                           in return $ Right (U.update y1 toUpdate)


---
-- | Interpolation of "bottle" data
--   optimized for pelagic ocean data
--   based on the method by Reinger and Ross (1968, DSR)
--   tuned by Jeff Dunn.
--
--   Matlab source is in  ..\/matlabGSW3.04\/gsw_geo_strf_dyn_height.m
--   Example below was taken from example.dat in the gamma_n package.
--   observed temperatures at xi are; 12.21,11.99,10.54,8.36,7.43,6.04
--   observed salinities   at xi are; 35.086,35.078,34.851,34.572,34.509,34.452
--
-- >>> let x = U.fromList [1,97,194,388,581,775,969] :: U.Vector Double
-- >>> let ytemp = U.fromList [12.25,12.09,11.69,9.35,7.86,6.87,5.5] :: U.Vector Double
-- >>> let ysalt = U.fromList [35.066,35.089,35.025,34.696,34.531,34.496,34.458] :: U.Vector Double
-- >>> let xi = U.fromList [48,145,291,485,678,872] :: U.Vector Double
-- >>> interpRRD x ytemp xi
-- Right (fromList [NaN,NaN,10.466597777964449,8.45058039567543,7.375825606234917,NaN])
-- >>> interpRRD x ysalt xi
-- Right (fromList [NaN,NaN,34.85724213411662,34.57730605075967,34.51419549586589,NaN])
--
interpRRD :: U.Vector Double -- ^ bottle pressure
          -> U.Vector Double -- ^ bottle data
          -> U.Vector Double -- ^ pressure to interpolate onto
          -> IO (Either String (U.Vector Double)) -- ^ output
interpRRD x y xi
    | U.length x < 4 = return $ Left "Input too short"
    | U.length x /= U.length y
                    = return $ Left "Inconsistent length"
    | not (U.all (>0) (U.zipWith (-) (U.tail x) x)
        || U.all (<0) (U.zipWith (-) (U.tail x) x)) = return $ Left "Unsorted input"
    | otherwise      = do
      let xfn = U.fromList [0, 300, 1800, 8000] :: U.Vector Double
          yfn = U.fromList [50, 200, 650, 1250] :: U.Vector Double

          interp1' :: Double -> Double
          interp1' x0 = case U.findIndex (x0<) xfn of
                         Nothing -> error $ "too deep at x=" ++ show x0
                         Just i0 -> (x0 - xfn U.! (i0-1)) * (yfn U.! i0 - yfn U.! (i0-1))
                                                         / (xfn U.! i0 - xfn U.! (i0-1))
                                       + yfn U.! (i0-1)
          maxDis = U.map interp1' xi
          ilen   = U.length xi
          ivec   = U.fromList [0 .. (ilen-1)] :: U.Vector Int

      yi <- UM.replicate ilen nan :: IO (UM.IOVector Double)

      let coincide :: (Double, Int) -> IO ()
          coincide (x0, j) =
            let i0 = U.minIndexBy (\a b -> abs (a - x0) `compare` abs (b - x0)) x
             in when (abs (x U.! i0 - x0) < U.minimum maxDis * 0.005) -- "concid_frac"
                   $ UM.write yi j (y U.! i0)

      U.mapM_ coincide $ U.zip xi ivec

    
      let interp :: (Double, Double, Int) -> IO ()
          interp (y0, x0, j)
            | not . isNaN $ y0 = return () -- already assigned
            | otherwise        = case inbetween x x0 of
                        Nothing -> return () -- not within the range
                        Just i0 -> if i0 < 1 || i0 > U.length x - 3
                                   then return () -- not within the range
                                   else do
                                   let sumin = abs $ x U.! i0 - x U.! (i0+1)
                                       outl  = abs $ x U.! (i0-1) - x U.! i0
                                       outr  = abs $ x U.! (i0+2) - x U.! (i0+1)
                                       md    = maxDis U.! j
                                   if sumin > md || outl+sumin+outr > md*3
                                        || min outr outl < sumin / 15.0
                                   -- sumin>maxdis(tdix,1) | outtest>maxdis(tidx,2) |
                                   --      minsep<cfrac
                                   then return ()
                                   else let intp = y U.! i0 + (x0 - x U.! i0)
                                                       * (y U.! (i0+1) - y U.! i0)
                                                       / (x U.! (i0+1) - x U.! i0)
                                            extl = y U.! i0 + (x0 - x U.! i0)
                                                       * (y U.! i0 - y U.! (i0-1))
                                                       / (x U.! i0 - x U.! (i0-1))
                                            extr = y U.! (i0+1) + (x0 - x U.! (i0+1))
                                                       * (y U.! (i0+2) - y U.! (i0+1))
                                                       / (x U.! (i0+2) - x U.! (i0+1))
                                            deno = abs (extl - intp) ** 1.7
                                                 + abs (intp - extr) ** 1.7
                                            nume = abs (extl - intp) ** 1.7 * extr
                                                 + abs (intp - extr) ** 1.7 * extl
                                            y9'' = if abs deno < 1.0e-4
                                                   then intp
                                                   else (intp + nume / deno) * 0.5
                                            -- do not overshoot
                                            y9'  = max y9'' $ min (y U.! i0) (y U.! (i0+1))
                                            y9   = min y9'  $ max (y U.! i0) (y U.! (i0+1))
                                         in UM.write yi j y9
      yi' <- U.freeze yi
      U.mapM_ interp $ U.zip3 yi' xi ivec
      Right `fmap` U.unsafeFreeze yi

--
-- | find the coodinates surrounding here on ruler
--
inbetween :: U.Vector Double -> Double -> Maybe Int
inbetween ruler here = let ruler' = U.map (subtract here) ruler
                        in if U.head ruler' * U.last ruler' > 0.0
                             then Nothing
                             else U.findIndex (\(p,q) -> p*q <= 0.0) $ U.zip ruler' (U.tail ruler')
 
--
-- | Simple linear interpolation
--
-- >>> let x = U.fromList [1,97,194,388,581,775,969]               :: U.Vector Double
-- >>> let y = U.fromList [12.25,12.09,11.69,9.35,7.86,6.87,5.5]  :: U.Vector Double
-- >>> let xi = U.fromList [48,145,291,485,678,872]              :: U.Vector Double
-- >>> linearInterp1 x y xi
-- Right (fromList [12.171666666666667,11.892061855670104,10.52,8.601139896373057,7.365,6.1850000000000005])
--
linearInterp1 :: U.Vector Double -- ^ bottle pressure
        -> U.Vector Double -- ^ bottle data
        -> U.Vector Double -- ^ pressure to interpolate onto
        -> IO (Either String (U.Vector Double)) -- ^ output
linearInterp1 x y xi
    | U.length x < 2 = return $ Left "Input too short"
    | U.length x /= U.length y
                    = return $ Left "Inconsistent length"
    | not (U.all (>0) (U.zipWith (-) (U.tail x) x)
        || U.all (<0) (U.zipWith (-) (U.tail x) x)) = return $ Left "Unsorted input"
    | otherwise      = do
        yi <- UM.replicate (U.length xi) nan :: IO (UM.IOVector Double)
        let interp :: (Double, Int) -> IO ()
            interp (x0, j) = case inbetween x x0 of
                        Nothing -> return () -- not within the range
                        Just i0 -> UM.write yi j $ y U.! i0 + (x0 - x U.! i0)
                                                    * (y U.! (i0+1) - y U.! i0)
                                                    / (x U.! (i0+1) - x U.! i0)
        U.mapM_ interp $ U.zip xi (U.fromList [0 .. (U.length xi - 1)])
        Right `fmap` U.unsafeFreeze yi
