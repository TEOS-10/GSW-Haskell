module Main where
import Oceanogr.GSW
import Oceanogr.GSWtools

import Control.Monad as C
import Data.Ord (comparing)
import Data.List (zip3, unzip3, zip4)
import qualified Data.Vector.Unboxed as U
import Text.Printf

main :: IO ()
main = do

    let matOut, matOutIntrp :: [Double]
        matOut = [24.935470, 24.206508, 23.452009, 22.715286, 21.962051, 20.496888,
                  19.089137, 17.738121, 16.435467, 15.148141, 13.901794, 12.721671,
                  11.607833, 10.558857, 8.6406232, 6.9477307, 5.4534273, 4.0812897,
                  2.9007551, 1.6775919, 0.50250733, (-0.62856791), (-2.0876251), (-2.7151673)]
        matOutIntrp = [23.405964,21.871422,20.366272,18.922001,17.534582,16.183310,
                       14.865464,13.597639,12.400151,11.272593,10.222672,9.2338914,
                       8.2915732,7.4189482,6.5915755,5.8173250,5.0785649,4.3712775,
                       3.6921749,3.0326584,2.4015048,1.7774420,1.1765406,0.58005098,
                       0.0055197216,(-0.56573040),(-1.1149203),(-1.6613837)]
        mpres = [100, 200 .. 2800] :: [Double]

    c <- readFile "Test/example.dat"

    let ls = lines c
        ws0 = words . head $ ls
        lon = read (ws0 !! 0) :: Double
        lat = read (ws0 !! 1) :: Double
        n   = read (ws0 !! 2) :: Int
        (ss, ts, ps) = unzip3 $ map (\l -> let ws = words l
                                               s  = read (ws !! 0) :: Double
                                               t  = read (ws !! 1) :: Double
                                               p  = read (ws !! 2) :: Double
                                            in (s,t,p)) $ tail ls -- skip the header line
        sa = U.fromList ss :: U.Vector Double
        ct = U.fromList ts :: U.Vector Double
        pr = U.fromList ps :: U.Vector Double
        pi'= U.fromList mpres :: U.Vector Double

    g <- gsw_geo_strf_dyn_height (U.toList sa) (U.toList ct) (U.toList pr) 2500.0 (U.length pr)

    case g of
        Nothing -> error "gsw_geo_strf_dyn_height: error"
        Just g1 -> do printf "Uninterpolated\n"
                      printf "Depth     Haskell    Matlab\n"
                      C.forM_ [0 .. (U.length pr - 1)] $ \i ->
                        printf "%10.0f %10.4f %10.4f\n" (pr U.! i) (g1 !! i) (matOut !! i)
   
                      gintp' <- interp1 pr (U.fromList g1) pi'
                      -- gintp' <- linearInterp1 pr (U.fromList g1) pi'
                      case gintp' of
                          Left e1 -> error e1
                          Right gintp -> do printf "Interpolated\n"
                                            printf "Depth     Haskell    Matlab\n"
                                            C.forM_ [0 .. (Prelude.length mpres - 1)] $ \i ->
                                              printf "%10.0f %10.4f %10.4f\n" (mpres !! i)
                                                              (gintp U.! i) (matOutIntrp !! i)
