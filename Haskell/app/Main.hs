module Main where
import Prelude hiding((<>))
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Data.Maybe

import Data.Time


phi13 = 5.371920351148152
phiPairs = [(3,0.01495585217958292),(5,0.2539398330063230),(7,0.9504178996162932),(9,2.097847961257068)]
bVector =  vector[64764752532480000, 32382376266240000, 7771770303897600,1187353796428800, 129060195264000, 10559470521600,670442572800, 33522128640, 1323241920,40840800, 960960, 16380, 182, 1]

matrixSize = 1000

checkNorm :: [(Integer,Double)] -> Matrix Double -> Integer
checkNorm [] _ = -1
checkNorm [x] t | snd x > norm_1 t = fst x 
                | otherwise = -1
checkNorm (x:xs) t | snd x > norm_1 t = fst x 
                   | otherwise = checkNorm xs t


multiplieMyMatrix :: Integer -> Matrix Double -> Matrix Double
multiplieMyMatrix a t | a == 0 = ident matrixSize
                      | a == 1 = t
                      | otherwise = t <> multiplieMyMatrix (a-1) t

calcPadeMU :: Integer -> Matrix Double -> Matrix Double
calcPadeMU a t | a == 0 = scale b2p1 mat2p1
               | otherwise = scale b2p1 mat2p1 + calcPadeMU (a-1) t
               where 
                  mat2p1 = multiplieMyMatrix k2p1 t
                  b2p1 = bVector ! fromIntegral k2p1
                  k2p1 = (2 * a) + 1
                  

calcPadeMV :: Integer -> Matrix Double -> Matrix Double
calcPadeMV a t | a == 0 = scale b2 (multiplieMyMatrix 0 t)
               | otherwise = scale b2 mat2 + calcPadeMV (a-1) t
                where  
                  mat2 = multiplieMyMatrix k2 t
                  b2 = bVector ! fromIntegral k2
                  k2 = 2 * a

calcPadeM :: Integer -> Integer ->  Matrix Double -> (Matrix Double , Matrix Double)
calcPadeM a s t | a == -1 = (calcPadeMU (13 `div` 2) (scale toScale t), calcPadeMV (13 `div` 2) (scale toScale t))
              | otherwise =  (calcPadeMU (a `div` 2) t, calcPadeMV (a `div` 2) t)
              where 
                toScale = 1 / 2^s

getS :: Matrix Double -> Integer
getS t = s
        where 
        s = ceiling logResult
        logResult = logBase 2 ( norm_1 t / phi13)

scaleBack :: Integer -> Matrix Double -> Matrix Double
scaleBack s amatrix | s == 0 = amatrix
              | otherwise = sm1<>sm1
              where
                sm1 = scaleBack (s-1) amatrix

shouldScaleBack m s amatrix | m == -1 = scaleBack s amatrix
                            | otherwise =  amatrix




kaczmarzLoopNumber = 1000
totalLoops = 100
runKaczmarzLoop :: Matrix Double -> Matrix Double -> Matrix Double -> Int -> Matrix Double
runKaczmarzLoop a b x i | i >= totalLoops = x
                          | otherwise = runKaczmarzLoop a b newX (i+1)
                          where
                              newX = fromColumns (kaczmarzForEach a b x matrixSize 0)

kaczmarzForEach :: Matrix Double -> Matrix Double -> Matrix Double -> Int -> Int -> [Vector Double]
kaczmarzForEach a b x n i | i >= n = []
                          | otherwise = kacI  : kaczmarzForEach a b x n (i+1)
                          where
                              kacI = kaczmarz' a bi xi kaczmarzLoopNumber 0
                              xi = flatten (x ¿ [i])
                              bi = flatten (b ¿ [i])
                              



changeV2L :: Int -> [Vector Double] -> [[Double]]
changeV2L index x | index < (matrixSize - 1) = myloop
                  | otherwise = [toList (x!!index)]
                  where 
                    myloop = toList (x!!index) : changeV2L (index+1) x

kaczmarz :: Matrix Double -> Vector Double -> Vector Double -> Int -> Vector Double
kaczmarz a b x n = kaczmarz' a b x n 0

kaczmarz' :: Matrix Double -> Vector Double -> Vector Double -> Int -> Int -> Vector Double
kaczmarz' a b x n i
    | i >= n = x
    | otherwise = kaczmarz' a b x' n (i+1)
    where
        r = a ! i -- get i-th row of A
        dot = r <.> x -- calculate the dot product of the row and x
        c = (b ! i - dot) / (r <.> r) -- calculate the correction term
        x' = x + scale c r -- update x

padeApproxemteLS :: Matrix Double -> Matrix Double
padeApproxemteLS aMatrix = finalAnswer 
                        where 
                          ourM = checkNorm phiPairs aMatrix
                          ourS = getS aMatrix
                          padePolynom = calcPadeM ourM ourS aMatrix
                          u = fst padePolynom
                          v = snd padePolynom
                          p = u + v
                          q = v - u
                          n =  1000 :: Int 
                          zero = 0 :: Double
                          x = konst zero (n,n)
                          kaczmarzAnswer = runKaczmarzLoop q p x 0
                          finalAnswer = shouldScaleBack ourM ourS kaczmarzAnswer


main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"

    let ourAnswer = padeApproxemteLS mat

    let haskellSol = expm mat

    let dif = norm_1 (ourAnswer - haskellSol)
    print dif
    

