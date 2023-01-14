module Main where
import Prelude hiding((<>))
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Data.Maybe

phi13 = 5.371920351148152
phiPairs = [(3,0.01495585217958292),(5,0.2539398330063230),(7,0.9504178996162932),(9,2.097847961257068)]
bVector =  vector[64764752532480000, 32382376266240000, 7771770303897600,1187353796428800, 129060195264000, 10559470521600,670442572800, 33522128640, 1323241920,40840800, 960960, 16380, 182, 1]

matrixSize = 1000

-- check matrix norm
checkNorm :: [(Integer,Double)] -> Matrix Double -> Integer
checkNorm [] _ = -1
checkNorm [x] t | snd x > norm_1 t = fst x 
                | otherwise = -1
checkNorm (x:xs) t | snd x > norm_1 t = fst x 
                   | otherwise = checkNorm xs t

-- mutiplie matrix can be optimized by having it ran from 0 to the int while passing the previous matrix forword
multiplieMyMatrix :: Integer -> Matrix Double -> Matrix Double
multiplieMyMatrix times myMatrix | times == 0 = ident matrixSize
                      | times == 1 = myMatrix
                      | otherwise = myMatrix <> multiplieMyMatrix (times-1) myMatrix

--calculate the U part of pade 
calcPadeMU :: Integer -> Matrix Double -> Matrix Double
calcPadeMU degree myMatrix | degree == 0 = scale b2p1 mat2p1
               | otherwise = scale b2p1 mat2p1 + calcPadeMU (degree-1) myMatrix
               where 
                  mat2p1 = multiplieMyMatrix k2p1 myMatrix
                  b2p1 = bVector ! fromIntegral k2p1
                  k2p1 = (2 * degree) + 1
                  
--calculate the V part of pade 
calcPadeMV :: Integer -> Matrix Double -> Matrix Double
calcPadeMV degree myMatrix | degree == 0 = scale b2 (multiplieMyMatrix 0 myMatrix)
               | otherwise = scale b2 mat2 + calcPadeMV (degree-1) myMatrix
                where  
                  mat2 = multiplieMyMatrix k2 myMatrix
                  b2 = bVector ! fromIntegral k2
                  k2 = 2 * degree

-- calculate pade of degree M returns a pair of matrices (U,V)
calcPadeM :: Integer -> Integer ->  Matrix Double -> (Matrix Double , Matrix Double)
calcPadeM degree s myMatrix | degree == -1 = (calcPadeMU pade13Index scaledMatrix, calcPadeMV pade13Index scaledMatrix)
              | otherwise =  (calcPadeMU padeIndex myMatrix, calcPadeMV padeIndex myMatrix)
              where 
                toScale = 1 / 2^s
                padeIndex = degree `div` 2
                pade13Index = 13 `div` 2
                scaledMatrix = scale toScale myMatrix

-- calculate the scaling factor 
getS :: Matrix Double -> Integer
getS myMatrix = s
        where 
        s = ceiling logResult
        logResult = logBase 2 ( norm_1 myMatrix / phi13)

-- scale back function
scaleBack :: Integer -> Matrix Double -> Matrix Double
scaleBack s myMatrix | s == 0 = myMatrix
              | otherwise = sm1<>sm1
              where
                sm1 = scaleBack (s-1) myMatrix

shouldScaleBack m s myMatrix | m == -1 = scaleBack s myMatrix
                            | otherwise =  myMatrix

numberOfSweeps = 100
-- run kaczmarz sweep
runKaczmarzSweeps :: Matrix Double -> Matrix Double -> Matrix Double -> Int -> Matrix Double
runKaczmarzSweeps a b x i | i >= numberOfSweeps = x
                          | otherwise = runKaczmarzSweeps a b newX (i+1)
                          where
                              newX = fromColumns (kaczmarzForEach a b x matrixSize 0)

kaczmarzLoopNumber = 1000
--for each column of x perform kaczmarz
kaczmarzForEach :: Matrix Double -> Matrix Double -> Matrix Double -> Int -> Int -> [Vector Double]
kaczmarzForEach a b x n i | i >= n = []
                          | otherwise = kacI  : kaczmarzForEach a b x n (i+1)
                          where
                              kacI = kaczmarz' a bi xi kaczmarzLoopNumber 0
                              xi = flatten (x ¿ [i])
                              bi = flatten (b ¿ [i])
-- perform kaczmarz algorithm                               
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
                          kaczmarzAnswer = runKaczmarzSweeps q p x 0
                          finalAnswer = shouldScaleBack ourM ourS kaczmarzAnswer


main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"

    let ourAnswer = padeApproxemteLS mat

    let haskellSol = expm mat

    let dif = norm_1 (ourAnswer - haskellSol)
    print dif
    

