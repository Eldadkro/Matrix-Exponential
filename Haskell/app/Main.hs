module Main where
import Prelude hiding((<>))
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
phi13 = 5.371920351148152
phiPairs = [(3,0.01495585217958292),(5,0.2539398330063230),(7,0.9504178996162932),(9,2.097847961257068)]
bVector =  vector[1,182,16380,960960,40840800,1323241920,33522128640,670442572800,10559470521600,129060195264000,1187353796428800,7771770303897600,32382376266240000,64764752532480000]

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
calcPadeMU a t | a == 0 = multiplieMyMatrix (a + 1) t
               | otherwise = scale (bVector ! (fromIntegral a + 1))  (multiplieMyMatrix (a+1) t) + calcPadeMU (a-1) t

calcPadeMV :: Integer -> Matrix Double -> Matrix Double
calcPadeMV a t | a == 0 = multiplieMyMatrix a t
               | otherwise = scale (bVector ! fromIntegral a )  (multiplieMyMatrix a t) + calcPadeMV (a-1) t

calcPadeM :: Integer -> Matrix Double -> (Matrix Double , Matrix Double)
calcPadeM a t | a == -1 = (calcPadeMU (13 `div` 2) (scale s t), calcPadeMV (13 `div` 2) (scale s t))
              | otherwise =  (calcPadeMU (a `div` 2) t, calcPadeMV (a `div` 2) t)
              where 
                s = fromInteger intS
                intS = floor logResult
                logResult = logBase 2 ( norm_1 t / phi13)

main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"
    let ourM = checkNorm phiPairs mat
    let answer = calcPadeM ourM mat
    print answer

    let b = fst answer + snd answer
    let a = snd answer - fst answer
    print (linearSolve a b)
