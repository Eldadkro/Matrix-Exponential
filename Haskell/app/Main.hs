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
calcPadeMU a t | a == 0 = scale b2p1 mat2p1
               | otherwise = scale b2p1 mat2p1 + calcPadeMU (a-1) t
               where 
                  mat2p1 = multiplieMyMatrix k2p1 t
                  b2p1 = bVector ! fromIntegral k2p1
                  k2p1 = (2 * a) + 1
                  

calcPadeMV :: Integer -> Matrix Double -> Matrix Double
calcPadeMV a t | a == 0 = multiplieMyMatrix 0 t
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

shouldScaleBack :: Integer -> Integer -> Matrix Double -> Matrix Double
shouldScaleBack a s t | a == -1 = t^toScale
              | otherwise =  t
              where 
                toScale = 2^s 

main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"
    let ourM = checkNorm phiPairs mat
    let ourS = getS mat
    let answer = calcPadeM ourM ourS mat
    let u = fst answer
    let v = snd answer
    let b = u + v
    let a = v - u
    let answer2 = linearSolveLS a b
    let finalAnswer = shouldScaleBack ourM ourS answer2
    saveMatrix "answer.txt" " %f ," finalAnswer
