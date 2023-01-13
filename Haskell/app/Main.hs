module Main where
import Prelude hiding((<>))
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Data.Maybe
import Text.ParserCombinators.ReadP (skipMany1)
phi13 = 5.371920351148152
phiPairs = [(3,0.01495585217958292),(5,0.2539398330063230),(7,0.9504178996162932),(9,2.097847961257068)]
bVector =  vector[64764752532480000, 32382376266240000, 7771770303897600,1187353796428800, 129060195264000, 10559470521600,670442572800, 33522128640, 1323241920,40840800, 960960, 16380, 182, 1]

matrixSize = 1000

checkNorm :: [(Integer,Double)] -> Matrix Double -> Integer
checkNorm [] _ = -1
checkNorm [x] amatrix | snd x > norm_1 amatrix = fst x 
                | otherwise = -1
checkNorm (x:xs) amatrix | snd x > norm_1 amatrix = fst x 
                   | otherwise = checkNorm xs amatrix


multiplieMyMatrix :: Integer -> Matrix Double -> Matrix Double
multiplieMyMatrix times amatrix | times == 0 = ident matrixSize
                      | times == 1 = amatrix
                      | otherwise = amatrix <> multiplieMyMatrix (times-1) amatrix

calcPadeMU :: Integer -> Matrix Double -> Matrix Double
calcPadeMU index amatrix | index == 0 = scale b2p1 mat2p1
               | otherwise = scale b2p1 mat2p1 + calcPadeMU (index-1) amatrix
               where 
                  mat2p1 = multiplieMyMatrix k2p1 amatrix
                  b2p1 = bVector ! fromIntegral k2p1
                  k2p1 = (2 * index) + 1
                  

calcPadeMV :: Integer -> Matrix Double -> Matrix Double
calcPadeMV index amatrix | index == 0 = scale b2 (multiplieMyMatrix 0  amatrix)
               | otherwise = scale b2 mat2 + calcPadeMV (index-1) amatrix
                where  
                  mat2 = multiplieMyMatrix k2 amatrix
                  b2 = bVector ! fromIntegral k2
                  k2 = index * index

calcPadeM :: Integer -> Integer ->  Matrix Double -> (Matrix Double , Matrix Double)
calcPadeM m s amatrix | m == -1 = (calcPadeMU (13 `div` 2) (scale toScale  amatrix), calcPadeMV (13 `div` 2) (scale toScale  amatrix))
              | otherwise =  (calcPadeMU (m `div` 2) amatrix, calcPadeMV (m `div` 2)  amatrix)
              where 
                toScale = 1 / 2^s

getS :: Matrix Double -> Integer
getS amatrix = s
        where 
        s = ceiling logResult
        logResult = logBase 2 ( norm_1 amatrix / phi13)

scaleBack :: Integer -> Matrix Double -> Matrix Double
scaleBack s amatrix | s == 0 = amatrix
              | otherwise = sm1<>sm1
              where
                sm1 = scaleBack (s-1) amatrix

shouldScaleBack :: Integer -> Integer -> Matrix Double -> Matrix Double
shouldScaleBack m s amatrix | m == -1 = scaleBack s amatrix
                            | otherwise =  amatrix

--kaczmarzStart :: Matrix Double -> Matrix Double -> Matrix Double
--kaczmarzStart a s amatrix | a == -1 = scaleBack s amatrix
--                     | otherwise =  amatrix

main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"
    let ourM = checkNorm phiPairs mat
    let ourS = getS mat
    print ourS
    let answer = calcPadeM ourM ourS mat
    let u = fst answer
    let v = snd answer
    saveMatrix "haskell-Utest.txt" "%0.10f" u
    saveMatrix "haskell-Vtest.txt" "%0.10f" v
    let p = u + v
    let q = v - u
    saveMatrix "haskell-Qtest.txt" "%0.10f" q
    saveMatrix "haskell-Ptest.txt" "%0.10f" p
    let answer2 = linearSolve q p 
    let answer3 = fromJust answer2
    saveMatrix "beforeFinalScale.txt" " %.20f ," answer3

--    let solutionKac = 
    let a2 = answer3<>answer3
    let a4 = a2<>a2
    let a8 = a4<>a4
    let a16 = a8<>a8
    let a32 = a16<>a16
    --print answer2
    let finalAnswer = shouldScaleBack ourM ourS answer3
    --print finalAnswer
    saveMatrix "answerDirectSquare.txt" " %0.20f ," a32
    saveMatrix "answer--InDirectSquare.txt" " %0.20f ," finalAnswer
