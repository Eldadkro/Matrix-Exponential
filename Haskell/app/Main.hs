module Main where
import Prelude hiding((<>))
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import GHC.Plugins (lookupPackageName)
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

shouldScaleBack :: Integer -> Integer -> Matrix Double -> Matrix Double
shouldScaleBack a s t | a == -1 = t^toScale
              | otherwise =  t
              where 
                toScale = 2^s 


kaczmarzLoopNumber = 50

kaczmarzFunc :: Int -> [Vector Double] -> Vector Double -> Vector Double -> [Vector Double]
kaczmarzFunc index rowsA colB xj | index == matrixSize = [xj]
                                 | otherwise = kaczmarzFunc (index + 1) rowsA colB xjp1
                                    where 
                                      xjp1 = xj + flatten proj
                                      proj = scale ((bi - dotAiXjp1)/lenAi) aColVector
                                      bi = colB ! index
                                      dotAiXjp1 = (aRowVector * asColumn xj) `atIndex` (0,0)
                                      lenAi = (aRowVector * aColVector) `atIndex` (0,0)
                                      aColVector = asColumn (rowsA!!index)
                                      aRowVector = asRow (rowsA!!index)

kaczmarzForEach :: Int -> [Vector Double] -> [Vector Double] -> [Vector Double] -> [Vector Double]
kaczmarzForEach index rows columns x | index < matrixSize = myloop
                                     | otherwise = x
                                      where 
                                        myloop = xj ++ xRest
                                        xj = kaczmarzFunc 0 rows (columns!!index) (x!!index)
                                        xRest = kaczmarzForEach (index+1) rows columns x

                              

--kaczmarzForEach index rows = rows ++ kaczmarzFunc index

kaczmarzLoop :: Int -> [Vector Double] -> [Vector Double] -> [Vector Double] -> [Vector Double]
kaczmarzLoop index rowsA columnsB x | index == kaczmarzLoopNumber = x
                                    | otherwise = next 
                                    where 
                                      newX = kaczmarzForEach 0 rowsA columnsB x
                                      next = kaczmarzLoop (index+1) rowsA columnsB newX

changeV2L :: Int -> [Vector Double] -> [[Double]]
changeV2L index x | index < (matrixSize - 1) = myloop
                  | otherwise = [toList (x!!index)]
                  where 
                    myloop = toList (x!!index) : changeV2L (index+1) x
main :: IO ()
main = do

    mat <- loadMatrix "ExpoMatrix1000x1000.txt"
    let ourM = checkNorm phiPairs mat
    let ourS = getS mat
    let answer = calcPadeM ourM ourS mat
    let u = fst answer
    let v = snd answer
    let p = u + v
    let q = v - u

    --Kaczmarz method
    let rowsA = toRows q
    let columnsB = toColumns p
    randx <- randn 1000 1000
    let columnsX = toColumns randx
    let listXCols = kaczmarzLoop 0 rowsA columnsB columnsX
    print "aaa"
    writeFile "file2.txt" (show (changeV2L 0 listXCols))
    print "aaa"
    --writeFile "file.txt" (show listXCols)
    
    print "aaa"
    --let ourSol = fromRows listXCols
    --print (rows ourSol)
    --print (cols ourSol)
    print "bbb1"
    --saveMatrix "OurSolve" " %f ," ourSol
    print "bbb"
    --let ourSolfinalAnswer = shouldScaleBack ourM ourS ourSol
    --saveMatrix "OurSolve" " %f ," ourSolfinalAnswer

    let answer2 = linearSolveLS q p 

    --print answer2
    let finalAnswer = shouldScaleBack ourM ourS answer2
    --print finalAnswer
    saveMatrix "answer.txt" " %f ," finalAnswer
