
    let ourM = checkNorm phiPairs mat
    let ourS = getS mat
    let answer = calcPadeM ourM ourS mat
    let u = fst answer
    let v = snd answer
    let p = u + v
    let q = v - u

    --Kaczmarz method

    --randx <- randn 1000 1000
    let randx = scale 0 (ident matrixSize)
    
    --result <- kaczmarz qMatrix pMatrix x n
    let columnsX = toColumns randx
    --let listXCols = kaczmarz2 0 rowsA columnsB columnsX
    --let answeraa = kaczmarz q 
    let answeraa = runKaczmarzLoop q p randx 0
    --saveMatrix "test22Kaz.txt" "%0.10f" answeraa
    let finalAnswerkaz = shouldScaleBack ourM ourS answeraa

    
    let answer2 = linearSolve q p 
    let answer23 = fromJust answer2
    --saveMatrix "beforeScaleLS.txt" " %f ," answer23
    --print answer2
    let finalAnswer = shouldScaleBack ourM ourS answer23
    --print finalAnswer
    let haskellSol = expm mat
    let dif = norm_1 (finalAnswerkaz - haskellSol)
    print dif
    --saveMatrix "answer.txt" " %f ," finalAnswer
