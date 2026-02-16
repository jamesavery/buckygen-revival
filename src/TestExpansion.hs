module Main where

import Seeds
import Expansion
import qualified Data.IntMap.Strict as IM
import Data.List (sort, nub, foldl')

---------------------------------------------------------------------
-- Test infrastructure
---------------------------------------------------------------------

data TestResult = Pass String | Fail String String
    deriving (Show)

passed :: TestResult -> Bool
passed (Pass _) = True
passed _        = False

runTests :: [TestResult] -> IO ()
runTests results = do
    let total = length results
        passes = length (filter passed results)
        fails  = filter (not . passed) results
    mapM_ (\r -> case r of
        Pass name -> putStrLn $ "  PASS: " ++ name
        Fail name msg -> putStrLn $ "  FAIL: " ++ name ++ " -- " ++ msg
        ) results
    putStrLn $ "\n" ++ show passes ++ "/" ++ show total ++ " tests passed."
    if null fails
        then putStrLn "All tests passed!"
        else do
            putStrLn $ show (length fails) ++ " FAILURES:"
            mapM_ (\(Fail n m) -> putStrLn $ "  " ++ n ++ ": " ++ m) fails

---------------------------------------------------------------------
-- CW rotation system validation (face computation)
---------------------------------------------------------------------

-- | Compute all faces of the embedding defined by the CW rotation system.
-- Traces the face to the right of each dart (u -> v): the next dart in
-- the face is (v, nextCW g v u).
-- Returns list of faces, where each face is a list of vertices.
computeFaces :: DualGraph -> [[Int]]
computeFaces g = go allDarts [] usedEmpty
  where
    nv = numVertices g
    allDarts = [(u, v) | u <- [0..nv - 1], v <- nbrs g u]

    -- Track used darts as a set of (u,v) pairs via IntMap of lists
    usedEmpty = IM.empty :: IM.IntMap [Int]
    usedMember (u,v) s = case IM.lookup u s of
        Nothing -> False
        Just vs -> v `elem` vs
    usedInsert (u,v) s = IM.insertWith (++) u [v] s
    usedInsertAll darts s = foldl' (\acc d -> usedInsert d acc) s darts

    go [] faces _ = faces
    go (dart:rest) faces used
        | usedMember dart used = go rest faces used
        | otherwise =
            let face = traceFace dart
                ds = dartsOfFace face
            in go rest (face : faces) (usedInsertAll ds used)

    -- Trace the face to the right of dart (u0, v0).
    -- Collects vertices: at dart (u, v), collect u, advance to (v, nextCW g v u).
    -- Stops when returning to the starting dart (u0, v0).
    traceFace (u0, v0) = go' u0 v0 []
      where
        maxSteps = 3 * nv  -- safety bound
        go' u v acc
            | length acc > maxSteps = acc  -- prevent infinite loop on bad CW data
            | not (null acc) && u == u0 && v == v0 = acc  -- closed the face
            | otherwise = go' v (nextCW g v u) (acc ++ [u])

    dartsOfFace vs =
        [(vs !! i, vs !! ((i + 1) `mod` length vs)) | i <- [0..length vs - 1]]

-- | Check that the CW rotation system defines a valid spherical triangulation.
-- Returns Right () if valid, Left errorMsg otherwise.
checkTriangulation :: DualGraph -> Either String ()
checkTriangulation g =
    let faces = computeFaces g
        nv = numVertices g
        nf = length faces
        -- Count directed edges / 2
        ne = sum (map (deg g) [0..nv-1]) `div` 2
        euler = nv - ne + nf
        nonTriFaces = filter (\f -> length f /= 3) faces
    in if euler /= 2
       then Left $ "Euler characteristic = " ++ show euler ++ " (expected 2)"
                 ++ ", V=" ++ show nv ++ " E=" ++ show ne ++ " F=" ++ show nf
       else if not (null nonTriFaces)
            then Left $ show (length nonTriFaces) ++ " non-triangular faces, sizes: "
                      ++ show (map length nonTriFaces)
            else Right ()

-- | Test helper: check that graph has valid CW triangulation
checkCWValid :: String -> DualGraph -> TestResult
checkCWValid name g = case checkTriangulation g of
    Right () -> Pass (name ++ " CW valid")
    Left err -> Fail (name ++ " CW valid") err

---------------------------------------------------------------------
-- Validation
---------------------------------------------------------------------

-- | Check that a DualGraph passes all structural invariants
checkValid :: String -> DualGraph -> TestResult
checkValid name g = case validateDualGraph g of
    Right () -> Pass (name ++ " validates")
    Left err -> Fail (name ++ " validates") err

-- | Check round-trip: expand then reduce should give back the original graph
checkRoundTrip :: String -> DualGraph -> Expansion -> TestResult
checkRoundTrip name g exp =
    let (g', pi) = applyExpansion exp g
        g''      = applyReduction exp pi g'
    in if g'' == g
       then Pass (name ++ " round-trip")
       else Fail (name ++ " round-trip")
            ("Graphs differ after round-trip.\n  Original nv=" ++ show (numVertices g)
            ++ "\n  Expanded nv=" ++ show (numVertices g')
            ++ "\n  Reduced nv=" ++ show (numVertices g''))

-- | Check that expansion produces a valid graph with correct vertex count
checkExpansionValid :: String -> DualGraph -> Expansion -> TestResult
checkExpansionValid name g exp =
    let (g', _) = applyExpansion exp g
    in case validateDualGraph g' of
        Right () -> Pass (name ++ " expansion valid")
        Left err -> Fail (name ++ " expansion valid") err

-- | Check expanded graph has correct vertex count
checkExpansionSize :: String -> DualGraph -> Expansion -> Int -> TestResult
checkExpansionSize name g exp expectedNv =
    let (g', _) = applyExpansion exp g
        actualNv = numVertices g'
    in if actualNv == expectedNv
       then Pass (name ++ " size=" ++ show expectedNv)
       else Fail (name ++ " size=" ++ show expectedNv)
            ("Expected " ++ show expectedNv ++ " vertices, got " ++ show actualNv)

---------------------------------------------------------------------
-- Tests
---------------------------------------------------------------------

-- | Test graph primitives on C20
testPrimitives :: [TestResult]
testPrimitives =
    [ -- C20 has 12 vertices, all degree-5
      check "c20 has 12 vertices" (numVertices c20 == 12)
    , check "c20 all degree-5" (all (\v -> deg c20 v == 5) [0..11])
    , check "c20 validates" (validateDualGraph c20 == Right ())
    -- nextCW/prevCW are inverses
    , check "nextCW . prevCW = id" $
        all (\u -> all (\v -> prevCW c20 u (nextCW c20 u v) == v) (nbrs c20 u)) [0..11]
    , check "prevCW . nextCW = id" $
        all (\u -> all (\v -> nextCW c20 u (prevCW c20 u v) == v) (nbrs c20 u)) [0..11]
    -- advanceCW 0 = identity
    , check "advanceCW 0 = id" $
        all (\u -> all (\v -> advanceCW c20 u v 0 == v) (nbrs c20 u)) [0..11]
    -- advanceCW deg = identity (full cycle)
    , check "advanceCW deg = id" $
        all (\u -> all (\v -> advanceCW c20 u v (deg c20 u) == v) (nbrs c20 u)) [0..11]
    ]
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test L0 expansions on C20
testL0OnC20 :: [TestResult]
testL0OnC20 =
    let exps = expansionsL0 c20
    in [ check ("c20 has " ++ show (length exps) ++ " L0 expansions")
              (length exps > 0)
       -- C20: all 12 vertices are degree-5. Each vertex has 5 neighbors,
       -- all of which are degree-5. Number of edges between degree-5 pairs:
       -- 30 edges (each vertex has 5 neighbors, 12*5/2 = 30). Each with 2 directions = 60.
       , check "c20 has 120 L0 expansions (bidirectional)" (length exps == 120)
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_L0_" ++ show i
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp 14  -- 12 + 2 = 14
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 exps))  -- test first 10 for speed
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test that all L0 round-trips work on C20
testAllL0RoundTrips :: [TestResult]
testAllL0RoundTrips =
    let exps = expansionsL0 c20
    in map (\(i, exp) -> checkRoundTrip ("c20_L0_all_" ++ show i) c20 exp)
           (zip [0..] exps)

-- | Test expansions on C28
testC28 :: [TestResult]
testC28 =
    let l0s = expansionsL0 c28
    in [ check "c28 validates" (validateDualGraph c28 == Right ())
       , check "c28 has 16 vertices" (numVertices c28 == 16)
       , check ("c28 has " ++ show (length l0s) ++ " L0 expansions") True
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c28_L0_" ++ show i
        in [ checkExpansionValid name c28 exp
           , checkExpansionSize name c28 exp 18
           , checkRoundTrip name c28 exp
           ]) (zip [0..] (take 5 l0s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test straight expansions (L_i, i >= 1) on expanded C20
testStraightOnExpandedC20 :: [TestResult]
testStraightOnExpandedC20 =
    -- First expand C20 via L0 to get a 14-vertex graph, then look for L1 expansions
    let exp0 = head (expansionsL0 c20)
        (g14, _) = applyExpansion exp0 c20
    in [ checkValid "c20+L0" g14
       , check "c20+L0 has 14 vertices" (numVertices g14 == 14)
       ]
    ++ let l1s = expansionsL 4 g14
       in [ check ("c20+L0 has " ++ show (length l1s) ++ " L_i expansions") True ]
       ++ concatMap (\(i, exp) ->
            let name = "c20+L0_L" ++ show (expI exp) ++ "_" ++ show i
            in [ checkExpansionValid name g14 exp
               , checkRoundTrip name g14 exp
               ]) (zip [0..] (take 5 l1s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"
    expI (Exp (L i) _ _) = i
    expI _ = -1

-- | Test B_{0,0} on C20
testBentZeroOnC20 :: [TestResult]
testBentZeroOnC20 =
    let b00s = [ Exp (B 0 0) (u, v) d
               | u <- degree5 c20
               , v <- nbrs c20 u
               , d <- [DLeft, DRight]
               , canBentPath c20 (u, v) d 0 0
               ]
    in [ check ("c20 has " ++ show (length b00s) ++ " B00 candidates") True ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_B00_" ++ show i
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp 15  -- 12 + 3 = 15
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 b00s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test general bent expansions on C20
testBentOnC20 :: [TestResult]
testBentOnC20 =
    let bents = expansionsB 4 c20
        -- Filter to just B_{i,j} with i+j > 0
        nonZeroBents = filter (\(Exp (B i j) _ _) -> i + j > 0) bents
    in [ check ("c20 has " ++ show (length nonZeroBents) ++ " B_{i,j} (i+j>0) expansions") True ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_B_" ++ show i
            Exp (B bi bj) _ _ = exp
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp (12 + bi + bj + 3)
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 nonZeroBents))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test on C30
testC30 :: [TestResult]
testC30 =
    let l0s = expansionsL0 c30
    in [ check "c30 validates" (validateDualGraph c30 == Right ())
       , check "c30 has 17 vertices" (numVertices c30 == 17)
       , check ("c30 has " ++ show (length l0s) ++ " L0 expansions") True
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c30_L0_" ++ show i
        in [ checkExpansionValid name c30 exp
           , checkRoundTrip name c30 exp
           ]) (zip [0..] (take 5 l0s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test bent expansions on a larger graph (expand C28 to get ~20 vertices)
testBentOnExpanded :: [TestResult]
testBentOnExpanded =
    -- Expand C28 via a couple of L0s to get a larger graph
    let exp0 = head (expansionsL0 c28)
        (g18, _) = applyExpansion exp0 c28
        exp1 = head (expansionsL0 g18)
        (g20, _) = applyExpansion exp1 g18
    in [ checkValid "c28+2L0" g20
       , check ("c28+2L0 has " ++ show (numVertices g20) ++ " vertices") (numVertices g20 == 20)
       ]
    ++ let bents = expansionsB 4 g20
       in [ check ("c28+2L0 has " ++ show (length bents) ++ " bent expansions") True ]
       ++ concatMap (\(i, exp) ->
            let name = "c28+2L0_B_" ++ show i
                Exp (B bi bj) _ _ = exp
            in [ checkExpansionValid name g20 exp
               , checkExpansionSize name g20 exp (numVertices g20 + bi + bj + 3)
               , checkRoundTrip name g20 exp
               ]) (zip [0..] (take 10 bents))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test all expansion types comprehensively on an expanded graph
testAllExpansionTypes :: [TestResult]
testAllExpansionTypes =
    let -- Build a larger graph by expanding C20 multiple times
        exp0 = head (expansionsL0 c20)
        (g14, _) = applyExpansion exp0 c20
        exp1 = head (expansionsL0 g14)
        (g16, _) = applyExpansion exp1 g14
        exp2 = head (expansionsL0 g16)
        (g18, _) = applyExpansion exp2 g16
        -- Now test ALL expansion types on g18
        allExps = expansions 4 g18
        l0s   = [e | e@(Exp (L 0) _ _) <- allExps]
        lis   = [e | e@(Exp (L i) _ _) <- allExps, i > 0]
        bents = [e | e@(Exp (B _ _) _ _) <- allExps]
    in [ checkValid "c20+3L0" g18
       , check ("c20+3L0 has " ++ show (numVertices g18) ++ " vertices") (numVertices g18 == 18)
       , check ("All expansions: " ++ show (length allExps)
                ++ " (L0=" ++ show (length l0s)
                ++ " L_i=" ++ show (length lis)
                ++ " B=" ++ show (length bents) ++ ")") True
       ]
    ++ concatMap (\(i, exp) ->
        let Exp k _ _ = exp
            name = "g18_" ++ show k ++ "_" ++ show i
        in [ checkExpansionValid name g18 exp
           , checkRoundTrip name g18 exp
           ]) (zip [0..] (take 20 allExps))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test F (nanotube ring) expansion/reduction round-trip on C30
testRingExpansion :: [TestResult]
testRingExpansion =
    case findNanotubeRing c30 of
        Nothing -> [Fail "c30 ring detection" "findNanotubeRing returned Nothing"]
        Just (ring, outer) ->
            let g22 = applyRing c30 ring outer
                g17 = reduceRing g22 ring outer
            in [ check "c30 has nanotube ring" True
               , check "c30 ring has 5 vertices" (length ring == 5)
               , check "c30 outer has 5 vertices" (length outer == 5)
               , check "c30+F has 22 vertices" (numVertices g22 == 22)
               , checkValid "c30+F" g22
               , checkCWValid "c30+F" g22
               , check "c30+F round-trip" (g17 == c30)
               -- Double expansion: apply F twice
               ] ++ case findNanotubeRing g22 of
                    Nothing -> [Fail "c30+F ring detection" "findNanotubeRing returned Nothing on expanded graph"]
                    Just (ring2, outer2) ->
                        let g27 = applyRing g22 ring2 outer2
                            g22' = reduceRing g27 ring2 outer2
                        in [ check "c30+2F has 27 vertices" (numVertices g27 == 27)
                           , checkValid "c30+2F" g27
                           , checkCWValid "c30+2F" g27
                           , check "c30+2F round-trip" (g22' == g22)
                           ]
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test L/B expansions applied to F-expanded graphs.
-- This catches CW ordering bugs in applyRing that round-trip tests miss
-- (same blind spot as Bug #7).
testExpansionsOnFExpanded :: [TestResult]
testExpansionsOnFExpanded =
    case findNanotubeRing c30 of
        Nothing -> [Fail "F+L/B: ring detection" "findNanotubeRing returned Nothing"]
        Just (ring, outer) ->
            let g22 = applyRing c30 ring outer
                -- L0 expansions on the F-expanded graph
                l0s = expansionsL0 g22
                -- L_i expansions (i >= 1)
                lis = expansionsL 4 g22
                -- Bent expansions
                bents = expansionsB 4 g22
            in [ check ("c30+F has " ++ show (length l0s) ++ " L0 expansions") (length l0s >= 0)
               , check ("c30+F has " ++ show (length lis) ++ " L_i expansions") True
               , check ("c30+F has " ++ show (length bents) ++ " bent expansions") True
               ]
            -- Test L0 expansions: validate, CW check, round-trip
            ++ concatMap (\(i, ex) ->
                let name = "c30+F_L0_" ++ show i
                    (g', pi) = applyExpansion ex g22
                    g'' = applyReduction ex pi g'
                in [ checkValid name g'
                   , checkCWValid name g'
                   , if g'' == g22
                     then Pass (name ++ " round-trip")
                     else Fail (name ++ " round-trip") "reduced /= original"
                   ]) (zip [0..] (take 10 l0s))
            -- Test L_i expansions
            ++ concatMap (\(i, ex) ->
                let Exp (L li) _ _ = ex
                    name = "c30+F_L" ++ show li ++ "_" ++ show i
                    (g', pi) = applyExpansion ex g22
                    g'' = applyReduction ex pi g'
                in [ checkValid name g'
                   , checkCWValid name g'
                   , if g'' == g22
                     then Pass (name ++ " round-trip")
                     else Fail (name ++ " round-trip") "reduced /= original"
                   ]) (zip [0..] (take 10 lis))
            -- Test bent expansions
            ++ concatMap (\(i, ex) ->
                let Exp (B bi bj) _ _ = ex
                    name = "c30+F_B" ++ show bi ++ show bj ++ "_" ++ show i
                    (g', pi) = applyExpansion ex g22
                    g'' = applyReduction ex pi g'
                in [ checkValid name g'
                   , checkCWValid name g'
                   , if g'' == g22
                     then Pass (name ++ " round-trip")
                     else Fail (name ++ " round-trip") "reduced /= original"
                   ]) (zip [0..] (take 10 bents))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test L/B expansions on double-F-expanded graph (c30 + 2F = 27 vertices)
testExpansionsOnDoubleFExpanded :: [TestResult]
testExpansionsOnDoubleFExpanded =
    case findNanotubeRing c30 of
        Nothing -> [Fail "2F+L/B: ring detection" "findNanotubeRing returned Nothing"]
        Just (ring, outer) ->
            let g22 = applyRing c30 ring outer
            in case findNanotubeRing g22 of
                Nothing -> [Fail "2F+L/B: second ring detection" "findNanotubeRing returned Nothing"]
                Just (ring2, outer2) ->
                    let g27 = applyRing g22 ring2 outer2
                        allExps = expansions 4 g27
                    in [ check ("c30+2F has " ++ show (length allExps) ++ " expansions") (length allExps > 0) ]
                    ++ concatMap (\(i, ex) ->
                        let Exp k _ _ = ex
                            name = "c30+2F_" ++ show k ++ "_" ++ show i
                            (g', pi) = applyExpansion ex g27
                            g'' = applyReduction ex pi g'
                        in [ checkValid name g'
                           , checkCWValid name g'
                           , if g'' == g27
                             then Pass (name ++ " round-trip")
                             else Fail (name ++ " round-trip") "reduced /= original"
                           ]) (zip [0..] (take 15 allExps))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

---------------------------------------------------------------------
-- Main
---------------------------------------------------------------------

main :: IO ()
main = do
    putStrLn "=== Graph Primitives ==="
    runTests testPrimitives

    putStrLn "\n=== L0 Expansions on C20 ==="
    runTests testL0OnC20

    putStrLn "\n=== All L0 Round-Trips on C20 ==="
    runTests testAllL0RoundTrips

    putStrLn "\n=== C28 Tests ==="
    runTests testC28

    putStrLn "\n=== Straight Expansions on Expanded C20 ==="
    runTests testStraightOnExpandedC20

    putStrLn "\n=== B_{0,0} on C20 ==="
    runTests testBentZeroOnC20

    putStrLn "\n=== General Bent on C20 ==="
    runTests testBentOnC20

    putStrLn "\n=== C30 Tests ==="
    runTests testC30

    putStrLn "\n=== Bent Expansions on Expanded C28 ==="
    runTests testBentOnExpanded

    putStrLn "\n=== All Expansion Types on Expanded C20 ==="
    runTests testAllExpansionTypes

    putStrLn "\n=== F (Nanotube Ring) Expansion ==="
    runTests testRingExpansion

    putStrLn "\n=== L/B Expansions on F-Expanded C30 ==="
    runTests testExpansionsOnFExpanded

    putStrLn "\n=== L/B Expansions on Double-F-Expanded C30 ==="
    runTests testExpansionsOnDoubleFExpanded
