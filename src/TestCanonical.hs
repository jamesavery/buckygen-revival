module Main where

import Seeds (DualGraph(..), c20, c28, c30, validateDualGraph)
import Expansion
import Canonical
import qualified Spiral
import qualified Data.IntMap.Strict as IM
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Set (Set)
import qualified Data.Set as Set
import Data.List (sort, foldl', intercalate)
import System.Environment (getArgs)
import Data.IORef
import System.IO (hFlush, stdout)

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
-- DualGraph -> Spiral.Graph conversion (from Generate.hs)
---------------------------------------------------------------------

toSpiralGraph :: DualGraph -> Spiral.Graph
toSpiralGraph g =
    Spiral.mkGraph [ reverse (neighbours g IM.! v)
                   | v <- [0 .. numVertices g - 1] ]

generalSpiralKey :: DualGraph -> Maybe Spiral.GeneralSpiral
generalSpiralKey = Spiral.canonicalGeneralSpiral . toSpiralGraph

---------------------------------------------------------------------
-- BFS code tests
---------------------------------------------------------------------

bfsTests :: [TestResult]
bfsTests =
    [ -- Same edge, same direction => same BFS code
      let code1 = bfsCanonicalForm c20 0 1 DRight
          code2 = bfsCanonicalForm c20 0 1 DRight
      in if code1 == code2
         then Pass "BFS deterministic: same edge same code"
         else Fail "BFS deterministic" "codes differ"

    , -- Different direction => may differ
      let code1 = bfsCanonicalForm c20 0 1 DRight
          code2 = bfsCanonicalForm c20 0 1 DLeft
      in Pass $ "BFS direction: DRight /= DLeft: " ++ show (code1 /= code2)

    , -- BFS code has correct length: sum(degree-1) + nv separators
      let BFS code = bfsCanonicalForm c20 0 1 DRight
          nv = numVertices c20
          expectedLen = sum [deg c20 v - 1 | v <- [0..nv-1]] + nv
      in if length code == expectedLen
         then Pass $ "BFS code length for C20: " ++ show (length code)
         else Fail "BFS code length" $
              "expected " ++ show expectedLen ++ ", got " ++ show (length code)

    , -- C28 BFS code length
      let BFS code = bfsCanonicalForm c28 0 1 DRight
          nv = numVertices c28
          expectedLen = sum [deg c28 v - 1 | v <- [0..nv-1]] + nv
      in if length code == expectedLen
         then Pass $ "BFS code length for C28: " ++ show (length code)
         else Fail "BFS code length C28" $
              "expected " ++ show expectedLen ++ ", got " ++ show (length code)
    ]

---------------------------------------------------------------------
-- Reduction enumeration tests
---------------------------------------------------------------------

reductionTests :: [TestResult]
reductionTests =
    [ -- C20 has NO L0 reductions: all vertices are degree-5, so there are
      -- no flanking degree-6 vertices. This matches the C code's
      -- has_L0_reduction which requires flanking degree-6.
      let reds = allReductions c20
          l0s = [r | r@(Red (L 0) _ _) <- reds]
      in if null l0s
         then Pass "C20 has no L0 reductions (correct: no degree-6 flanking)"
         else Fail "C20 L0 reductions" ("expected none, got " ++ show (length l0s))

    , -- C20 should have reductions (it's a valid fullerene)
      let reds = allReductions c20
      in if not (null reds)
         then Pass $ "C20 total reductions: " ++ show (length reds)
         else Fail "C20 reductions" "none found"

    , -- Expand C20 by L0, check expanded graph has reductions
      let exps = expansionsL0 c20
      in if null exps
         then Fail "C20 L0 expansions" "none found"
         else let e = head exps
                  (g', _) = applyExpansion e c20
                  reds = allReductions g'
              in if not (null reds)
                 then Pass $ "C24 (C20+L0) has reductions: " ++ show (length reds)
                 else Fail "C24 reductions" "none found"
    ]

---------------------------------------------------------------------
-- Per-size statistics
---------------------------------------------------------------------

data SizeStats = SS
    { ssParents          :: !Int  -- parent graphs expanded at this size
    , ssExpansionsAll    :: !Int  -- total expansion sites before Rule 2
    , ssExpansionsR2     :: !Int  -- expansion sites after Rule 2 filtering
    , ssCanonPasses      :: !Int  -- passed canonical test
    , ssDedupCatches     :: !Int  -- passed canonical but caught by spiral dedup
    , ssHasL0            :: !Int  -- parents with an L0 reduction
    , ssExpsIfBounded    :: !Int  -- expansions that would survive Lemma 3 bound
    , ssRingExps         :: !Int  -- F (nanotube ring) expansions
    } deriving (Show)

emptySS :: SizeStats
emptySS = SS 0 0 0 0 0 0 0 0

addSS :: SizeStats -> SizeStats -> SizeStats
addSS a b = SS (ssParents a + ssParents b)
                (ssExpansionsAll a + ssExpansionsAll b)
                (ssExpansionsR2 a + ssExpansionsR2 b)
                (ssCanonPasses a + ssCanonPasses b)
                (ssDedupCatches a + ssDedupCatches b)
                (ssHasL0 a + ssHasL0 b)
                (ssExpsIfBounded a + ssExpsIfBounded b)
                (ssRingExps a + ssRingExps b)

---------------------------------------------------------------------
-- Canonical generation with Rule 2 + instrumentation
---------------------------------------------------------------------

data GenState = GS
    { gsSeen   :: !(Set Spiral.GeneralSpiral)
    , gsGraphs :: !(Map Int [DualGraph])
    , gsFailed :: !Int
    , gsStats  :: !(Map Int SizeStats)
    }

emptyState :: GenState
emptyState = GS Set.empty Map.empty 0 Map.empty

addStat :: Int -> SizeStats -> GenState -> GenState
addStat nv ss st = st { gsStats = Map.insertWith addSS nv ss (gsStats st) }

generate :: Int -> GenState
generate maxDV =
    -- Compute automorphism groups for seeds
    let seedsWithAuts = [ (g, canonAuts (canonicalBFSAndGroup g))
                        | g <- [c20, c28, c30] ]
    in foldl' (\st (g, auts) -> processTree maxDV st g auts) emptyState seedsWithAuts

-- | Process a graph: record it, then expand its children.
-- Takes the graph's automorphism group for Rule 2 filtering.
processTree :: Int -> GenState -> DualGraph -> [Automorphism] -> GenState
processTree maxDV st g auts
    | numVertices g > maxDV = st
    | otherwise =
        let nv = numVertices g
        in case generalSpiralKey g of
            Just gs | Set.member gs (gsSeen st) ->
                -- Dedup catch: should be 0 with correct Rule 2
                addStat nv (emptySS { ssDedupCatches = 1 }) st
            Just gs ->
                let st' = st { gsSeen = Set.insert gs (gsSeen st)
                             , gsGraphs = Map.insertWith (++) nv [g] (gsGraphs st)
                             }
                in expandChildren maxDV st' g auts
            Nothing ->
                st { gsFailed = gsFailed st + 1 }

-- | Expand a graph's children using Rule 2 to filter equivalent expansions.
expandChildren :: Int -> GenState -> DualGraph -> [Automorphism] -> GenState
expandChildren maxDV st g auts =
    let nv = numVertices g
        allReds = allReductions g
        hasL0 = any (\(Red k _ _) -> k == L 0) allReds

        -- Current (unbounded) max length
        maxLen = min (maxDV - nv + 2) 5

        -- What Lemma 3 would give
        boundedMaxLen = if hasL0 then min maxLen 3 else maxLen

        -- All expansions before Rule 2
        expsAll = expansions maxLen g
        numExpsAll = length expsAll

        -- Rule 2: filter to one per equivalence class
        expsR2 = filterByRule2 g auts expsAll
        numExpsR2 = length expsR2

        -- Count how many are within the Lemma 3 bound
        expsInBound = length [ e | e <- expsAll,
                               reductionLength (expKind e) <= boundedMaxLen ]

        -- Apply canonical test to Rule-2-filtered expansions
        -- For each child that passes, also compute its automorphism group
        childrenWithAuts =
            [ (g', childAuts)
            | e <- expsR2
            , let (g', _) = applyExpansion e g
            , numVertices g' <= maxDV
            , isCanonical e nv g'
            , let childAuts = canonAuts (canonicalBFSAndGroup g')
            ]
        numCanonPasses = length childrenWithAuts

        -- F expansion for (5,0) nanotubes (always canonical)
        ringChildrenWithAuts = case findNanotubeRing g of
            Just (ring, outer) | nv + 5 <= maxDV ->
                let g' = applyRing g ring outer
                in [(g', canonAuts (canonicalBFSAndGroup g'))]
            _ -> []
        numRing = length ringChildrenWithAuts

        -- Record stats for this parent
        stats = emptySS { ssParents = 1
                        , ssExpansionsAll = numExpsAll
                        , ssExpansionsR2 = numExpsR2
                        , ssCanonPasses = numCanonPasses
                        , ssHasL0 = if hasL0 then 1 else 0
                        , ssExpsIfBounded = expsInBound
                        , ssRingExps = numRing
                        }
        st' = addStat nv stats st

    in foldl' (\s (g', childAuts) -> processTree maxDV s g' childAuts)
              st' (childrenWithAuts ++ ringChildrenWithAuts)

---------------------------------------------------------------------
-- Known isomer counts
---------------------------------------------------------------------

knownCounts :: [(Int, Int)]
knownCounts =
    [ (20, 1), (24, 1), (26, 1), (28, 2), (30, 3)
    , (32, 6), (34, 6), (36, 15), (38, 17), (40, 40)
    , (42, 45), (44, 89), (46, 116), (48, 199), (50, 271)
    , (52, 437), (54, 580), (56, 924), (58, 1205), (60, 1812)
    ]

dualVerts :: Int -> Int
dualVerts c = c `div` 2 + 2

carbonAtoms :: Int -> Int
carbonAtoms dv = (dv - 2) * 2

---------------------------------------------------------------------
-- Generation count tests
---------------------------------------------------------------------

generationTests :: Int -> (GenState, [TestResult])
generationTests maxDV =
    let st = generate maxDV
        results = gsGraphs st
        relevant = [ (dualVerts c, c, expected)
                   | (c, expected) <- knownCounts
                   , dualVerts c <= maxDV
                   ]
        failedSpirals = gsFailed st
        tests = map (\(dv, cn, expected) ->
            let found = length $ Map.findWithDefault [] dv results
            in if found == expected
               then Pass $ "C" ++ show cn ++ ": " ++ show found ++ " isomers"
               else Fail ("C" ++ show cn) $
                    "expected " ++ show expected ++ ", got " ++ show found
            ) relevant
         ++ [if failedSpirals == 0
             then Pass "No spiral failures"
             else Fail "Spiral failures" $ show failedSpirals ++ " graphs failed"
            ]
    in (st, tests)

---------------------------------------------------------------------
-- Stats reporting
---------------------------------------------------------------------

printStats :: Int -> GenState -> IO ()
printStats maxDV st = do
    let stats = gsStats st
        allSizes = sort $ Map.keys stats

    putStrLn "\n=== Performance Statistics ==="
    putStrLn ""

    -- Header
    putStrLn $ padR 8 "DualV" ++ padR 8 "C_n" ++ padR 10 "Parents"
            ++ padR 12 "ExpsAll" ++ padR 10 "Rule2" ++ padR 12 "Bounded"
            ++ padR 10 "Canon" ++ padR 10 "Dedup" ++ padR 8 "HasL0"
            ++ padR 8 "Ring"
    putStrLn $ replicate 96 '-'

    -- Per-size rows
    let relevantSizes = filter (<= maxDV) allSizes
    mapM_ (\nv -> do
        let ss = Map.findWithDefault emptySS nv stats
            cn = carbonAtoms nv
        putStrLn $ padR 8 (show nv)
               ++ padR 8 ("C" ++ show cn)
               ++ padR 10 (show (ssParents ss))
               ++ padR 12 (show (ssExpansionsAll ss))
               ++ padR 10 (show (ssExpansionsR2 ss))
               ++ padR 12 (show (ssExpsIfBounded ss))
               ++ padR 10 (show (ssCanonPasses ss))
               ++ padR 10 (show (ssDedupCatches ss))
               ++ padR 8 (show (ssHasL0 ss))
               ++ padR 8 (show (ssRingExps ss))
        ) relevantSizes

    -- Totals
    let totals = foldl' addSS emptySS (Map.elems stats)
    putStrLn $ replicate 96 '-'
    putStrLn $ padR 8 "TOTAL" ++ padR 8 ""
            ++ padR 10 (show (ssParents totals))
            ++ padR 12 (show (ssExpansionsAll totals))
            ++ padR 10 (show (ssExpansionsR2 totals))
            ++ padR 12 (show (ssExpsIfBounded totals))
            ++ padR 10 (show (ssCanonPasses totals))
            ++ padR 10 (show (ssDedupCatches totals))
            ++ padR 8 (show (ssHasL0 totals))
            ++ padR 8 (show (ssRingExps totals))

    -- Summary ratios
    let totalExpsAll = ssExpansionsAll totals
        totalExpsR2 = ssExpansionsR2 totals
        boundedExps = ssExpsIfBounded totals
        canonPasses = ssCanonPasses totals
        dedupCatches = ssDedupCatches totals
        r2Saved = totalExpsAll - totalExpsR2
    putStrLn ""
    putStrLn $ "Expansion sites before Rule 2: " ++ show totalExpsAll
    putStrLn $ "Expansion sites after Rule 2:  " ++ show totalExpsR2
            ++ " (Rule 2 eliminated " ++ show r2Saved
            ++ ", " ++ showPct r2Saved totalExpsAll ++ "%)"
    putStrLn $ "Would survive Lemma 3 bound:   " ++ show boundedExps
            ++ " (saves " ++ showPct (totalExpsAll - boundedExps) totalExpsAll
            ++ "% of total)"
    putStrLn $ "Canon acceptance rate: " ++ show canonPasses ++ "/"
            ++ show totalExpsR2
            ++ " (" ++ showPct canonPasses totalExpsR2 ++ "% of Rule-2-filtered)"
    putStrLn $ "Dedup catches (should be 0 with correct Rule 2): "
            ++ show dedupCatches
    if dedupCatches > 0
        then putStrLn $ "WARNING: Dedup caught " ++ show dedupCatches
                     ++ " duplicates — Rule 2 filtering incomplete"
        else putStrLn "CONFIRMED: Rule 2 eliminated all duplicates"

padR :: Int -> String -> String
padR n s = s ++ replicate (max 0 (n - length s)) ' '

showPct :: Int -> Int -> String
showPct _ 0 = "N/A"
showPct n d = show (round (100.0 * fromIntegral n / fromIntegral d :: Double) :: Int)

---------------------------------------------------------------------
-- Main
---------------------------------------------------------------------

main :: IO ()
main = do
    args <- getArgs
    let maxDV = case args of
                  (s:_) -> read s
                  []    -> 22  -- C40

    putStrLn "=== BFS Code Tests ==="
    runTests bfsTests

    putStrLn "\n=== Reduction Enumeration Tests ==="
    runTests reductionTests

    putStrLn $ "\n=== Canonical Generation Tests (up to C"
             ++ show (carbonAtoms maxDV) ++ ") ==="
    let (st, genTests) = generationTests maxDV
    runTests genTests

    -- Print instrumentation
    printStats maxDV st

    putStrLn ""
    let allTests = bfsTests ++ reductionTests ++ genTests
        allOk = all passed allTests
    if allOk
        then putStrLn "ALL TESTS PASSED"
        else putStrLn "SOME TESTS FAILED"
