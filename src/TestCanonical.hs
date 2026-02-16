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
import Data.List (sort, foldl', intercalate, isPrefixOf, partition)
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
-- Configuration
---------------------------------------------------------------------

data Config = Config
    { cfgMaxDV       :: !Int   -- max dual vertices (= C_n/2 + 2)
    , cfgDebugSpiral :: !Bool  -- compute spiral keys for dedup checking
    } deriving (Show)

---------------------------------------------------------------------
-- Per-size statistics
---------------------------------------------------------------------

data SizeStats = SS
    { ssParents          :: !Int  -- parent graphs expanded at this size
    , ssExpansionsAll    :: !Int  -- total expansion sites (after bounding, before Rule 2)
    , ssExpansionsR2     :: !Int  -- expansion sites after Rule 2 filtering
    , ssCanonPasses      :: !Int  -- passed canonical test
    , ssDedupCatches     :: !Int  -- passed canonical but caught by spiral dedup
    , ssMaxLen           :: !Int  -- max expansion length used (from bounding lemmas)
    , ssRingExps         :: !Int  -- F (nanotube ring) expansions
    } deriving (Show)

emptySS :: SizeStats
emptySS = SS 0 0 0 0 0 0 0

addSS :: SizeStats -> SizeStats -> SizeStats
addSS a b = SS (ssParents a + ssParents b)
                (ssExpansionsAll a + ssExpansionsAll b)
                (ssExpansionsR2 a + ssExpansionsR2 b)
                (ssCanonPasses a + ssCanonPasses b)
                (ssDedupCatches a + ssDedupCatches b)
                (max (ssMaxLen a) (ssMaxLen b))
                (ssRingExps a + ssRingExps b)

---------------------------------------------------------------------
-- Canonical generation with Rule 2 + instrumentation
---------------------------------------------------------------------

data GenState = GS
    { gsSeen   :: !(Set Spiral.GeneralSpiral)  -- only populated with --spirals
    , gsCounts :: !(Map Int Int)                -- isomer count per dual vertex count
    , gsGraphs :: !(Map Int [DualGraph])        -- all generated graphs (always populated)
    , gsFailed :: !Int
    , gsStats  :: !(Map Int SizeStats)
    }

emptyState :: GenState
emptyState = GS Set.empty Map.empty Map.empty 0 Map.empty

addStat :: Int -> SizeStats -> GenState -> GenState
addStat nv ss st = st { gsStats = Map.insertWith addSS nv ss (gsStats st) }

generate :: Config -> GenState
generate cfg =
    -- Compute automorphism groups for seeds
    let seedsWithAuts = [ (g, canonAuts (canonicalBFSAndGroup g))
                        | g <- [c20, c28, c30] ]
    in foldl' (\st (g, auts) -> processTree cfg st g auts) emptyState seedsWithAuts

-- | Process a graph: record it, then expand its children.
-- Takes the graph's automorphism group for Rule 2 filtering.
processTree :: Config -> GenState -> DualGraph -> [Automorphism] -> GenState
processTree cfg st g auts
    | numVertices g > cfgMaxDV cfg = st
    | otherwise =
        let nv = numVertices g
            record s = s { gsCounts = Map.insertWith (+) nv 1 (gsCounts s)
                         , gsGraphs = Map.insertWith (++) nv [g] (gsGraphs s)
                         }
        in if cfgDebugSpiral cfg
           then case generalSpiralKey g of
               Just gs | Set.member gs (gsSeen st) ->
                   addStat nv (emptySS { ssDedupCatches = 1 }) st
               Just gs ->
                   let st' = record st { gsSeen = Set.insert gs (gsSeen st) }
                   in expandChildren cfg st' g auts
               Nothing ->
                   st { gsFailed = gsFailed st + 1 }
           else expandChildren cfg (record st) g auts

-- | Expand a graph's children using Rule 2 to filter equivalent expansions.
expandChildren :: Config -> GenState -> DualGraph -> [Automorphism] -> GenState
expandChildren cfg st g auts =
    let maxDV = cfgMaxDV cfg
        nv = numVertices g
        allReds = allReductions g

        -- Apply bounding lemmas to get maximum expansion length
        maxLen = maxExpansionLength maxDV g allReds

        -- All expansions within the bound
        expsAll = expansions maxLen g
        numExpsAll = length expsAll

        -- Rule 2: filter to one per equivalence class
        expsR2 = filterByRule2 g auts expsAll
        numExpsR2 = length expsR2

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
                        , ssMaxLen = maxLen
                        , ssRingExps = numRing
                        }
        st' = addStat nv stats st

    in foldl' (\s (g', childAuts) -> processTree cfg s g' childAuts)
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

generationTests :: Config -> (GenState, [TestResult])
generationTests cfg =
    let maxDV = cfgMaxDV cfg
        st = generate cfg
        counts = gsCounts st
        relevant = [ (dualVerts c, c, expected)
                   | (c, expected) <- knownCounts
                   , dualVerts c <= maxDV
                   ]
        tests = map (\(dv, cn, expected) ->
            let found = Map.findWithDefault 0 dv counts
            in if found == expected
               then Pass $ "C" ++ show cn ++ ": " ++ show found ++ " isomers"
               else Fail ("C" ++ show cn) $
                    "expected " ++ show expected ++ ", got " ++ show found
            ) relevant
         ++ if cfgDebugSpiral cfg
            then let failedSpirals = gsFailed st
                 in [if failedSpirals == 0
                     then Pass "No spiral failures"
                     else Fail "Spiral failures" $ show failedSpirals ++ " graphs failed"
                    ]
            else []
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
            ++ padR 12 "Exps" ++ padR 10 "Rule2"
            ++ padR 10 "Canon" ++ padR 10 "Dedup"
            ++ padR 8 "MaxL" ++ padR 8 "Ring"
    putStrLn $ replicate 84 '-'

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
               ++ padR 10 (show (ssCanonPasses ss))
               ++ padR 10 (show (ssDedupCatches ss))
               ++ padR 8 (show (ssMaxLen ss))
               ++ padR 8 (show (ssRingExps ss))
        ) relevantSizes

    -- Totals
    let totals = foldl' addSS emptySS (Map.elems stats)
    putStrLn $ replicate 84 '-'
    putStrLn $ padR 8 "TOTAL" ++ padR 8 ""
            ++ padR 10 (show (ssParents totals))
            ++ padR 12 (show (ssExpansionsAll totals))
            ++ padR 10 (show (ssExpansionsR2 totals))
            ++ padR 10 (show (ssCanonPasses totals))
            ++ padR 10 (show (ssDedupCatches totals))
            ++ padR 8 (show (ssMaxLen totals))
            ++ padR 8 (show (ssRingExps totals))

    -- Summary ratios
    let totalExpsAll = ssExpansionsAll totals
        totalExpsR2 = ssExpansionsR2 totals
        canonPasses = ssCanonPasses totals
        dedupCatches = ssDedupCatches totals
        r2Saved = totalExpsAll - totalExpsR2
    putStrLn ""
    putStrLn $ "Expansion sites (bounded): " ++ show totalExpsAll
    putStrLn $ "After Rule 2:              " ++ show totalExpsR2
            ++ " (Rule 2 eliminated " ++ show r2Saved
            ++ ", " ++ showPct r2Saved totalExpsAll ++ "%)"
    putStrLn $ "Canon acceptance rate: " ++ show canonPasses ++ "/"
            ++ show totalExpsR2
            ++ " (" ++ showPct canonPasses totalExpsR2 ++ "%)"
    putStrLn $ "Dedup catches: " ++ show dedupCatches
    if dedupCatches > 0
        then putStrLn $ "WARNING: Dedup caught " ++ show dedupCatches
                     ++ " duplicates — Rule 2 filtering incomplete"
        else putStrLn "OK: zero duplicates"

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
    let (flags, positional) = Data.List.partition ("--" `isPrefixOf`) args
        spirals = "--spirals" `elem` flags
        maxDV = case positional of
                  (s:_) -> read s
                  []    -> 22  -- C40
        cfg = Config { cfgMaxDV = maxDV, cfgDebugSpiral = spirals }

    putStrLn $ "Mode: " ++ if spirals then "debug (spiral dedup ON)" else "fast (spiral dedup OFF)"
    putStrLn "=== BFS Code Tests ==="
    runTests bfsTests

    putStrLn "\n=== Reduction Enumeration Tests ==="
    runTests reductionTests

    putStrLn $ "\n=== Canonical Generation Tests (up to C"
             ++ show (carbonAtoms maxDV) ++ ") ==="
    let (st, genTests) = generationTests cfg
    runTests genTests

    -- Print instrumentation
    printStats maxDV st

    putStrLn ""
    let allTests = bfsTests ++ reductionTests ++ genTests
        allOk = all passed allTests
    if allOk
        then putStrLn "ALL TESTS PASSED"
        else putStrLn "SOME TESTS FAILED"
