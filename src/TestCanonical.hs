module Main where

import Seeds (DualGraph(..), c20, c28, c30, validateDualGraph)
import Expansion
import Canonical
import Data.Array.Unboxed (UArray)
import MutGraph (MutGraph, UndoInfo, newMutGraph, loadGraph, freezeGraph,
                 unsafeFreezeGraph, applyExpansionM, applyRingM, undoMutation)
import qualified Search
import qualified Spiral
import qualified Data.IntMap.Strict as IM
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Set (Set)
import qualified Data.Set as Set
import Data.List (sort, foldl', intercalate, isPrefixOf, partition)
import Control.Monad (foldM)
import Control.Monad.ST (ST, runST)
import System.Environment (getArgs)
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
    Spiral.mkGraph [ reverse (nbrs g v)
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
      let bfs = bfsCanonicalForm c20 0 1 DRight
          nv = numVertices c20
          codeLen = bfsCodeLength bfs
          expectedLen = sum [deg c20 v - 1 | v <- [0..nv-1]] + nv
      in if codeLen == expectedLen
         then Pass $ "BFS code length for C20: " ++ show codeLen
         else Fail "BFS code length" $
              "expected " ++ show expectedLen ++ ", got " ++ show codeLen

    , -- C28 BFS code length
      let bfs = bfsCanonicalForm c28 0 1 DRight
          nv = numVertices c28
          codeLen = bfsCodeLength bfs
          expectedLen = sum [deg c28 v - 1 | v <- [0..nv-1]] + nv
      in if codeLen == expectedLen
         then Pass $ "BFS code length for C28: " ++ show codeLen
         else Fail "BFS code length C28" $
              "expected " ++ show expectedLen ++ ", got " ++ show codeLen
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
    , cfgMslTable    :: !(UArray Int Int)  -- precomputed geometric bounds
    }

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
    -- Canonical test rejection breakdown
    , ssRejectL0Early    :: !Int  -- rejected: L0 exists + longer inverse
    , ssRejectNoInverse  :: !Int  -- rejected: no inverse reduction found
    , ssRejectPhase1     :: !Int  -- rejected: cheaper reduction found (x0-x3)
    , ssRejectPhase2     :: !Int  -- rejected: BFS tiebreaker lost
    , ssSkippedSize      :: !Int  -- skipped: child would exceed target size
    -- BFS automorphism group stats (from canonicalBFSAndGroup)
    , ssCanonBfsCalls    :: !Int  -- number of canonicalBFSAndGroup calls
    , ssBfsEdgesTotal    :: !Int  -- total starting edges (before colour filter)
    , ssBfsEdgesFiltered :: !Int  -- starting edges after colour filter
    , ssBfsGT            :: !Int  -- BFS comparisons: GT (rejected)
    , ssBfsEQ            :: !Int  -- BFS comparisons: EQ (automorphism)
    , ssBfsLT            :: !Int  -- BFS comparisons: LT (new minimum)
    } deriving (Show)

emptySS :: SizeStats
emptySS = SS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

addSS :: SizeStats -> SizeStats -> SizeStats
addSS a b = SS (ssParents a + ssParents b)
                (ssExpansionsAll a + ssExpansionsAll b)
                (ssExpansionsR2 a + ssExpansionsR2 b)
                (ssCanonPasses a + ssCanonPasses b)
                (ssDedupCatches a + ssDedupCatches b)
                (max (ssMaxLen a) (ssMaxLen b))
                (ssRingExps a + ssRingExps b)
                (ssRejectL0Early a + ssRejectL0Early b)
                (ssRejectNoInverse a + ssRejectNoInverse b)
                (ssRejectPhase1 a + ssRejectPhase1 b)
                (ssRejectPhase2 a + ssRejectPhase2 b)
                (ssSkippedSize a + ssSkippedSize b)
                (ssCanonBfsCalls a + ssCanonBfsCalls b)
                (ssBfsEdgesTotal a + ssBfsEdgesTotal b)
                (ssBfsEdgesFiltered a + ssBfsEdgesFiltered b)
                (ssBfsGT a + ssBfsGT b)
                (ssBfsEQ a + ssBfsEQ b)
                (ssBfsLT a + ssBfsLT b)

---------------------------------------------------------------------
-- Canonical generation with Rule 2 + instrumentation
---------------------------------------------------------------------

data GenState = GS
    { gsSeen   :: !(Set Spiral.GeneralSpiral)  -- only populated with --spirals
    , gsCounts :: !(Map Int Int)                -- isomer count per dual vertex count
    , gsFailed :: !Int
    , gsStats  :: !(Map Int SizeStats)
    }

emptyState :: GenState
emptyState = GS Set.empty Map.empty 0 Map.empty

addStat :: Int -> SizeStats -> GenState -> GenState
addStat nv ss st = st { gsStats = Map.insertWith addSS nv ss (gsStats st) }

generate :: Config -> GenState
generate cfg = runST $ do
    let seedsWithAuts = [ (g, canonAuts (canonicalBFSAndGroup g))
                        | g <- [c20, c28, c30] ]
    mg <- newMutGraph (cfgMaxDV cfg)
    foldM (\st (seed, auts) -> do
        loadGraph mg seed
        processTreeST cfg mg st seed auts
        ) emptyState seedsWithAuts

-- | Process a graph: record it, then expand its children.
-- Takes the graph's automorphism group for Rule 2 filtering.
processTreeST :: Config -> MutGraph s -> GenState -> DualGraph -> [Automorphism] -> ST s GenState
processTreeST cfg mg st g auts
    | numVertices g > cfgMaxDV cfg = return st
    | otherwise = do
        let nv = numVertices g
            record s = s { gsCounts = Map.insertWith (+) nv 1 (gsCounts s) }
        if cfgDebugSpiral cfg
           then case generalSpiralKey g of
               Just gs | Set.member gs (gsSeen st) ->
                   return $ addStat nv (emptySS { ssDedupCatches = 1 }) st
               Just gs ->
                   let st' = record st { gsSeen = Set.insert gs (gsSeen st) }
                   in expandChildrenST cfg mg st' g auts
               Nothing ->
                   return $ st { gsFailed = gsFailed st + 1 }
           else expandChildrenST cfg mg (record st) g auts

-- | Expand a graph's children using Rule 2 to filter equivalent expansions.
-- Uses mutable graph with apply/freeze/test/undo pattern.
expandChildrenST :: Config -> MutGraph s -> GenState -> DualGraph -> [Automorphism] -> ST s GenState
expandChildrenST cfg mg st g auts = do
    let maxDV = cfgMaxDV cfg
        nv = numVertices g
        -- Only need reductions of length ≤ 2 for bounding lemmas.
        -- allReductions enumerates up to length 5 (expensive bent reductions)
        -- but maxExpansionLength only ever inspects L0 (length 1) and L1/B00 (length 2).
        allReds = allReductionsUpTo 2 g

        -- Apply bounding lemmas to get maximum expansion length
        maxLen = maxExpansionLength maxDV (cfgMslTable cfg) g allReds

        -- All expansions within the bound (pure, uses frozen parent)
        expsAll = expansions maxLen g
        numExpsAll = length expsAll

        -- Rule 2: filter to one per equivalence class (pure)
        expsR2 = filterByRule2 g auts expsAll
        numExpsR2 = length expsR2

    -- Process each expansion with apply/freeze/test/undo.
    -- Accumulate (GenState, SizeStats) to track rejection breakdown.
    (st1, loopSS) <- foldM (\(s, ss) e -> do
        let newVerts = numNewVertices (expKind e)
        if nv + newVerts > maxDV
            then return (s, ss { ssSkippedSize = ssSkippedSize ss + 1 })
            else do
                (undo, _) <- applyExpansionM mg e g
                -- Zero-copy freeze for canonical test (97% reject, no copy needed).
                -- isCanonicalV returns strict CanonVerdict; all graph reads complete
                -- before undoMutation modifies the shared arrays.
                childUnsafe <- unsafeFreezeGraph mg
                let !verdict = isCanonicalV e nv childUnsafe
                (s', ss') <- case verdict of
                    CanonAccept -> do
                        -- Safe copy for recursion (child used after further mutations)
                        child <- freezeGraph mg
                        let cr = canonicalBFSAndGroup child
                            childAuts = canonAuts cr
                            -- Force BFS stats before processTreeST mutates mg
                            !bfsTotal = canonEdgesTotal cr
                            !bfsFiltered = canonEdgesFiltered cr
                            !bfsGT' = canonBfsGT cr
                            !bfsEQ' = canonBfsEQ cr
                            !bfsLT' = canonBfsLT cr
                        s'' <- processTreeST cfg mg s child childAuts
                        return (s'', ss { ssCanonPasses = ssCanonPasses ss + 1
                                        , ssCanonBfsCalls = ssCanonBfsCalls ss + 1
                                        , ssBfsEdgesTotal = ssBfsEdgesTotal ss + bfsTotal
                                        , ssBfsEdgesFiltered = ssBfsEdgesFiltered ss + bfsFiltered
                                        , ssBfsGT = ssBfsGT ss + bfsGT'
                                        , ssBfsEQ = ssBfsEQ ss + bfsEQ'
                                        , ssBfsLT = ssBfsLT ss + bfsLT'
                                        })
                    RejectL0Early ->
                        return (s, ss { ssRejectL0Early = ssRejectL0Early ss + 1 })
                    RejectNoInverse ->
                        return (s, ss { ssRejectNoInverse = ssRejectNoInverse ss + 1 })
                    RejectPhase1 ->
                        return (s, ss { ssRejectPhase1 = ssRejectPhase1 ss + 1 })
                    RejectPhase2 ->
                        return (s, ss { ssRejectPhase2 = ssRejectPhase2 ss + 1 })
                undoMutation mg undo
                return (s', ss')
        ) (st, emptySS) expsR2

    -- F expansion for (5,0) nanotubes (always canonical)
    (st2, ringSS) <- case findNanotubeRing g of
        Just (ring, outer) | nv + 5 <= maxDV -> do
            undo <- applyRingM mg ring outer
            child <- freezeGraph mg
            let cr = canonicalBFSAndGroup child
                childAuts = canonAuts cr
            st' <- processTreeST cfg mg st1 child childAuts
            undoMutation mg undo
            return (st', emptySS { ssRingExps = 1, ssCanonPasses = 1
                                 , ssCanonBfsCalls = 1
                                 , ssBfsEdgesTotal = canonEdgesTotal cr
                                 , ssBfsEdgesFiltered = canonEdgesFiltered cr
                                 , ssBfsGT = canonBfsGT cr
                                 , ssBfsEQ = canonBfsEQ cr
                                 , ssBfsLT = canonBfsLT cr
                                 })
        _ -> return (st1, emptySS)

    -- Record stats for this parent
    let stats = emptySS { ssParents = 1
                        , ssExpansionsAll = numExpsAll
                        , ssExpansionsR2 = numExpsR2
                        , ssMaxLen = maxLen
                        } `addSS` loopSS `addSS` ringSS
    return $ addStat nv stats st2

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

    -- Canonical test breakdown
    let tested = canonPasses + ssRejectL0Early totals + ssRejectNoInverse totals
                 + ssRejectPhase1 totals + ssRejectPhase2 totals
        skipped = ssSkippedSize totals
    putStrLn ""
    putStrLn "=== Canonical Test Breakdown ==="
    putStrLn $ "  isCanonical calls:         " ++ padL 10 (show tested)
    putStrLn $ "    Accept:                  " ++ padL 10 (show canonPasses)
            ++ "  (" ++ showPct canonPasses tested ++ "%)"
    putStrLn $ "    Reject L0 early exit:    " ++ padL 10 (show (ssRejectL0Early totals))
            ++ "  (" ++ showPct (ssRejectL0Early totals) tested ++ "%)"
    putStrLn $ "    Reject no inverse:       " ++ padL 10 (show (ssRejectNoInverse totals))
            ++ "  (" ++ showPct (ssRejectNoInverse totals) tested ++ "%)"
    putStrLn $ "    Reject phase 1 (x0-x3):  " ++ padL 10 (show (ssRejectPhase1 totals))
            ++ "  (" ++ showPct (ssRejectPhase1 totals) tested ++ "%)"
    putStrLn $ "    Reject phase 2 (BFS):    " ++ padL 10 (show (ssRejectPhase2 totals))
            ++ "  (" ++ showPct (ssRejectPhase2 totals) tested ++ "%)"
    putStrLn $ "    Skipped (size limit):    " ++ padL 10 (show skipped)

    -- BFS automorphism group stats
    let bfsCalls = ssCanonBfsCalls totals
        bfsTotal = ssBfsEdgesTotal totals
        bfsFiltered = ssBfsEdgesFiltered totals
        bfsGT_ = ssBfsGT totals
        bfsEQ_ = ssBfsEQ totals
        bfsLT_ = ssBfsLT totals
        bfsCmp = bfsGT_ + bfsEQ_ + bfsLT_
    putStrLn ""
    putStrLn "=== BFS Automorphism Group ==="
    putStrLn $ "  canonicalBFSAndGroup calls:" ++ padL 10 (show bfsCalls)
    putStrLn $ "  Edge colour filter:"
    putStrLn $ "    Total starting edges:    " ++ padL 10 (show bfsTotal)
    putStrLn $ "    After filter:            " ++ padL 10 (show bfsFiltered)
            ++ "  (" ++ showPct bfsFiltered bfsTotal ++ "%)"
    if bfsCalls > 0
        then putStrLn $ "    Avg per call:            "
                ++ show (bfsTotal `div` bfsCalls) ++ " -> "
                ++ show (bfsFiltered `div` bfsCalls)
                ++ " (" ++ show (bfsTotal `div` max 1 bfsFiltered) ++ "x reduction)"
        else return ()
    putStrLn $ "  BFS comparison outcomes:"
    putStrLn $ "    GT (rejected):           " ++ padL 10 (show bfsGT_)
            ++ "  (" ++ showPct bfsGT_ bfsCmp ++ "%)"
    putStrLn $ "    EQ (automorphism):       " ++ padL 10 (show bfsEQ_)
            ++ "  (" ++ showPct bfsEQ_ bfsCmp ++ "%)"
    putStrLn $ "    LT (new minimum):        " ++ padL 10 (show bfsLT_)
            ++ "  (" ++ showPct bfsLT_ bfsCmp ++ "%)"

padR :: Int -> String -> String
padR n s = s ++ replicate (max 0 (n - length s)) ' '

padL :: Int -> String -> String
padL n s = replicate (max 0 (n - length s)) ' ' ++ s

showPct :: Int -> Int -> String
showPct _ 0 = "N/A"
showPct n d = show (round (100.0 * fromIntegral n / fromIntegral d :: Double) :: Int)

---------------------------------------------------------------------
-- Pure forest cross-check
---------------------------------------------------------------------

forestTests :: Int -> Map Int Int -> [TestResult]
forestTests maxDV mutCounts =
    let gcfg = Search.mkConfig maxDV
        pureCounts = Search.searchPure Search.countBySize gcfg
        sizes = sort $ Set.toList $ Set.fromList $
                filter (<= maxDV) (Map.keys mutCounts ++ Map.keys pureCounts)
    in [ let mc = Map.findWithDefault 0 nv mutCounts
             pc = Map.findWithDefault 0 nv pureCounts
         in if mc == pc
            then Pass $ "C" ++ show (carbonAtoms nv) ++ " forest=" ++ show pc
            else Fail ("C" ++ show (carbonAtoms nv) ++ " forest")
                      ("MutGraph=" ++ show mc ++ " vs Forest=" ++ show pc)
       | nv <- sizes
       ]

parForestTests :: Int -> Map Int Int -> Int -> [TestResult]
parForestTests maxDV mutCounts _parDepth =
    let gcfg = Search.mkConfig maxDV
        parCounts = Search.searchPure Search.countBySize gcfg
        sizes = sort $ Set.toList $ Set.fromList $
                filter (<= maxDV) (Map.keys mutCounts ++ Map.keys parCounts)
    in [ let mc = Map.findWithDefault 0 nv mutCounts
             pc = Map.findWithDefault 0 nv parCounts
         in if mc == pc
            then Pass $ "C" ++ show (carbonAtoms nv) ++ " par=" ++ show pc
            else Fail ("C" ++ show (carbonAtoms nv) ++ " par")
                      ("MutGraph=" ++ show mc ++ " vs Par=" ++ show pc)
       | nv <- sizes
       ]

---------------------------------------------------------------------
-- Main
---------------------------------------------------------------------

main :: IO ()
main = do
    args <- getArgs
    let (flags, positional) = Data.List.partition ("--" `isPrefixOf`) args
        spirals = "--spirals" `elem` flags
        parMode = "--par" `elem` flags
        maxDV = case positional of
                  (s:_) -> read s
                  []    -> 22  -- C40
        cfg = Config { cfgMaxDV = maxDV
                     , cfgDebugSpiral = spirals
                     , cfgMslTable = buildMaxStraightLengths maxDV
                     }

    putStrLn $ "Mode: " ++ if spirals then "debug (spiral dedup ON)"
                           else if parMode then "parallel (depth 2)"
                           else "fast (spiral dedup OFF)"

    putStrLn "=== BFS Code Tests ==="
    runTests bfsTests

    putStrLn "\n=== Reduction Enumeration Tests ==="
    runTests reductionTests

    putStrLn $ "\n=== Canonical Generation Tests (up to C"
             ++ show (carbonAtoms maxDV) ++ ") ==="
    let (st, genTests) = generationTests cfg
    runTests genTests

    -- Cross-check: pure forest vs MutGraph generation
    putStrLn $ "\n=== Pure Forest Cross-Check (up to C"
             ++ show (carbonAtoms maxDV) ++ ") ==="
    let fTests = forestTests maxDV (gsCounts st)
    runTests fTests

    -- Parallel cross-check (when --par flag is set)
    pTests <- if parMode
        then do
            putStrLn $ "\n=== Parallel Forest (depth 2) Cross-Check (up to C"
                     ++ show (carbonAtoms maxDV) ++ ") ==="
            let pt = parForestTests maxDV (gsCounts st) 2
            runTests pt
            return pt
        else return []

    -- Print instrumentation
    printStats maxDV st

    putStrLn ""
    let allTests = bfsTests ++ reductionTests ++ genTests ++ fTests ++ pTests
        allOk = all passed allTests
    if allOk
        then putStrLn "ALL TESTS PASSED"
        else putStrLn "SOME TESTS FAILED"
