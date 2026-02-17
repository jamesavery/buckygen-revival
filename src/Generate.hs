{-# LANGUAGE StrictData #-}
-- | Generate the expansion tree from seeds, deduplicate by canonical
-- generalized spiral, and compare counts against known fullerene
-- isomer numbers.
--
-- Cross-checks spiral dedup against backtracking isomorphism for
-- small sizes (up to C40) where backtracking is tractable.
--
-- Usage: runghc Generate.hs [maxDualVerts]
--   default maxDualVerts = 22 (= C40, 40 known isomers)

module Main where

import Seeds (DualGraph(..), initEdgeList, c20, c28, c30, validateDualGraph)
import Expansion hiding (intersect)
import qualified Spiral
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IM
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Set (Set)
import qualified Data.Set as Set
import Data.List (sort, foldl', delete, sortBy)
import Data.Ord (comparing)
import System.Environment (getArgs)

---------------------------------------------------------------------
-- DualGraph -> Spiral.Graph conversion
---------------------------------------------------------------------

-- | Convert a DualGraph (CW adjacency lists in IntMap) to a
-- Spiral.Graph (CCW adjacency lists in Vector). We reverse each
-- neighbor list to flip CW -> CCW.
toSpiralGraph :: DualGraph -> Spiral.Graph
toSpiralGraph g =
    Spiral.mkGraph [ reverse (nbrs g v)
                   | v <- [0 .. numVertices g - 1] ]

-- | Canonical general spiral key for a DualGraph.
-- Uses canonicalGeneralSpiral which handles both regular and
-- generalized spirals (with jumps). Should succeed for all
-- 3-connected planar graphs.
generalSpiralKey :: DualGraph -> Maybe Spiral.GeneralSpiral
generalSpiralKey = Spiral.canonicalGeneralSpiral . toSpiralGraph

---------------------------------------------------------------------
-- Graph isomorphism by backtracking (kept for cross-checking)
---------------------------------------------------------------------

-- | Check if two dual graphs are isomorphic.
isIsomorphic :: DualGraph -> DualGraph -> Bool
isIsomorphic g1 g2
    | numVertices g1 /= numVertices g2 = False
    | deg5Sig g1 /= deg5Sig g2 = False
    | otherwise = tryMapping g1 g2
  where
    deg5Sig g = sort [ sort (map (deg g) (nbrs g v)) | v <- degree5 g ]

tryMapping :: DualGraph -> DualGraph -> Bool
tryMapping g1 g2 = go (sortedVerts g1) IM.empty [] (allVerts g2)
  where
    allVerts g = [0 .. numVertices g - 1]
    vertSig g v = (deg g v, sort (map (deg g) (nbrs g v)))
    sortedVerts g = map snd $ sortBy (comparing fst)
                    [(vertSig g v, v) | v <- allVerts g]

    go [] _mapping _ _ = True
    go (v1:rest1) mapping _used available =
        let d1 = deg g1 v1
            candidates = [ v2 | v2 <- available
                         , deg g2 v2 == d1
                         , compatible v1 v2 mapping ]
        in any (\v2 -> go rest1 (IM.insert v1 v2 mapping) (v2:_used)
                           (delete v2 available))
               candidates

    compatible v1 v2 mapping =
        let n1 = nbrs g1 v1
            n2Set = nbrs g2 v2
            mappedN1 = [ (w1, IM.lookup w1 mapping) | w1 <- n1 ]
        in all (\(_, mw2) -> case mw2 of
                    Nothing -> True
                    Just w2 -> w2 `elem` n2Set)
               mappedN1
        && all (\(w1, w2) -> if w2 `elem` n2Set
                             then w1 `elem` nbrs g1 v1
                             else True)
               [(w1, w2) | (w1, Just w2) <- [ (w, IM.lookup w mapping) | w <- [0..numVertices g1 - 1] ]]

---------------------------------------------------------------------
-- Generation: DFS expansion tree with general-spiral deduplication
---------------------------------------------------------------------

data GenState = GS
    { gsSeen     :: !(Set Spiral.GeneralSpiral)
    , gsGraphs   :: !(Map Int [DualGraph])
    , gsFailed   :: !Int  -- graphs where even general spiral failed
    }

emptyState :: GenState
emptyState = GS Set.empty Map.empty 0

generate :: Int -> GenState
generate maxDV = foldl' (processTree maxDV) emptyState [c20, c28, c30]

processTree :: Int -> GenState -> DualGraph -> GenState
processTree maxDV st g
    | numVertices g > maxDV = st
    | otherwise =
        let nv = numVertices g
        in case generalSpiralKey g of
            Just gs | Set.member gs (gsSeen st) -> st
            Just gs ->
                let st' = st { gsSeen = Set.insert gs (gsSeen st)
                             , gsGraphs = Map.insertWith (++) nv [g] (gsGraphs st)
                             }
                in expandChildren maxDV st' g
            Nothing ->
                -- Should not happen for valid fullerene duals
                let st' = st { gsFailed = gsFailed st + 1 }
                in st'

expandChildren :: Int -> GenState -> DualGraph -> GenState
expandChildren maxDV st g =
    let nv = numVertices g
        maxLen = min 10 (maxDV - nv + 2)
        exps = expansions maxLen g
        children = map (\e -> fst (applyExpansion e g)) exps
        ringChildren = case findNanotubeRing g of
            Just (ring, outer) | nv + 5 <= maxDV ->
                [applyRing g ring outer]
            _ -> []
    in foldl' (processTree maxDV) st (children ++ ringChildren)

---------------------------------------------------------------------
-- Cross-check: spiral vs backtracking on all accepted pairs
-- Only run for small sizes where backtracking is tractable.
---------------------------------------------------------------------

crossCheck :: Int -> Map Int [DualGraph] -> IO Bool
crossCheck maxDV graphsBySize = do
    let sizes = [ (nv, gs) | (nv, gs) <- Map.toAscList graphsBySize
                            , nv <= maxDV ]
    if null sizes then do
        putStrLn "(skipped, no sizes in range)"
        return True
    else do
        results <- mapM checkSize sizes
        return (and results)
  where
    checkSize (nv, gs) = do
        let n = length gs
            pairs = [(i, j) | i <- [0..n-1], j <- [i+1..n-1]]
            keys = map generalSpiralKey gs
        errors <- mapM (\(i, j) -> do
            let gi = gs !! i
                gj = gs !! j
                ki = keys !! i
                kj = keys !! j
                spiralSame = ki == kj && ki /= Nothing
                btSame = isIsomorphic gi gj
            if spiralSame /= btSame
                then do
                    putStrLn $ "  MISMATCH at C" ++ show (carbonAtoms nv)
                             ++ " graphs " ++ show i ++ "," ++ show j
                             ++ ": spiral=" ++ show spiralSame
                             ++ " backtrack=" ++ show btSame
                    return False
                else return True
            ) pairs
        return (and errors)

---------------------------------------------------------------------
-- Known isomer counts
---------------------------------------------------------------------

knownCounts :: [(Int, Int)]
knownCounts =
    [ (20, 1), (24, 1), (26, 1), (28, 2), (30, 3)
    , (32, 6), (34, 6), (36, 15), (38, 17), (40, 40)
    , (42, 45), (44, 89), (46, 116), (48, 199), (50, 271)
    , (52, 437), (54, 580), (56, 924), (58, 1205), (60, 1812)
    , (62, 2385), (64, 3465), (66, 4478), (68, 6332), (70, 8149)
    , (72, 11190), (74, 14246), (76, 19151), (78, 24109), (80, 31924)
    ]

dualVerts :: Int -> Int
dualVerts c = c `div` 2 + 2

carbonAtoms :: Int -> Int
carbonAtoms dv = (dv - 2) * 2

---------------------------------------------------------------------
-- Main
---------------------------------------------------------------------

main :: IO ()
main = do
    args <- getArgs
    let maxDV = case args of
                  (s:_) -> read s
                  []    -> 22

    putStrLn $ "Generating fullerene duals up to " ++ show maxDV
             ++ " vertices (C" ++ show (carbonAtoms maxDV) ++ ")"
    putStrLn ""

    let st = generate maxDV
        results = gsGraphs st
        relevant = [ (dualVerts c, c, expected)
                   | (c, expected) <- knownCounts
                   , dualVerts c <= maxDV
                   ]

    putStrLn $ padR 8 "Dual V" ++ padR 8 "C_n" ++ padR 12 "Found"
             ++ padR 12 "Expected" ++ "Status"
    putStrLn $ replicate 52 '-'

    mapM_ (\(dv, cn, expected) -> do
        let found = length $ Map.findWithDefault [] dv results
            status | found == expected = "OK"
                   | found < expected  = "MISSING " ++ show (expected - found)
                   | otherwise         = "EXTRA " ++ show (found - expected)
        putStrLn $ padR 8 (show dv) ++ padR 8 ("C" ++ show cn)
                 ++ padR 12 (show found) ++ padR 12 (show expected)
                 ++ status
        ) relevant

    putStrLn ""
    let allOk = all (\(dv, _, expected) ->
                    length (Map.findWithDefault [] dv results) == expected) relevant
    if allOk
        then putStrLn "All counts match!"
        else putStrLn "MISMATCH detected"

    if gsFailed st > 0
        then putStrLn $ "WARNING: " ++ show (gsFailed st)
                      ++ " graphs where general spiral failed (should be 0)"
        else return ()

    -- Cross-check against backtracking for small sizes
    let crossCheckMax = min maxDV 22  -- up to C40
    putStrLn ""
    putStrLn $ "Cross-checking spiral vs backtracking (up to C"
             ++ show (carbonAtoms crossCheckMax) ++ ")..."
    ok <- crossCheck crossCheckMax results
    if ok
        then putStrLn "Cross-check passed."
        else putStrLn "Cross-check FAILED!"

padR :: Int -> String -> String
padR n s = s ++ replicate (max 0 (n - length s)) ' '
