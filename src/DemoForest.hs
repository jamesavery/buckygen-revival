module Main where

import GenForest
import Seeds (DualGraph(..))
import Expansion (deg, nbrAt)
import Canonical (Automorphism(..))
import qualified Data.Map.Strict as Map
import Data.List (foldl', sort, group, intercalate)
import System.IO (hFlush, stdout)
import Control.DeepSeq (NFData(..))

main :: IO ()
main = do
    let maxDV = 32  -- C60
        cfg = mkConfig maxDV

    putStrLn "=== Forest Traversal Schemes (C60) ===\n"

    -- 1. Sequential DFS (baseline)
    putStrLn "1. Sequential DFS (dfsCountBySize)"
    let dfsCounts = dfsCountBySize cfg
    putStrLn $ "   Total: " ++ show (sum (Map.elems dfsCounts))
    hFlush stdout

    -- 2. BFS level-by-level
    putStrLn "\n2. BFS level-by-level (bfsByLevel)"
    let levels = bfsByLevel cfg
        bfsCounts' = bfsCountBySize cfg
        nLevels = length levels
    putStrLn $ "   Tree depth: " ++ show nLevels ++ " levels"
    putStrLn $ "   Level sizes: " ++ show (map length levels)
    putStrLn $ "   Total: " ++ show (sum (Map.elems bfsCounts'))
    putStrLn $ "   Match DFS: " ++ show (bfsCounts' == dfsCounts)
    hFlush stdout

    -- 3. Generic fold (reimplements dfsCountBySize)
    putStrLn "\n3. Generic fold (foldForest)"
    let foldCounts = foldForest
            (\m g _ -> Map.insertWith (+) (numVertices g) 1 m) Map.empty cfg
    putStrLn $ "   Total: " ++ show (sum (Map.elems foldCounts))
    putStrLn $ "   Match DFS: " ++ show (foldCounts == dfsCounts)
    hFlush stdout

    -- 4. Lazy DFS stream
    putStrLn "\n4. Lazy DFS stream (allGraphs)"
    let stream = allGraphs cfg
        streamCounts = foldl' (\m (g, _) ->
            Map.insertWith (+) (numVertices g) 1 m) Map.empty stream
    putStrLn $ "   Total: " ++ show (sum (Map.elems streamCounts))
    putStrLn $ "   Match DFS: " ++ show (streamCounts == dfsCounts)
    -- Demonstrate laziness: take first 10
    let first10 = take 10 stream
    putStrLn $ "   First 10 sizes: " ++ show (map (numVertices . fst) first10)
    hFlush stdout

    -- 5. Parallel flat-fork (pure)
    putStrLn "\n5. Parallel flat-fork (parCountBySize, depth=4)"
    let parCounts = parCountBySize cfg 4
    putStrLn $ "   Total: " ++ show (sum (Map.elems parCounts))
    putStrLn $ "   Match DFS: " ++ show (parCounts == dfsCounts)
    hFlush stdout

    -- 6. Parallel flat-fork (MutGraph)
    putStrLn "\n6. Parallel flat-fork + MutGraph (parMutCountBySize, depth=4)"
    let mutCounts = parMutCountBySize cfg 4
    putStrLn $ "   Total: " ++ show (sum (Map.elems mutCounts))
    putStrLn $ "   Match DFS: " ++ show (mutCounts == dfsCounts)
    hFlush stdout

    -- 7. Target-size generation
    putStrLn "\n7. Target-size generation (graphsOfSize)"
    let targets = [12, 16, 22, 27, 32]
    mapM_ (\t -> do
        let gs = graphsOfSize cfg t
            cn = (t - 2) * 2
            expected = Map.findWithDefault 0 t dfsCounts
        putStrLn $ "   C" ++ show cn ++ " (dv=" ++ show t ++ "): "
                 ++ show (length gs) ++ " isomers"
                 ++ (if length gs == expected then " OK" else " MISMATCH (expected " ++ show expected ++ ")")
        ) targets
    hFlush stdout

    -- 8. Parallel map: compute automorphism group sizes
    putStrLn "\n8. Parallel map: automorphism group size distribution"
    let autSizes = parMapForest cfg 4
            (\_ auts -> length auts) :: [Int]
        autDist = Map.fromListWith (+) [(s, 1 :: Int) | s <- autSizes]
    putStrLn $ "   Graphs processed: " ++ show (length autSizes)
    putStrLn $ "   |Aut| distribution: "
            ++ intercalate ", " [show s ++ ":" ++ show n | (s, n) <- Map.toAscList autDist]
    hFlush stdout

    -- 9. Parallel map: degree-5 vertex distance histogram
    putStrLn "\n9. Parallel map: min degree-5 distance per graph"
    let minDists = parMapForest cfg 4
            (\g _ -> minDeg5Distance g) :: [Int]
        distHist = Map.fromListWith (+) [(d, 1 :: Int) | d <- minDists]
    putStrLn $ "   Distance histogram: "
            ++ intercalate ", " ["d=" ++ show d ++ ":" ++ show n
                                | (d, n) <- Map.toAscList distHist]
    hFlush stdout

    -- 10. Generic fold: compute max automorphism group size per graph size
    putStrLn "\n10. Generic fold: max |Aut| per C_n"
    let maxAuts = foldForest
            (\m g auts -> Map.insertWith max (numVertices g) (length auts) m)
            Map.empty cfg
    putStrLn $ "   " ++ intercalate ", "
            ["C" ++ show ((nv - 2) * 2) ++ ":" ++ show a
            | (nv, a) <- Map.toAscList maxAuts]
    hFlush stdout

    -- Summary
    putStrLn "\n=== All cross-checks ==="
    let checks = [ ("BFS == DFS", bfsCounts' == dfsCounts)
                 , ("fold == DFS", foldCounts == dfsCounts)
                 , ("stream == DFS", streamCounts == dfsCounts)
                 , ("par == DFS", parCounts == dfsCounts)
                 , ("mut == DFS", mutCounts == dfsCounts)
                 ]
    mapM_ (\(name, ok) ->
        putStrLn $ "   " ++ (if ok then "PASS" else "FAIL") ++ ": " ++ name) checks
    let allOk = all snd checks
    putStrLn $ if allOk then "\nALL SCHEMES AGREE" else "\nMISMATCH DETECTED"

-- | Minimum distance between any two degree-5 vertices in the graph.
-- Uses simple BFS from each degree-5 vertex (graphs are small).
minDeg5Distance :: DualGraph -> Int
minDeg5Distance g =
    let d5 = degree5 g
        nv = numVertices g
    in minimum [ bfsDist g nv a b | a <- d5, b <- d5, a < b ]

bfsDist :: DualGraph -> Int -> Int -> Int -> Int
bfsDist g nv src tgt = go [src] (Map.singleton src 0)
  where
    go [] _ = nv  -- unreachable
    go (u:queue) visited =
        let d = visited Map.! u
        in if u == tgt then d
           else let nbrs = [nbrAt g u i | i <- [0 .. deg g u - 1]]
                    new = [v | v <- nbrs, not (Map.member v visited)]
                    visited' = foldl' (\m v -> Map.insert v (d+1) m) visited new
                in go (queue ++ new) visited'
