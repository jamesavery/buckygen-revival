-- | Benchmark exercising all Search schedule backends.
--
-- Usage: buckygen-bench-search <maxDV> [runs] [nWorkers]
--   maxDV:    max dual vertices (default 42 = C80)
--   runs:     number of timing runs per test (default 3)
--   nWorkers: threads for work-queue/hier/bfsdfs (default 10)
--
-- Run with: +RTS -N10 -RTS  (match nWorkers to capabilities)

module Main where

import Search
import qualified Data.Map.Strict as Map
import System.Environment (getArgs)
import System.IO (hFlush, stdout)
import Data.List (sort)
import Control.DeepSeq (NFData, rnf)
import GHC.Clock (getMonotonicTimeNSec)

main :: IO ()
main = do
    args <- getArgs
    let maxDV    = case args of { (s:_) -> read s; [] -> 42 }
        runs     = case drop 1 args of { (s:_) -> read s; [] -> 3 }
        nWorkers = case drop 2 args of { (s:_) -> read s; [] -> 10 }
        cfg = mkConfig maxDV
        cn  = (maxDV - 2) * 2
        depth = 9  -- standard fork depth for C80+

    putStrLn $ "=== Search Benchmark ==="
    putStrLn $ "    C" ++ show cn ++ " (dv=" ++ show maxDV
             ++ "), " ++ show runs ++ " runs, "
             ++ show nWorkers ++ " workers, depth=" ++ show depth ++ "\n"

    -- 1. Sequential DFS (baseline)
    putStrLn "--- 1. Sequential DFS ---"
    refCounts <- bench runs "  searchPure countBySize" $
        let !c = searchPure countBySize cfg in return c

    -- 2. Sequential BFS
    putStrLn "--- 2. Sequential BFS ---"
    bfsCounts <- bench runs "  search SeqBFS" $
        search countBySize SeqBFS cfg
    check "BFS" refCounts bfsCounts

    -- 3. Parallel Pure (flat-fork)
    putStrLn $ "--- 3. Parallel Pure (depth=" ++ show depth ++ ") ---"
    parCounts <- bench runs "  search (ParFlat)" $
        search countBySize (ParFlat depth) cfg
    check "ParFlat" refCounts parCounts

    -- 4. Parallel MutGraph (flat-fork)
    putStrLn $ "--- 4. Parallel MutGraph (depth=" ++ show depth ++ ") ---"
    mutCounts <- bench runs "  search (ParMut)" $
        search countBySize (ParMut depth) cfg
    check "ParMut" refCounts mutCounts

    -- 5. Work-queue, pure
    putStrLn $ "--- 5. Work-queue Pure (" ++ show nWorkers ++ "w, depth=" ++ show depth ++ ") ---"
    wqCounts <- bench runs "  search (WorkQ)" $
        search countBySize (WorkQ nWorkers depth) cfg
    check "WorkQ" refCounts wqCounts

    -- 6. Work-queue, MutGraph
    putStrLn $ "--- 6. Work-queue MutGraph (" ++ show nWorkers ++ "w, depth=" ++ show depth ++ ") ---"
    wqmCounts <- bench runs "  search (WorkQMut)" $
        search countBySize (WorkQMut nWorkers depth) cfg
    check "WorkQMut" refCounts wqmCounts

    -- 7. Hierarchical MutGraph
    let d1 = 4; d2 = 5
    putStrLn $ "--- 7. Hierarchical MutGraph (" ++ show nWorkers ++ "w, d1=" ++ show d1 ++ " d2=" ++ show d2 ++ ") ---"
    hierCounts <- bench runs "  search (HierMut)" $
        search countBySize (HierMut nWorkers d1 d2) cfg
    check "HierMut" refCounts hierCounts

    -- 8. BFS+DFS
    let bfsDepth = 7
    putStrLn $ "--- 8. BFS+DFS MutGraph (" ++ show nWorkers ++ "w, depth=" ++ show bfsDepth ++ ") ---"
    bdCounts <- bench runs "  search (BfsDfs)" $
        search countBySize (BfsDfs nWorkers nWorkers bfsDepth) cfg
    check "BfsDfs" refCounts bdCounts

    -- Final summary
    let total = sum (Map.elems refCounts)
    putStrLn $ "\n=== Summary: C" ++ show cn ++ " (" ++ show total ++ " isomers) ==="
    let allMatch = and [ bfsCounts == refCounts, parCounts == refCounts
                       , mutCounts == refCounts, wqCounts == refCounts
                       , wqmCounts == refCounts, hierCounts == refCounts
                       , bdCounts == refCounts ]
    putStrLn $ "  All results match: " ++ show allMatch
    if allMatch
        then putStrLn "  ALL CORRECT"
        else putStrLn "  *** MISMATCH ***"

check :: String -> Map.Map Int Int -> Map.Map Int Int -> IO ()
check label ref new = do
    let ok = ref == new
    putStrLn $ "  Match: " ++ show ok ++ if ok then "" else " *** " ++ label ++ " MISMATCH ***"
    putStrLn ""

-- | Run an IO action n times using wall-clock time, report median.
bench :: (NFData a) => Int -> String -> IO a -> IO a
bench n label action = do
    putStr $ label ++ ": "
    hFlush stdout
    times <- mapM (\_ -> do
        t0 <- getMonotonicTimeNSec
        !result <- action
        _ <- return $! rnf result
        t1 <- getMonotonicTimeNSec
        let ms = fromIntegral (t1 - t0) / 1e6 :: Double
        return (ms, result)
        ) [1..n]
    let ms = sort (map fst times)
        median = ms !! (n `div` 2)
        result = snd (last times)
    putStrLn $ showMs median ++ "  (runs: "
             ++ show (map showMs ms) ++ ")"
    hFlush stdout
    return result

showMs :: Double -> String
showMs ms
    | ms < 1000 = show (round ms :: Int) ++ "ms"
    | otherwise  = let s = ms / 1000 in show (fromIntegral (round (s * 10) :: Int) / 10 :: Double) ++ "s"
