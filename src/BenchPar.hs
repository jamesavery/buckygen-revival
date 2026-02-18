module Main where

import GenForest
import qualified Data.Map.Strict as Map
import System.Environment (getArgs)
import System.IO (hFlush, stdout)
import GHC.Conc (getNumCapabilities)

-- Usage: buckygen-par <maxDV> <depth> <mode> [arg4]
--   maxDV:  max dual vertices (default 32 = C60)
--   depth:  fork depth (default 6). For hier mode: d1,d2 (e.g. "4,5")
--   mode:   pure | mut | wq | wqmut | hier | analyze (default pure)
--   arg4:   nWorkers for wq/wqmut/hier (default = capabilities)

main :: IO ()
main = do
    args <- getArgs
    let maxDV = case args of { (s:_) -> read s; [] -> 32 }
        cfg = mkConfig maxDV
        depthStr = case drop 1 args of { (s:_) -> s; [] -> "6" }
        depth = read depthStr :: Int   -- for non-hier modes
        mode  = case drop 2 args of
                    ("mut":_)     -> "mut"
                    ("pure":_)    -> "pure"
                    ("wq":_)      -> "wq"
                    ("wqmut":_)   -> "wqmut"
                    ("hier":_)    -> "hier"
                    ("bfsdfs":_)  -> "bfsdfs"
                    ("analyze":_) -> "analyze"
                    _             -> "pure"

    case mode of
        "analyze" -> do
            putStrLn "=== Tree Depth Profile ==="
            let dp = depthProfile cfg
            mapM_ (\(d,n) -> putStrLn $ "  depth " ++ show d ++ ": " ++ show n ++ " nodes") dp
            let total = sum (map snd dp)
            putStrLn $ "  total: " ++ show total ++ " nodes"
            hFlush stdout

            putStrLn $ "\n=== Subtree sizes at depth " ++ show depth ++ " ==="
            let sizes = subtreeSizes cfg depth
                nSubs = length sizes
                biggest = if null sizes then 0 else head sizes
                smallest = if null sizes then 0 else last sizes
                belowFork = sum sizes
                aboveFork = total - belowFork
            putStrLn $ "  " ++ show nSubs ++ " subtrees"
            putStrLn $ "  above fork: " ++ show aboveFork ++ " nodes"
            putStrLn $ "  below fork: " ++ show belowFork ++ " nodes"
            putStrLn $ "  biggest:  " ++ show biggest
            putStrLn $ "  smallest: " ++ show smallest
            putStrLn $ "  top 10: " ++ show (take 10 sizes)
            let avg = fromIntegral belowFork / fromIntegral nSubs :: Double
            putStrLn $ "  ratio biggest/avg: " ++ show (fromIntegral biggest / avg)

        "mut" -> do
            let counts = parMutCountBySize cfg depth
                isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (MutGraph, d=" ++ show depth ++ "): " ++ show isomers

        "wq" -> do
            nWorkers <- getWorkerCount args
            counts <- workQueueCountBySize cfg nWorkers depth
            let isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (wq pure, d=" ++ show depth
                     ++ ", " ++ show nWorkers ++ "w): " ++ show isomers

        "wqmut" -> do
            nWorkers <- getWorkerCount args
            counts <- wqMutCountBySize cfg nWorkers depth
            let isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (wq MutGraph, d=" ++ show depth
                     ++ ", " ++ show nWorkers ++ "w): " ++ show isomers

        "hier" -> do
            nWorkers <- getWorkerCount args
            let (d1, d2) = parseDepthPair depthStr
            counts <- hierarchicalMutCount cfg nWorkers d1 d2
            let isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (hier MutGraph, d1=" ++ show d1
                     ++ " d2=" ++ show d2 ++ ", " ++ show nWorkers ++ "w): "
                     ++ show isomers

        "bfsdfs" -> do
            nWorkers <- getWorkerCount args
            -- depth = fork depth; BFS workers = DFS workers = nWorkers
            counts <- parBfsDfsCount cfg nWorkers nWorkers depth
            let isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (BFS+DFS MutGraph, d=" ++ show depth
                     ++ ", " ++ show nWorkers ++ "w): " ++ show isomers

        _ -> do
            let counts = parCountBySize cfg depth
                isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (pure, d=" ++ show depth ++ "): " ++ show isomers

getWorkerCount :: [String] -> IO Int
getWorkerCount args = case drop 3 args of
    (s:_) -> return (read s)
    []    -> getNumCapabilities

-- | Parse "d1,d2" or just "d" (treated as d1=4, d2=d).
parseDepthPair :: String -> (Int, Int)
parseDepthPair s = case break (== ',') s of
    (a, ',':b) -> (read a, read b)
    (a, _)     -> (4, read a)      -- default d1=4
