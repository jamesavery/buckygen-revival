module Main where

import GenForest
import qualified Data.Map.Strict as Map
import System.Environment (getArgs)
import System.IO (hFlush, stdout)
import GHC.Conc (getNumCapabilities)

main :: IO ()
main = do
    args <- getArgs
    let maxDV = case args of { (s:_) -> read s; [] -> 32 }
        cfg = mkConfig maxDV
        depth = case drop 1 args of { (s:_) -> read s; [] -> 4 }
        mode  = case drop 2 args of
                    ("mut":_)     -> "mut"
                    ("pure":_)    -> "pure"
                    ("wq":_)      -> "wq"
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
            putStrLn $ "Total isomers (MutGraph): " ++ show isomers

        "wq" -> do
            nCaps <- getNumCapabilities
            -- depth arg repurposed as worker count; default to capabilities
            let nWorkers = if length args >= 2 then depth else nCaps
            counts <- workQueueCountBySize cfg nWorkers
            let isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (work-queue, " ++ show nWorkers ++ " workers): " ++ show isomers

        _ -> do
            let counts = parCountBySize cfg depth
                isomers = sum (Map.elems counts)
            putStrLn $ "Total isomers (pure): " ++ show isomers
