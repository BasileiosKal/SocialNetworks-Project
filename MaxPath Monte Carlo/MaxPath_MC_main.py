from MaxPath_MC_utils import maxPath_MonteCarlo
from utils import *
import multiprocessing as mp
import sqlite3
import time
import datetime
from tqdm import tqdm


if __name__ == "__main__":
    CATEGORY_TABLE = "SG_2"    # Table corresponding to the category we will create a graph for
    RUNS = 1000                # Monte Carlo simulations
    THRESHOLD = 0.01           # Threshold of Monte Carlo algorithm
    SUBGRAPH_SIZE = 10000
    MAX_PATH_LENGTH = 10       # The maximum distance from the source node that a user can be and be
                               #  considered to be in the source nodes reachable set.

    pbar = tqdm(total=10)
    def collect_result(result):
        global results
        global pbar
        results.append(result)
        pbar.update()

    # ====================================================================================================== #
    # Create the graphs
    # ====================================================================================================== #
    SGdb = sqlite3.connect('C:/Users/gauss/Downloads/SG_2.db', uri=True)
    # Create main graph
    MainGraph = get_category_graph(CATEGORY_TABLE, SGdb)
    SGdb.close()
    # Get subgraph
    supgraph = get_supgraph(MainGraph, supgraph_size=SUBGRAPH_SIZE)

    # Create users list here for faster iteration later
    users_in_graph = [uid for uid in supgraph.nodes]

    print("Number of processors: ", mp.cpu_count())


    # ====================================================================================================== #
    # Run Monte Carlo
    # ====================================================================================================== #
    print("Starting at: ", datetime.datetime.now())

    # ************************************************************** #
    # Transfer the Graph from Networkx to dictionaries to be passed
    # on the Monte_carlo_parallel function. That is necessary as a
    # result of the way the multiprocessing module works
    # ************************************************************** #
    # Adjacency list of the sup graph
    dict_subgraph = {}
    for node in users_in_graph:
        dict_subgraph[node] = list(supgraph.neighbors(node))

    # puvi (nested) dictionary.
    dict_puvi = {}
    for user1 in users_in_graph:
        dict_puvi[user1] = {}
        for user2 in supgraph[user1]:
            dict_puvi[user1][user2] = supgraph[user1][user2]['weight']


    start_time = time.time()
    print("Running Monte Carlo...", end=" ")

    pool = mp.Pool(mp.cpu_count())

    prob_threshold = RUNS * THRESHOLD

    results = []
    for user in users_in_graph:
        pool.apply_async(maxPath_MonteCarlo, args=(user, RUNS, dict_puvi, dict_subgraph, prob_threshold, MAX_PATH_LENGTH), callback = collect_result)


    pool.close()
    pool.join()

    end_time = time.time()

    print()
    print("Done")
    print("===================================")
    print("Running Time: ", end_time - start_time)

    time_in_min = end_time - start_time
    time_in_min = time_in_min / 60

    print("Required time in min: ", time_in_min)
    print("Required time in hours: ", time_in_min / 60)

    # ====================================================================================================== #
    # RResults Stats
    # ====================================================================================================== #
    print("=======================================================================================================")
    print("=======================================================================================================")
    print("=======================================================================================================")
    ReliableSets2 = {}
    for res in results:
        user = res[0]
        RS = res[1]
        ReliableSets2[user] = RS

    get_reliable_sets_stats(ReliableSets2)
