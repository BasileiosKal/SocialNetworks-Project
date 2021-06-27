from collections import defaultdict
import time
import numpy as np



def reachable_set(graph, source):
    """ Get the reachable set of the source user in the graph using DFS """
    R = []
    R = find_path(graph, source, R)
    return R


def find_path(graph, source, R):
    """ Helping function to run DFS on the Graph """
    R.append(source)
    if source not in graph.keys():  # Use graph.nodes if graph uses networkx
        return R

    for node in graph[source]:
        if node not in R:
            find_path(graph, node, R)
    return R


def Monte_Carlo(userid, T, dict_puvi, prob_threshold=0.02):
    """Run <T> Monte Carlo simulations, from the source node <userid>, in the whole
    graph defined by <dict_puvi>"""
    edges_with_weights = []
    for user1 in dict_puvi.keys():
        for user2 in dict_puvi[user1]:
            edges_with_weights.append((user1, user2, dict_puvi[user1][user2]))


    rnd_list = np.random.uniform(0, 1, size=((len(edges_with_weights)) * T,))
    i = 0

    c = defaultdict(lambda: 0)  # dictionary to store the counters for each userid
    for h in range(T):  # T = Number of samples

        G_tmp = {}  # G_tmp = the graph  that is created in h iteration as a DICTIONARY

        for user1, user2, p_uvi in edges_with_weights:

            rnd = rnd_list[i]
            i += 1
            if rnd <= p_uvi:
                # G_tmp.add_edge(user1, user2)  # To be used if G_tmp uses networkx
                if user1 not in G_tmp.keys():
                    G_tmp[user1] = [user2]
                else:
                    G_tmp[user1].append(user2)

        R_tmp = reachable_set(G_tmp, userid)

        for v in R_tmp:
            c[v] += 1

    Ru = []
    for v in c.keys():
        if c[v] >= prob_threshold:
            Ru.append(v)
    return (userid, Ru)

