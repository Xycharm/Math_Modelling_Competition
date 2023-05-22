import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import json
from matplotlib.ticker import FormatStrFormatter


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def export2json(filename, sol_best):
    coord = []
    for line in open(filename, "r").readlines():
        x = line.strip("\r\n").split(",")
        coord.append({'lat': x[0], 'lng': x[1]})

    export_data = []
    for i in range(len(sol_best)):
        export_data.append(coord[int(sol_best[i])])

    file = open("path.json", 'w')
    file.write(json.dumps(export_data))
    file.close()


def sum_distmat(p, distmat):
    dist = 0
    num_location = p.shape[0]
    for i in range(num_location - 1):
        dist += distmat[p[i]][p[i + 1]]
    dist += distmat[p[0]][p[num_location - 1]]
    return dist


def get_distmat(p):
    num_location = p.shape[0]
    # 1 degree of lat/lon ~ 111km = 111000m in Taiwan
    p *= 111000
    distmat = np.zeros((num_location, num_location))
    for i in range(num_location):
        for j in range(i, num_location):
            distmat[i][j] = distmat[j][i] = np.linalg.norm(p[i] - p[j])
    return distmat


def swap(sol_new):
    while True:
        n1 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        n2 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        if n1 != n2:
            break
    sol_new[n1], sol_new[n2] = sol_new[n2], sol_new[n1]
    return sol_new


def reverse(sol_new):
    while True:
        n1 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        n2 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        if n1 != n2:
            break
    sol_new[n1:n2] = sol_new[n1:n2][::-1]

    return sol_new


def transpose(sol_new):
    while True:
        n1 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        n2 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        n3 = np.int(np.floor(np.random.uniform(0, sol_new.shape[0])))
        if n1 != n2 != n3 != n1:
            break
    # Let n1 < n2 < n3
    n1, n2, n3 = sorted([n1, n2, n3])

    # Insert data between [n1,n2) after n3
    tmplist = sol_new[n1:n2].copy()
    sol_new[n1: n1 + n3 - n2 + 1] = sol_new[n2: n3 + 1].copy()
    sol_new[n3 - n2 + 1 + n1: n3 + 1] = tmplist.copy()
    return sol_new


def accept(cost_new, cost_current, T):
    # If new cost better than current, accept it
    # If new cost not better than current, accept it by probability P(dE)
    # P(dE) = exp(dE/(kT)), defined by Metropolis
    return (cost_new < cost_current or
            np.random.rand() < np.exp(-(cost_new - cost_current) / T))


def save_sqlite(payloads):
    conn = sqlite3.connect('data/tsp.db')
    c = conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS TSP (costs REAL, route TEXT, markov_step INTEGER) ")
    c.execute('INSERT INTO TSP VALUES (?,?,?)', payloads)
    conn.commit()
    conn.close()
