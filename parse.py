#!/usr/bin/env python

import os, sys, glob, subprocess
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# plot total total ns since we started keeping track (11/15/16):
totals = []
for i in glob.glob("*.dat"):
    cmd = "tail -n 1 " + i + " | awk '{print $2}'"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    totals.append(float(ps.communicate()[0].strip()))

totals = sorted(totals)
total_ns = max(totals)
first_ns = min(totals)
ns_since_start = max(totals)-min(totals)
days = len(totals)-1
avg_ns = ns_since_start / days
x = range(len(glob.glob("*.dat")))

plt.figure(1)
plt.xlabel("Days since November 15th, 2016")
plt.ylabel("Total ns")
plt.title("Avg. ns per day: " + str(avg_ns))
plt.plot(x, totals)
plt.savefig("overall.png")

# plot all projects separately
plt.figure(2)
labels = []
for i in glob.glob("*dat"):
    with open(i) as f:
        for line in f:
            if "total" not in line:
                labels.append(line.split()[0])

labels = list(set(labels))
labels = [[i] for i in labels]

for i in labels:
    for j in glob.glob("*.dat"):
        with open(j) as f:
            for line in f:
                if i[0] in line:                
                    columns = line.split()
                    i.append(float(columns[1]))
    
labels = [sorted(i) for i in labels]
days_since_start = max([len(i) for i in labels])

best_ns = max([i[-2]-i[-3] for i in labels])
for i in labels:
    if i[-2]-i[-3] == best_ns:
        best_job = i[-1] 

for i in range(len(labels)):
    starting_ns=labels[i][0]
    for j in range(len(labels[i][:-1])):
        labels[i][j] = labels[i][j]-starting_ns
    if labels[i][-1] != "total":
        plt.plot(range(days_since_start-len(labels[i]), days_since_start-1), labels[i][:-1])
plt.xlabel("Days since November 15th, 2016")
plt.ylabel("nanoseconds")
plt.title("Best job in last 24 hours: " + best_job + " with " + str(best_ns) + " ns in the last day.")
plt.savefig("all_proj.png")
    
