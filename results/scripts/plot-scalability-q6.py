import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

num_rows = [32768, 131072, 524288, 2097152, 8388608]
rows = {32768 : 0, 131072 : 1, 524288 : 2, 2097152 : 3, 8388608 : 4}

q = [0, 0, 0, 0, 0]

num_runs = [0, 0, 0, 0, 0]

input_q6="./results-sysx/tpch/tpch-q6.out"

# load results
with open(input_q6,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running'):
            if int(row[1]) in rows.keys():
                index = rows[int(row[1])]
                q[index] += float(row[2])
                num_runs[index] +=1


# compute average values
for i in range(len(num_rows)):
    if (num_runs[i] > 0):
        q[i] = round(q[i]/num_runs[i], 4)

# plot
labels = ['32K', '128K', '512K', '2M', '8M']

N = len(labels)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.3       # the width of the bars

p1 = ax.bar(ind-width/2, q, width, color='goldenrod', hatch=".")

ax.set_ylabel('Time (s)', fontsize=20)
ax.set_xlabel('Rows', fontsize=20)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#ax.set_title(title, fontsize=18)
ax.set_xticks(ind - width/2)

ax.set_xticklabels(labels)
ax.legend((p1[0],), ('TPC-H Q6',), fontsize=16)

ax.autoscale_view()
fig.tight_layout()
plt.show()