import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

### order-by, group-by, dostinct ###

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [1000, 10000, 100000, 1000000, 10000000]
num_rows_order = [1024, 2048, 4096, 8196, 16384, 32768, 65536, 131072, 262144, 524228, 1048576, 2097152]

num_runs_group = [0, 0, 0, 0, 0]
num_runs_distinct = [0, 0, 0, 0, 0]
num_runs_order = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

col_distinct = [0, 0, 0, 0, 0]
col_groupby = [0, 0, 0, 0, 0]
col_order = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

try:
    input_file_distinct = sys.argv[1]
    input_file_groupby = sys.argv[2]
    input_file_order = sys.argv[3]
    title = sys.argv[4]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_distinct> <input_file_groupby> <input_file_order> <title>")

# load distinct results
with open(input_file_distinct,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == '1' and int(row[2]) < 100000000:
            index = int(math.log10(int(row[2])/1000))
            col_distinct[index] += float(row[3])
            num_runs_distinct[index] +=1

# compute average values
for i in range(len(col_distinct)):
    col_distinct[i] = col_distinct[i]/num_runs_distinct[i]
    col_distinct[i] = round(col_distinct[i], 3)

# load groupby results
with open(input_file_groupby,'r') as csvfile:
    rows = csv.reader(csvfile, delimiter='\t')
    for row in rows:
        if row[0].startswith('1') and int(row[0]) < 100000000:
            index = int(math.log10(int(row[0])/1000))
            col_groupby[index] += float(row[2])
            num_runs_group[index] +=1

# compute average values
for i in range(len(col_groupby)):
    if (num_runs_group[i] > 0):
        col_groupby[i] = col_groupby[i]/num_runs_group[i]
        col_groupby[i] = round(col_groupby[i], 3)


# load order-by results
with open(input_file_order,'r') as csvfile:
    rows = csv.reader(csvfile, delimiter='\t')
    for row in rows:
        if not row[0].startswith('Running') and row[0] == '1':
            index = int(math.log2(int(row[2]))) - 10
            #print(row[2], index, len(col_order), len(row))
            col_order[index] += float(row[3])
            num_runs_order[index] +=1

# compute average values
for i in range(len(col_groupby)):
    if (num_runs_order[i] > 0):
        col_order[i] = col_order[i]/num_runs_order[i]
        col_order[i] = round(col_order[i], 3)

# plot
fig, ax = plt.subplots()

p1 = ax.plot(num_rows, zero_to_nan(col_distinct), color='green', marker='.', linestyle='dashed', linewidth=2, markersize=12)
p2 = ax.plot(num_rows, zero_to_nan(col_groupby), color='darkred', marker='x', linewidth=2, markersize=12)
p3 = ax.plot(num_rows_order, zero_to_nan(col_order) , color='darkblue', marker='+', linestyle='dotted', linewidth=2, markersize=12)

ax.set_ylabel('Time (s)', fontsize=16)
ax.set_xlabel('# ROWS', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

ax.set_title(title, fontsize=18)


ax.legend((p1[0], p2[0], p3[0]), ('Distinct', 'Group-By', 'Order-By'), fontsize=16)

ax.autoscale_view()
fig.tight_layout()
plt.show()