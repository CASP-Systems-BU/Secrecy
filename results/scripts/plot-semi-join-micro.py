import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math
import statistics

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [1024, 4096, 16384, 65536, 262144, 1048576]
num_runs = [0, 0, 0, 0, 0, 0]

col_left_fixed = {}
col_right_fixed = {}

try:
    input_file_semijoin = sys.argv[1]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_semi_join>")

# load exchange results
with open(input_file_semijoin,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if len(row) > 1 and 'IN-BATCH' in row[1]:
            print(row)
            new_key = int(row[2])
            new_key2 = int(row[3])
            if new_key < 1024:
                continue
            new_val = float(row[4])
            key = "{}-{}".format(new_key,new_key2)
            val =  col_left_fixed.get(key)
            if val is not None:
                col_left_fixed[key].append(new_val)
            else:
                l = [new_val]
                col_left_fixed.setdefault(key, l)

# print(col_left_fixed)

left_fixed = [0, 0, 0, 0, 0, 0]
right_fixed = [0, 0, 0, 0, 0, 0]

for key, val in col_left_fixed.items():
    keys = key.split('-')
    left = int(keys[0])
    right = int(keys[1])
    m = statistics.mean(val)
    if left==1024:
        left_fixed[num_rows.index(right)] = m
    if right==1024:
        right_fixed[num_rows.index(left)] = m

# print(left_fixed)
# print(right_fixed)

# plot
N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.35         # the width of the bars

p1 = ax.plot(num_rows, zero_to_nan(left_fixed), color='green', marker='.', linestyle='dashed', linewidth=2, markersize=12)
p2 = ax.plot(num_rows, zero_to_nan(right_fixed), color='darkred', marker='x', linewidth=2, markersize=12)

ax.set_ylabel('Time (s)', fontsize=20)
ax.set_xlabel('Rows', fontsize=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#ax.set_title(title, fontsize=18)
ax.legend((p1[0], p2[0]), ('Left fixed','Right fixed'), fontsize=16)
ax.set_xticklabels(['0', '0', '200K', '400K', '600K', '800K', '1M'])

ax.autoscale_view()
fig.tight_layout()
plt.show()
