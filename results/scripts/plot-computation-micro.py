import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [1000, 10000, 100000, 1000000, 10000000]
num_runs = [0, 0, 0, 0, 0]

col_equality = [0, 0, 0, 0, 0]
col_greater = [0, 0, 0, 0, 0]
col_rca = [0, 0, 0, 0, 0]

try:
    in_equality = sys.argv[1]
    in_greater = sys.argv[2]
    in_rca = sys.argv[3]
    title = sys.argv[4]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_equality> <input_greater> <input_rca>")

# load equality results
with open(in_equality,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == 'ASYNC-ARRAY':
            if int(row[1]) < 100000000:
                index = int(math.log10(int(row[1])/1000))
                col_equality[index] += float(row[2])
                num_runs[index] +=1

# compute average values
for i in range(len(col_equality)):
    if (num_runs[i] > 0):
        col_equality[i] = col_equality[i]/num_runs[i]
        col_equality[i] = round(col_equality[i], 3)
        num_runs[i] = 0

# load greater results
with open(in_greater,'r') as csvfile:
    rows = csv.reader(csvfile, delimiter='\t')
    for row in rows:
        if row[0].startswith('1'):
            index = int(math.log10(int(row[0])/1000))
            col_greater[index] += float(row[1])
            num_runs[index] +=1

# compute average values
for i in range(len(col_greater)):
    if (num_runs[i] > 0):
        col_greater[i] = col_greater[i]/num_runs[i]
        col_greater[i] = round(col_greater[i], 3)

# load rca results
with open(in_rca,'r') as csvfile:
    rows = csv.reader(csvfile, delimiter='\t')
    for row in rows:
        if row[0].startswith('1'):
            index = int(math.log10(int(row[0])/1000))
            col_rca[index] += float(row[1])
            num_runs[index] +=1

# compute average values
for i in range(len(col_rca)):
    if (num_runs[i] > 0):
        col_rca[i] = col_rca[i]/num_runs[i]
        col_rca[i] = round(col_rca[i], 3)

# plot
N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.20         # the width of the bars

p1 = ax.plot(num_rows, zero_to_nan(col_equality), color='green', marker='.', linestyle='dashed', linewidth=2, markersize=12)
p2 = ax.plot(num_rows, zero_to_nan(col_greater), color='darkred', marker='x', linewidth=2, markersize=12)
p3 = ax.plot(num_rows, zero_to_nan(col_rca) , color='darkblue', marker='+', linestyle='dotted', linewidth=2, markersize=12)

ax.set_ylabel('Time (s)', fontsize=16)
ax.set_xlabel('# ROWS', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

ax.set_title(title, fontsize=18)

ax.legend((p1[0], p2[0], p3[0]), ('Equality', 'Inequality', 'Addition'), fontsize=16)

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

#autolabel(p1)

ax.autoscale_view()
fig.tight_layout()
plt.show()