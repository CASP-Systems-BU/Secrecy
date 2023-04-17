import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [1000, 10000, 100000, 1000000, 10000000, 100000000]
num_runs = [0, 0, 0, 0, 0, 0]

col_exchange = [0, 0, 0, 0, 0, 0]
col_reveal = [0, 0, 0, 0, 0, 0]

try:
    input_file_exchange = sys.argv[1]
    input_file_reveal = sys.argv[2]
    title = sys.argv[3]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_exchange> <input_file_reveal> <title>")

# load exchange results
with open(input_file_exchange,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == 'BATCHED':
            index = int(math.log10(int(row[1])/1000))
            col_exchange[index] += float(row[2])
            num_runs[index] +=1

# compute average values
for i in range(len(col_exchange)):
    col_exchange[i] = col_exchange[i]/num_runs[i]
    num_runs[i] = 0

# load reveal results
with open(input_file_reveal,'r') as csvfile:
    rows = csv.reader(csvfile, delimiter='\t')
    for row in rows:
        if row[0] == 'ASYNC':
            index = int(math.log10(int(row[1])/1000))
            col_reveal[index] += float(row[2])
            num_runs[index] +=1

# compute average values
for i in range(len(col_reveal)):
    col_reveal[i] = col_reveal[i]/num_runs[i]

# plot
N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.35         # the width of the bars

#p1 = ax.bar(ind - width/2, col_exchange, width, color='green')
#p2 = ax.bar(ind + width/2, col_reveal, width, color='salmon')

p1 = ax.plot(num_rows, zero_to_nan(col_exchange), color='green', marker='.', linestyle='dashed', linewidth=2, markersize=12)
p2 = ax.plot(num_rows, zero_to_nan(col_reveal), color='darkred', marker='x', linewidth=2, markersize=12)

ax.set_ylabel('Time (s)', fontsize=16)
ax.set_xlabel('# ROWS', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

ax.set_title(title, fontsize=18)
ax.legend((p1[0], p2[0]), ('Exchange', 'Reveal'), fontsize=16)


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