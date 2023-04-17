import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

num_rows = []

time_sync = []
time_async = []
time_batched = []
time_batched_inter = []

try:
    input_file = sys.argv[1]
    title = sys.argv[2]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_path> <title>")

with open(input_file,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == 'SYNC':
            num_rows.append(int(int(row[1])/1000))
            time_sync.append(float(row[2]))
        if row[0] == 'ASYNC':
            time_async.append(float(row[2]))
        if row[0] == 'ASYNC-ARRAY':
            time_batched.append(float(row[2]))
        if row[0] == 'ASYNC-INTER':
            time_batched_inter.append(float(row[2]))

N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.2         # the width of the bars

x_pos = [i + width/2 for i, _ in enumerate(num_rows)]

p1 = ax.bar(ind - 2*width, time_sync, width, color='green')
p2 = ax.bar(ind - width, time_async, width, color='salmon')
p3 = ax.bar(ind, time_batched, width, color='darkblue')
p4 = ax.bar(ind + width, time_batched_inter, width, color='gray')

ax.set_ylabel('Time (s)')
ax.set_xlabel('K ROWS')
ax.set_yscale('log')

ax.set_title(title)
ax.set_xticks(ind - width/2)
ax.set_xticklabels(num_rows)

ax.legend((p1[0], p2[0], p3[0], p4[0]), ('Element-wise SYNC', 'Element-wise ASYNC', 'Batched', 'Batched-interleaved'))


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(p3)

ax.autoscale_view()
fig.tight_layout()
plt.show()
