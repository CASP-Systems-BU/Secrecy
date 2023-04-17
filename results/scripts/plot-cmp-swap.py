import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

num_rows = []

col_1 = []
col_2 = []
col_4 = []
col_8 = []

try:
    input_file = sys.argv[1]
    title = sys.argv[2]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_path> <title>")

with open(input_file,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if (not row[0].startswith("Running")):
            if (int(row[2]) >= 131072 and int(row[2]) < 33554432):
                if row[0] == '1':
                    num_rows.append(int(int(row[2])/1024))
                    col_1.append(float(row[3]))
                if row[0] == '2':
                    col_2.append(float(row[3]))
                if row[0] == '4':
                    col_4.append(float(row[3]))
                if row[0] == '8':
                    col_8.append(float(row[3]))

N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.2         # the width of the bars

x_pos = [i + width/2 for i, _ in enumerate(num_rows)]

p1 = ax.bar(ind - 2*width, col_1, width, color='green', hatch="/")
p2 = ax.bar(ind - width, col_2, width, color='salmon', hatch="-")
p3 = ax.bar(ind, col_4, width, color='skyblue', hatch="+")
p4 = ax.bar(ind + width, col_8, width, color='gray', hatch="x")

ax.set_ylabel('Time (s)', fontsize=16)
ax.set_xlabel('ROWS', fontsize=16)
ax.set_yscale('log')
ax.grid(True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

ax.set_title(title, fontsize=18)
ax.set_xticks(ind - width/2)
ax.set_xticklabels(['128K', '256K', '512K', '1M', '2M', '4M', '8M', '16M'])

ax.legend((p1[0], p2[0], p3[0], p4[0]), ('1 COLS', '2 COLS', '4 COLS', '8 COLS'), fontsize=16)


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(p2)

ax.autoscale_view()
fig.tight_layout()
plt.show()
