import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [1024, 4096, 8192, 16384, 32768]
rows = {1024 : 0, 4096 : 1, 8192 : 2, 16384 : 3, 32768 : 4}

q_asp = [0, 0, 0, 0, 0]
q4 = [0, 0, 0, 0, 0]
q13 = [0, 0, 0, 0, 0]

num_runs_1 = [0, 0, 0, 0, 0]
num_runs_2 = [0, 0, 0, 0, 0]
num_runs_3 = [0, 0, 0, 0, 0]

input_asp="./results-sysx/medical/q3.out"
input_q4="./results-sysx/tpch/tpch-q4.out"
input_q13="./results-sysx/tpch/tpch-q13.out"

# load results
with open(input_asp,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running'):
            if row[1] == 'Q3-BATCH' and int(row[2]) != 2048:
                index = rows[int(row[2])]
                q_asp[index] += float(row[5])
                num_runs_1[index] +=1

with open(input_q4,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running') and int(row[2]) > 512 and int(row[2]) != 2048:
            if row[1] == 'TPCH-Q4':
                index = rows[int(row[2])]
                q4[index] += float(row[4])
                num_runs_2[index] +=1

with open(input_q13,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running'):
            if row[1] == 'TPCH-Q13' and int(row[2]) != 2048:
                index = rows[int(row[2])]
                q13[index] += float(row[5])
                num_runs_3[index] +=1


# compute average values
for i in range(len(num_rows)):
    if (num_runs_1[i] > 0):
        q_asp[i] = round(q_asp[i]/num_runs_1[i], 2)
    if (num_runs_2[i] > 0):
        q4[i] = round(q4[i]/num_runs_2[i], 2)
    if (num_runs_3[i] > 0):
        q13[i] = round(q13[i]/num_runs_3[i], 2)

# plot
labels = ['1x', '4x', '8x', '16x', '32x']
N = len(labels)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.17        # the width of the bars
width2 = 0.2
x_pos = [i + width/2 for i, _ in enumerate(num_rows)]

p1 = ax.bar(ind - width2, zero_to_nan(q4), width, color='seagreen', hatch="+")
p2 = ax.bar(ind, zero_to_nan(q13), width, color='yellowgreen', hatch="\\")
p3 = ax.bar(ind + width2, zero_to_nan(q_asp), width, color='lightseagreen', hatch="/")

ax.set_ylabel('Time (s)', fontsize=20)
ax.set_xlabel('Scaling factor', fontsize=20)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#ax.set_title(title, fontsize=18)
ax.set_xticks(ind - width/2)

ax.set_xticklabels(labels)
ax.legend((p1[0], p2[0], p3[0]), ('TPC-H Q4', 'TPC-H Q13', 'Aspirin Count'), fontsize=16)

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=16)

#autolabel(p3)

ax.autoscale_view()
fig.tight_layout()
plt.show()