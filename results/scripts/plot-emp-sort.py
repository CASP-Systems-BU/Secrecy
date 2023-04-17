import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

#num_rows = [1024, 5000, 10000, 20000, 30000, 40000, 60000, 80000, 100000]
#rows = {1024 : 0, 4096 : 1, 16384 : 2, 65536 : 3, 262144 : 4, 1048576 : 5, 4194304 : 6}

rows = {1024 : 0, 2048: 1, 4096 : 2, 8192: 3, 16384 : 4, 32768 : 5, 65536 : 6, 131072 : 7, 262144 : 8, 524288 : 9, 1048576 : 10, 2097152 : 11, 4194304 : 12}

num_runs_emp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
num_runs_sysx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

col_emp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
col_sysx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# try:
#     input_file_emp = sys.argv[1]
#     input_file_sysx = sys.argv[2]
#     title = sys.argv[3]
# except IndexError:
#     raise SystemExit(f"Usage: {sys.argv[0]} <input_file_emp> <input_file_sysx> <title>")

input_file_emp="./results-emp/sort.out"
input_file_sysx="./results-sysx/micro/sort.out"

# load emp results
with open(input_file_emp,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        index = rows[int(row[0])]
        col_emp[index] += float(row[1])
        num_runs_emp[index] +=1

# load sysX results
with open(input_file_sysx,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        index = rows[int(row[2])]
        col_sysx[index] += float(row[3])
        num_runs_sysx[index] +=1

# compute average values
for i in range(len(num_runs_emp)):
    if (num_runs_emp[i] > 0):
        col_emp[i] = round(col_emp[i]/num_runs_emp[i], 2)

for i in range(len(num_runs_sysx)):
    if (num_runs_sysx[i] > 0):
        col_sysx[i] = round(col_sysx[i]/num_runs_sysx[i], 2)

# plot
N = len(num_runs_emp)
fig, ax = plt.subplots()
ind = np.arange(7)    # x locations
width = 0.2         # the width of the bars

#x_pos = [i + width/2 for i, _ in enumerate(num_runs_emp)]

# 64K, 256K, 1M, 4M
col_emp2 = [e for (i, e) in enumerate(col_emp) if i>=6 ]#and i%2!=1]
col_sysx2 = [e for (i, e) in enumerate(col_sysx) if i>=6 ]#and i%2!=1]

p1 = ax.bar(ind - width, zero_to_nan(col_emp2), width, color='bisque', hatch="/")
p2 = ax.bar(ind, zero_to_nan(col_sysx2), width, color='darkred', hatch="-")

#ax.set_ylabel('Time (s)', fontsize=18)
ax.set_xlabel('Rows', fontsize=18)
ax.set_yscale('log')
ax.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax.set_ylim([0,100000])

ax.set_title("Sort", fontsize=18, fontweight='bold')
ax.set_xticks(ind - width/2)

#ax.set_xticklabels(['1K', '2K', '4K', '8K', '16K', '32K', '64K', '128K', '256K', '512K', '1M', '2M', '4M'])
ax.set_xticklabels(['64K', '128K', '256K', '512K', '1M', '2M', '4M'])
ax.legend((p1[0], p2[0]), ('EMP', 'Secrecy'), fontsize=16, loc='upper left')


def autolabel(rects1, rects2):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for (rec1, rec2) in zip(rects1, rects2):
        height1 = rec1.get_height()
        height2 = rec2.get_height()
        ax.annotate('{}x'.format(round(height1/height2, 2)),
                    xy=(rec1.get_x() + rec1.get_width(), height1),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=16, rotation=0, fontweight='bold')

autolabel(p1, p2)

ax.autoscale_view()
fig.tight_layout()
plt.show()