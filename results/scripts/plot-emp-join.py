import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [10000, 20000, 30000, 40000, 60000, 80000, 100000]
rows = {10000 : 0, 20000 : 1, 30000 : 2, 40000 : 3, 60000 : 4, 80000 : 5, 100000 : 6}

num_runs_emp = [0, 0, 0, 0, 0, 0, 0]
num_runs_sysx = [0, 0, 0, 0, 0, 0, 0]

col_emp = [0, 0, 0, 0, 0, 0, 0]
col_sysx = [0, 0, 0, 0, 0, 0, 0]

# try:
#     input_file_emp = sys.argv[1]
#     input_file_sysx = sys.argv[2]
#     title = sys.argv[3]
# except IndexError:
#     raise SystemExit(f"Usage: {sys.argv[0]} <input_file_emp> <input_file_sysx> <title>")

input_file_emp="./results-emp/join_all.out"
input_file_sysx="./results-sysx/micro/join.out"

# load emp results
with open(input_file_emp,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if int(row[0]) > 5000:
            index = rows[int(row[0])]
            col_emp[index] += float(row[2])
            num_runs_emp[index] +=1

# load sysX results
with open(input_file_sysx,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if int(row[2]) > 5000:
            index = rows[int(row[2])]
            col_sysx[index] += float(row[3])
            num_runs_sysx[index] +=1

# compute average values
for i in range(len(num_rows)):
    if (num_runs_emp[i] > 0):
        col_emp[i] = round(col_emp[i]/num_runs_emp[i], 2)

for i in range(len(num_runs_sysx)):
    if (num_runs_sysx[i] > 0):
        col_sysx[i] = round(col_sysx[i]/num_runs_sysx[i], 2)

# plot
N = len(num_rows)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.2         # the width of the bars

x_pos = [i + width/2 for i, _ in enumerate(num_rows)]

#col_emp[3] = 100000

emp_over = [0, 0, 0, 0, 100000, 100000, 100000]

p1 = ax.bar(ind - width, zero_to_nan(col_emp), width, color='bisque', hatch="/")
p2 = ax.bar(ind, zero_to_nan(col_sysx), width, color='darkred', hatch="-")
p3 = ax.bar(ind - width, zero_to_nan(emp_over), width, color='bisque')

ax.set_ylabel('Time (s)', fontsize=18)
ax.set_xlabel('Rows per input', fontsize=18)
ax.set_yscale('log')
ax.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax.set_ylim([0,100000])

ax.set_title("Equi-Join", fontsize=18, fontweight='bold')
ax.set_xticks(ind - width/2)

ax.set_xticklabels(['10K', '20K', '30K', '40K', '60K', '80K', '100K'])
ax.legend((p1[0], p2[0]), ('EMP', 'Secrecy'), fontsize=16)

#(xmin, xmax) = ax.get_xlim()
#ax.hlines(y=43200, xmin=2.5, xmax=xmax, linewidth=2, color='r', linestyle='--')
#ax.annotate('12h', xy=(2.2, 42000), ha='center', va='bottom', color='r', fontstyle='italic', fontweight='bold', fontsize=12)

ax.annotate('over 15h', xy=(p1[4].get_x()+width/2, 1100), ha='center', va='bottom', color='darkslategrey', fontstyle='italic', fontweight='bold', fontsize=14, rotation=90)
ax.annotate('over 15h', xy=(p1[5].get_x()+width/2, 1100), ha='center', va='bottom', color='darkslategrey', fontstyle='italic', fontweight='bold', fontsize=14, rotation=90)
ax.annotate('over 15h', xy=(p1[6].get_x()+width/2, 1100), ha='center', va='bottom', color='darkslategrey', fontstyle='italic', fontweight='bold', fontsize=14, rotation=90)

def autolabel(rects1, rects2):
    for (rec1, rec2) in zip(rects1[0:4], rects2[0:4]):
        height1 = rec1.get_height()
        height2 = rec2.get_height()
        ax.annotate('{}x'.format(round(height1/height2, 2)),
                    xy=(rec2.get_x() + rec2.get_width(), height2),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=16, rotation=0, fontweight='bold')

autolabel(p1, p2)

def label_bar(rect, hours):
    """Attach a text label above each bar in *rects*, displaying its height."""
    height = rect.get_height()
    ax.annotate(hours,
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=14, fontweight='bold', color='darkred')

#label_bar(p1[3], '14.3h')
label_bar(p2[6], '12h')

ax.autoscale_view()
fig.tight_layout()
plt.show()