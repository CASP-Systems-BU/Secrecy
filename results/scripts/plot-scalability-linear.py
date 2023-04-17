import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

num_rows = [16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152]
rows = {16384 : 0, 32768 : 1, 65536 : 2, 131072 : 3, 262144 : 4, 524288 : 5, 1048576 : 6, 2097152 : 7}

q1 = [0, 0, 0, 0, 0, 0, 0, 0]
q2 = [0, 0, 0, 0, 0, 0, 0, 0]
q3 = [0, 0, 0, 0, 0, 0, 0, 0]
q4 = [0, 0, 0, 0, 0, 0, 0, 0]

num_runs_1 = [0, 0, 0, 0, 0, 0, 0, 0]
num_runs_2 = [0, 0, 0, 0, 0, 0, 0, 0]
num_runs_3 = [0, 0, 0, 0, 0, 0, 0, 0]
num_runs_4 = [0, 0, 0, 0, 0, 0, 0, 0]

input_comrb="./results-sysx/medical/q1.out"
input_cdiff="./results-sysx/medical/q2.out"
input_pwd="./results-sysx/password/scalability.out"
input_credit="./results-sysx/creditscore/scalability.out"

# load results
with open(input_comrb,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running'):
            if row[1] == 'Q1' and int(row[3]) == 256 and int(row[2]) > 16000:
                index = rows[int(row[2])]
                q1[index] += float(row[4])
                num_runs_1[index] +=1

with open(input_cdiff,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running') and int(row[2]) > 16000:
            if row[1] == 'Q2':
                index = rows[int(row[2])]
                q2[index] += float(row[3])
                num_runs_2[index] +=1

with open(input_pwd,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running') and int(row[0]) > 16000:
            index = rows[int(row[0])]
            q3[index] += float(row[1])
            num_runs_3[index] +=1


with open(input_credit,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if not row[0].startswith('Running') and int(row[0]) > 16000:
            index = rows[int(row[0])]
            q4[index] += float(row[1])
            num_runs_4[index] +=1


# compute average values
for i in range(len(num_rows)):
    if (num_runs_1[i] > 0):
        q1[i] = round(q1[i]/num_runs_1[i], 2)
    if (num_runs_2[i] > 0):
        q2[i] = round(q2[i]/num_runs_2[i], 2)
    if (num_runs_3[i] > 0):
        q3[i] = round(q3[i]/num_runs_3[i], 2)
    if (num_runs_4[i] > 0):
        q4[i] = round(q4[i]/num_runs_4[i], 2)

# plot
# 32K, 128K, 512K, 2M
q1_2 = [e for (i, e) in enumerate(q1) if i%2==1]
q2_2 = [e for (i, e) in enumerate(q2) if i%2==1]
q3_2 = [e for (i, e) in enumerate(q3) if i%2==1]
q4_2 = [e for (i, e) in enumerate(q4) if i%2==1]
labels = ['16K', '32K', '64K', '128K', '256K','512K', '1M', '2M']
labels_2 = [e for (i, e) in enumerate(labels) if i%2==1]

N = len(labels_2)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.17        # the width of the bars
width2 = 0.2
x_pos = [i + width/2 for i, _ in enumerate(num_rows)]

p1 = ax.bar(ind - 2*width2, zero_to_nan(q3_2), width, color='sandybrown', hatch="+")
p2 = ax.bar(ind - width2, zero_to_nan(q4_2), width, color='indianred', hatch="\\")
p3 = ax.bar(ind, zero_to_nan(q1_2), width, color='darkorange', hatch="/")
p4 = ax.bar(ind + width2, zero_to_nan(q2_2), width, color='firebrick', hatch="-")

#p1 = ax.plot(num_rows, zero_to_nan(q3), color='green', marker='.', linestyle='dashed', linewidth=2, markersize=12)
#p2 = ax.plot(num_rows, zero_to_nan(q4), color='darkblue', marker='+', linestyle='dotted', linewidth=2, markersize=12)
#p3 = ax.plot(num_rows, zero_to_nan(q1), color='darkorange', marker='x', linewidth=2, markersize=12)
#p4 = ax.plot(num_rows, zero_to_nan(q2), color='firebrick', marker='*', linewidth=2, markersize=12)


ax.set_ylabel('Time (s)', fontsize=20)
ax.set_xlabel('Rows', fontsize=20)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#ax.set_title(title, fontsize=18)
ax.set_xticks(ind - width/2)

ax.set_xticklabels(labels_2)
ax.legend((p1[0], p2[0], p3[0], p4[0]), ('Password Reuse', 'Credit Score', 'Comorbidity', 'Recurrent C. Diff.'), fontsize=16)

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