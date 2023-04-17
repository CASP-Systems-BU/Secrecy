import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

q_base = [0, 0, 0, 0, 0, 0, 0, 0]
q_opt = [0, 0, 0, 0, 0, 0, 0, 0]

num_runs_base = [0, 0, 0, 0, 0, 0, 0, 0]
num_runs_opt = [0, 0, 0, 0, 0, 0, 0, 0]

#try:
#    input_file = sys.argv[1]
#    title = sys.argv[2]
#except IndexError:
#    raise SystemExit(f"Usage: {sys.argv[0]} <input_file_in> <title>")

input_med="./results-sysx/medical/baseline.out"
input_psw="./results-sysx/password/baseline.out"
input_credit="./results-sysx/creditscore/qcredit100.out"
input_tpch="./results-sysx/tpch/baseline.out"

# load med results
with open(input_med,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0].startswith('Q1-BASELINE'):
            q_base[0] += float(row[3])
            num_runs_base[0] +=1
        if row[0] == ('Q1'):
            q_opt[0] += float(row[3])
            num_runs_opt[0] +=1
        if row[0].startswith('Q2-baseline'):
            q_base[1] += float(row[2])
            num_runs_base[1] +=1
        if row[0] == ('Q2'):
            q_opt[1] += float(row[2])
            num_runs_opt[1] +=1
        if row[0].startswith('Q3-BASELINE'):
            q_base[2] += float(row[3])
            num_runs_base[2] +=1
        if row[0] == ('Q3-BATCH'):
            q_opt[2] += float(row[4])
            num_runs_opt[2] +=1


# load passwd results
with open(input_psw,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0].startswith('BASELINE'):
            q_base[3] += float(row[2])
            num_runs_base[3] +=1
        if row[0] == ('SYSX'):
            q_opt[3] += float(row[2])
            num_runs_opt[3] +=1

# load credit results
with open(input_credit,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0].startswith('BASELINE') and not (row[1].startswith('Running')):
            q_base[4] += float(row[2])
            num_runs_base[4] +=1
        if row[0].startswith('1024'):
            q_opt[4] += float(row[1])
            num_runs_opt[4] +=1

# load TPC-H Q4 results
with open(input_tpch,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == ('TPCH-Q4-BASELINE'):
            q_base[5] += float(row[3])
            num_runs_base[5] +=1
        if row[0] == ('TPCH-Q4'):
            q_opt[5] += float(row[3])
            num_runs_opt[5] +=1

# load TPC-H Q6 results
with open(input_tpch,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == ('TPCH-Q6-BASELINE'):
            q_base[6] += float(row[2])
            num_runs_base[6] +=1
        if row[0] == ('TPCH-Q6'):
            q_opt[6] += float(row[2])
            print(float(row[2]))
            num_runs_opt[6] +=1

# load TPC-H Q13 results
with open(input_tpch,'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        if row[0] == ('TPCH-Q13-BASELINE'):
            q_base[7] += float(row[3])
            num_runs_base[7] +=1
        if row[0] == ('TPCH-Q13'):
            q_opt[7] += float(row[4])
            num_runs_opt[7] +=1


# compute average values
for i in range(len(q_base)):
    q_base[i] = round(q_base[i]/num_runs_base[i], 4)
    q_opt[i] = round(q_opt[i]/num_runs_opt[i], 4)

# plot
N = len(q_base)
fig, ax = plt.subplots()
ind = np.arange(N)    # x locations
width = 0.1        # the width of the bars

p1 = ax.bar(ind, zero_to_nan(q_base), width, color='lightgray', hatch="/")
p2 = ax.bar(ind+width, zero_to_nan(q_opt), width, color='darkred', hatch="-")

ax.set_ylabel('Time (s)', fontsize=16)
ax.set_yscale('log')
ax.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=18)

ax.set_xticks(ind + width/2)
ax.set_xticklabels(['Comorb', 'RC. Diff', 'Aspirin', 'Pwd', 'Credit', 'TPC-H Q4', 'TPC-H Q6', 'TPC-H Q13'], rotation=10)
for tick in ax.xaxis.get_major_ticks():
    tick.set_pad(-1)
ax.legend((p1[0], p2[0]), ('No-opt', 'Secrecy'), fontsize=16)

def autolabel(rects1, rects2):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for (rec1, rec2) in zip(rects1, rects2):
        height1 = rec1.get_height()
        height2 = rec2.get_height()
        ax.annotate('{}x'.format(int(round(height1/height2, 0))),
                    xy=(rec2.get_x() + rec2.get_width() / 2, height2),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=16, rotation=0, fontweight='bold')

autolabel(p1, p2)

ax.autoscale_view()
fig.tight_layout()
plt.show()

