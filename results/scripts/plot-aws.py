import matplotlib.pyplot as plt
import numpy as np
import yaml
import math

import sys
import argparse

plot_files_directory = "../results-sysx/AWS/figure/"
output_files_directory = "../results-sysx/AWS/output/same-region/"


def getdata(experimentName):
    global input_files_directory, output_files_directory
    with open(output_files_directory + experimentName + ".yaml",
                'r') as output_stream:
        try:
            outputData = yaml.safe_load(output_stream)

            xData = outputData["verticalAxisTags"]
            titles = [
                list(experiment.keys())[0]
                for experiment in outputData["experiments"]
            ]
            yData = [[
                (run["inputOne"] + run["inputTwo"] + run["inputThree"])
                / 3.0 for run in list(experiment.values())[0]
            ] for experiment in outputData["experiments"]]

            return xData, titles, yData

        except yaml.YAMLError as exc:
            print(exc)


def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x == 0 else x for x in values]


def autolabel(rects1, rects2, ax, maxRects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for (rec1, rec2) in zip(rects1[0:maxRects], rects2[0:maxRects]):
        height1 = rec1.get_height()
        height2 = rec2.get_height()
        ax.annotate(
            '{}x'.format(int(round(height1 / height2, 0))),
            xy=(rec2.get_x() + rec2.get_width() / 2, height2),
            xytext=(0, 2),  # 2 points vertical offset
            textcoords="offset points",
            ha='center',
            va='bottom',
            fontsize=16,
            rotation=0,
            fontweight='bold')


def plotData(experimentName,
             subName,
             xData,
             titles,
             yData,
             xLog,
             yLog,
             bar,
             label,
             hatches,
             colors,
             lineStyle=None):

    global plot_files_directory

    # plot
    N = len(xData)
    fig, ax = plt.subplots()
    ind = np.arange(N)  # x locations
    width = 0.2  # the width of the bars
    ps = []

    if bar:
        for i in range(0, len(yData)):
            ps += [
                ax.bar(np.arange(len(yData[i])) + i * width,
                       zero_to_nan(yData[i]),
                       width,
                       color=colors[i],
                       hatch=hatches[i])
            ]
        ax.legend(ps, titles, fontsize=16)

        if label and len(yData) > 1:
            autolabel(ps[0], ps[1], ax, min(len(yData[0]), len(yData[1])))
    else:
        for i in range(0, len(yData)):
            pt, = ax.plot(xData[0:len(yData[i])],
                        zero_to_nan(yData[i]),
                        color=colors[i],
                        marker=hatches[i],
                        linestyle=lineStyle[i],
                        linewidth=2,
                        markersize=12)
            ps += [pt]
        ax.legend(ps, titles, fontsize=16)

    ax.set_ylabel('Time (s)', fontsize=16)

    if yLog:
        ax.set_yscale('log')
    
    if xLog:
        ax.set_xscale('log')

    ax.grid(True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=18)

    if bar:
        ax.set_xticks(ind + width / 2)
        ax.set_xticklabels(xData, rotation=10)
        for tick in ax.xaxis.get_major_ticks():
            tick.set_pad(-1)

    ax.set_title(experimentName + " " + subName, fontsize=15)
    ax.autoscale_view()
    fig.tight_layout()
    # plt.show()
    plt.savefig( plot_files_directory + experimentName + '.png')


def processExperiment(experiment_name,
                      subName,
                      xlabels,
                      xLog,
                      yLog,
                      bar,
                      label,
                      hatches,
                      colors,
                      lineStyle=None):

    print(experiment_name)

    xData, titles, yData = getdata(experiment_name)
    if xlabels == None:
        xlabels = xData
    plotData(experiment_name, subName, xlabels, titles, yData, xLog, yLog, bar, label, hatches, colors,
             lineStyle)


def getPlotConfig(configFile):
    with open(plot_files_directory + configFile + ".yaml",
              'r') as input_stream:
        try:
            inputData = yaml.safe_load(input_stream)

            return [list(plot.values())[0] for plot in inputData["plots"]]
        except yaml.YAMLError as exc:
            print(exc)


def main():

    my_parser = argparse.ArgumentParser(description='Plots Secrecy figures')

    my_parser.add_argument('-p',
                           '--path',
                           default="plot",
                           metavar='path',
                           type=str,
                           help='the path to list')

    args = my_parser.parse_args()

    plotConfigs = getPlotConfig(args.path)

    # xLog, yLog, bar, labelhatches, colors, lineStyle
    for plotConfig in plotConfigs:
        experimentName = subName = labels = xLog = bar = label = hatches = colors = lineStyle = None
        if "name" in plotConfig:
            experimentName = plotConfig["name"]
        
        if "subName" in plotConfig:
            subName = plotConfig["subName"]

        if "labels" in plotConfig:
            labels = plotConfig["labels"]

        if "xLog" in plotConfig:
            xLog = plotConfig["xLog"]

        if "yLog" in plotConfig:
            yLog = plotConfig["yLog"]

        if "bar" in plotConfig:
            bar = plotConfig["bar"]

        if "label" in plotConfig:
            label = plotConfig["label"]

        if "hatches" in plotConfig:
            hatches = plotConfig["hatches"]

        if "colors" in plotConfig:
            colors = plotConfig["colors"]

        if "lineStyle" in plotConfig:
            lineStyle = plotConfig["lineStyle"]

        processExperiment(experimentName, subName, labels, xLog, yLog, bar, label,
                          hatches, colors, lineStyle)


if __name__ == "__main__":
    main()