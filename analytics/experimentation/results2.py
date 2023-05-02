
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0, 2, 9, 10542, 1],
    [1, 9.0, 3, 17974, 1],
    [0, 3.0, 4.65, 11986, 2],
    [1, 7.18, 0.65, 17237, 2],
    [0, 13.15, 3.8, 24494, 3],
    [1, 7.5,0.15 , 18129, 3]
])

# Extract the unique theta and battery levels
theta_levels = np.unique(data[:, 0])
battery_levels = np.unique(data[:, 1])

# Define the bar chart parameters
bar_width = 0.35
space_width = 0.1

# Get the number of policies and calculate the total bar width
num_policies = len(np.unique(data[:, 4]))
total_bar_width = num_policies * len(battery_levels) * (bar_width + space_width)

data_plot = []

# Iterate over the battery levels and theta levels to create the bars
for i, theta_level in enumerate(theta_levels):
    for j, battery_level in enumerate(battery_levels):
        for k, policy in enumerate(np.unique(data[:, 4])):
            values = data[(data[:, 0] == theta_level) & (data[:, 1] == battery_level) & (data[:, 4] == policy), 2]
            if (not np.any(values)): continue
            if policy == 1:
                if theta_level == 0.0:
                    local_data = ["WFH T1", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T1", "bidirectional", values[0]]
            elif policy == 2:
                if theta_level == 0.0:
                    local_data = ["WFH T2", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T2", "bidirectional", values[0]]
            elif policy == 3:
                if theta_level == 0.0:
                    local_data = ["WFH T3", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T3", "bidirectional", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'Battery (kWh)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['Battery (kWh)'] = df['Battery (kWh)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Battery (kWh)", hue="Connectivity", hue_order=["unidirectional","bidirectional"], order=["WFH T1", "WFH T2", "WFH T3"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 10)
sns_plot.set_xlabel('')



plt.show()













import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0, 2, 5.14, 10542, 1],
    [1, 9.0, 3.94, 17974, 1],
    [0, 3.0, 4.51, 11986, 2],
    [1, 7.18, 3.71, 17237, 2],
    [0, 13.15, 4.31, 24494, 3],
    [1, 7.5, 3.8, 18129, 3]
])

# Extract the unique theta and battery levels
theta_levels = np.unique(data[:, 0])
battery_levels = np.unique(data[:, 1])

# Define the bar chart parameters
bar_width = 0.35
space_width = 0.1

# Get the number of policies and calculate the total bar width
num_policies = len(np.unique(data[:, 4]))
total_bar_width = num_policies * len(battery_levels) * (bar_width + space_width)

data_plot = []

# Iterate over the battery levels and theta levels to create the bars
for i, theta_level in enumerate(theta_levels):
    for j, battery_level in enumerate(battery_levels):
        for k, policy in enumerate(np.unique(data[:, 4])):
            values = data[(data[:, 0] == theta_level) & (data[:, 1] == battery_level) & (data[:, 4] == policy), 2]
            if (not np.any(values)): continue
            if policy == 1:
                if theta_level == 0.0:
                    local_data = ["WFH T1", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T1", "bidirectional", values[0]]
            elif policy == 2:
                if theta_level == 0.0:
                    local_data = ["WFH T2", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T2", "bidirectional", values[0]]
            elif policy == 3:
                if theta_level == 0.0:
                    local_data = ["WFH T3", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T3", "bidirectional", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'PV (kW)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['PV (kW)'] = df['PV (kW)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="PV (kW)", hue="Connectivity", hue_order=["unidirectional","bidirectional"], order=["WFH T1", "WFH T2", "WFH T3"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 7)
sns_plot.set_xlabel('')



plt.show()













import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0, 2, 17974, 10542, 1],
    [1, 9.0, 11986, 17974, 1],
    [0, 3.0, 14282, 11986, 2],
    [1, 7.18, 10290, 17237, 2],
    [0, 13.15, 13353, 24494, 3],
    [1, 7.5,10291 , 18129, 3]
])

# Extract the unique theta and battery levels
theta_levels = np.unique(data[:, 0])
battery_levels = np.unique(data[:, 1])

# Define the bar chart parameters
bar_width = 0.35
space_width = 0.1

# Get the number of policies and calculate the total bar width
num_policies = len(np.unique(data[:, 4]))
total_bar_width = num_policies * len(battery_levels) * (bar_width + space_width)

data_plot = []

# Iterate over the battery levels and theta levels to create the bars
for i, theta_level in enumerate(theta_levels):
    for j, battery_level in enumerate(battery_levels):
        for k, policy in enumerate(np.unique(data[:, 4])):
            values = data[(data[:, 0] == theta_level) & (data[:, 1] == battery_level) & (data[:, 4] == policy), 2]
            if (not np.any(values)): continue
            if policy == 1:
                if theta_level == 0.0:
                    local_data = ["WFH T1", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T1", "bidirectional", values[0]]
            elif policy == 2:
                if theta_level == 0.0:
                    local_data = ["WFH T2", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T2", "bidirectional", values[0]]
            elif policy == 3:
                if theta_level == 0.0:
                    local_data = ["WFH T3", "unidirectional", values[0]]
                elif theta_level == 1.0:
                    local_data = ["WFH T3", "bidirectional", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'Cost ($)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['Cost ($)'] = df['Cost ($)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Cost ($)", hue="Connectivity", hue_order=["unidirectional","bidirectional"], order=["WFH T1", "WFH T2", "WFH T3"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 23000)
sns_plot.set_xlabel('')



plt.show()











