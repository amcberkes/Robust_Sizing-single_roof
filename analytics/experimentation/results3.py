import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [1, 2, 9.00, 10542, 0],
    [1, 9.0, 4.00, 17974, 1],
    [1, 2, 3.50, 10542, 2],
    [1, 3.0, 3.00, 11986, 3],
    [2, 7.18, 4.10, 17237, 0],
    [2, 13.15, 0.55, 15171, 1],
    [2, 7.5,0.70 , 11134, 2],
    [2, 2,0.55 , 10542, 3],
    [3, 2, 4.05, 10542, 0],
    [3, 2, 0.15, 10542, 1],
    [3, 2, 0.25, 10542, 2],
    [3, 2, 0.15, 10542, 3]
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
            if policy == 0:
                if theta_level == 1:
                    local_data = ["unidirectional", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["unidirectional", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["unidirectional", "WFH T3", values[0]]
            if policy == 1:
                if theta_level == 1:
                    local_data = ["min_storage", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["min_storage", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["min_storage", "WFH T3", values[0]]
            elif policy == 2:
                if theta_level == 1:
                    local_data = ["r_degradation", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["r_degradation", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["r_degradation", "WFH T3", values[0]]
            elif policy == 3:
                if theta_level == 1:
                    local_data = ["most_sustainable", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["most_sustainable", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["most_sustainable", "WFH T3", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'Battery (kWh)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['Battery (kWh)'] = df['Battery (kWh)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Battery (kWh)", hue="Connectivity", hue_order=["WFH T1","WFH T2", "WFH T3"], order=["unidirectional", "min_storage", "r_degradation", "most_sustainable"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 11)
sns_plot.set_xlabel('')



plt.show()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [1, 2, 5.14, 10542, 0],
    [1, 9.0, 3.85, 17974, 1],
    [1, 2, 3.85, 10542, 2],
    [1, 3.0, 3.94, 11986, 3],
    [2, 7.18, 4.45, 17237, 0],
    [2, 13.15, 3.71, 15171, 1],
    [2, 7.5, 3.68, 11134, 2],
    [2, 2, 3.68, 10542, 3],
    [3, 2,4.20 , 10542, 0],
    [3, 2, 3.80, 10542, 1],
    [3, 2, 3.77, 10542, 2],
    [3, 2, 3.77, 10542, 3]
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
            if policy == 0:
                if theta_level == 1:
                    local_data = ["unidirectional", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["unidirectional", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["unidirectional", "WFH T3", values[0]]
            if policy == 1:
                if theta_level == 1:
                    local_data = ["min_storage", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["min_storage", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["min_storage", "WFH T3", values[0]]
            elif policy == 2:
                if theta_level == 1:
                    local_data = ["r_degradation", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["r_degradation", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["r_degradation", "WFH T3", values[0]]
            elif policy == 3:
                if theta_level == 1:
                    local_data = ["most_sustainable", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["most_sustainable", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["most_sustainable", "WFH T3", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'PV (kW)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['PV (kW)'] = df['PV (kW)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="PV (kW)", hue="Connectivity", hue_order=["WFH T1","WFH T2", "WFH T3"], order=["unidirectional", "min_storage", "r_degradation", "most_sustainable"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 6)
sns_plot.set_xlabel('')



plt.show()




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [1, 2, 17974, 10542, 0],
    [1, 9.0, 12215, 17974, 1],
    [1, 2, 11985, 10542, 2],
    [1, 3.0, 11986, 11986, 3],
    [2, 7.18, 13875, 17237, 0],
    [2, 13.15, 10244, 15171, 1],
    [2, 7.5, 10236, 11134, 2],
    [2, 2, 10167, 10542, 3],
    [3, 2,13161 , 10542, 0],
    [3, 2, 10291, 10542, 1],
    [3, 2, 10260, 10542, 2],
    [3, 2, 10214, 10542, 3]
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
            if policy == 0:
                if theta_level == 1:
                    local_data = ["unidirectional", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["unidirectional", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["unidirectional", "WFH T3", values[0]]
            if policy == 1:
                if theta_level == 1:
                    local_data = ["min_storage", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["min_storage", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["min_storage", "WFH T3", values[0]]
            elif policy == 2:
                if theta_level == 1:
                    local_data = ["r_degradation", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["r_degradation", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["r_degradation", "WFH T3", values[0]]
            elif policy == 3:
                if theta_level == 1:
                    local_data = ["most_sustainable", "WFH T1", values[0]]
                elif theta_level == 2:
                    local_data = ["most_sustainable", "WFH T2", values[0]]
                elif theta_level == 3:
                    local_data = ["most_sustainable", "WFH T3", values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas


colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Connectivity', 'Cost ($)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Connectivity'] = df['Connectivity'].astype('category')
df['Cost ($)'] = df['Cost ($)'].astype(np.float32)


plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Cost ($)", hue="Connectivity", hue_order=["WFH T1","WFH T2", "WFH T3"], order=["unidirectional", "min_storage", "r_degradation", "most_sustainable"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 19000)
sns_plot.set_xlabel('')



plt.show()




