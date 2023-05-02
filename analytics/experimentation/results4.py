# UNIDIRECTIONAL PART
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0.60, 2.86, 18.06, 10542, 0],
    [0.75, 9.0, 11.30, 17974, 0],
    [0.90, 9.0,4.75 , 17974, 0],
    [0.60, 3.0, 12.68, 11986, 1],
    [0.75, 7.18, 6.40, 17237, 1],
    [0.90, 9.0, 1.25, 17974, 1],
    [0.60, 13.15, 12.81, 24494, 2],
    [0.75, 7.5, 5.90, 18129, 2],
    [0.90, 9.0, 0.95, 17974, 2]
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
                local_data = ["WFH T1", theta_level, values[0]]
            elif policy == 1:
                local_data = ["WFH T2", theta_level, values[0]]
            elif policy == 2:
                local_data = ["WFH T3", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'Battery (kWh)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['Battery (kWh)'] = df['Battery (kWh)'].astype(np.float32)

plt.figure(figsize=(8, 7))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Battery (kWh)", hue="Theta", hue_order=theta_levels, order=["WFH T1", "WFH T2", "WFH T3"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 20)
sns_plot.set_xlabel('')



plt.show()










import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0.60, 2.86, 8.60, 10542, 0],
    [0.75, 9.0, 5.94, 17974, 0],
    [0.90, 9.0, 3.62, 17974, 0],
    [0.60, 3.0, 7.88, 11986, 1],
    [0.75, 7.18, 5.17, 17237, 1],
    [0.90, 9.0, 2.65, 17974, 1],
    [0.60, 13.15, 7.74, 24494, 2],
    [0.75, 7.5, 5.14, 18129, 2],
    [0.90, 9.0, 2.60, 17974, 2]
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
                local_data = ["WFH T1", theta_level, values[0]]
            elif policy == 1:
                local_data = ["WFH T2", theta_level, values[0]]
            elif policy == 2:
                local_data = ["WFH T3", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'PV (kW)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['PV (kW)'] = df['PV (kW)'].astype(np.float32)

plt.figure(figsize=(8, 7))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="PV (kW)", hue="Theta", hue_order=theta_levels, order=["WFH T1", "WFH T2", "WFH T3"])

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
    [0.60, 2.86, 31442, 10542, 0],
    [0.75, 9.0, 21184, 17974, 0],
    [0.90, 9.0, 11945, 17974, 0],
    [0.60, 3.0, 27048, 11986, 1],
    [0.75, 7.18, 16855, 17237, 1],
    [0.90, 9.0, 7722, 17974, 1],
    [0.60, 13.15, 26722, 24494, 2],
    [0.75, 7.5, 16543, 18129, 2],
    [0.90, 9.0, 7431, 17974, 2]
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
                local_data = ["WFH T1", theta_level, values[0]]
            elif policy == 1:
                local_data = ["WFH T2", theta_level, values[0]]
            elif policy == 2:
                local_data = ["WFH T3", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8", "#b284be"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'Cost ($)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['Cost ($)'] = df['Cost ($)'].astype(np.float32)

plt.figure(figsize=(8, 7))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Cost ($)", hue="Theta", hue_order=theta_levels, order=["WFH T1", "WFH T2", "WFH T3"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 33000)
sns_plot.set_xlabel('')



plt.show()










# BIDIRECTIONAL PART 