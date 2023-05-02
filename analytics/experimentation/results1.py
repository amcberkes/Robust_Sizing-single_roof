# impact of ev connectivity type

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0.8, 2.86, 2.86, 10542, 0],
    [0.8, 9.0, 9.0, 17974, 1],
    [0.8, 3.0, 3.0, 11986, 2],
    [0.7, 7.18, 7.18, 17237, 0],
    [0.7, 13.15, 13.15, 24494, 1],
    [0.7, 7.5, 7.5, 18129, 2]
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
    for j, battery_level in enumerate([2.86, 9.0, 3.0, 7.18, 13.15, 7.5]):
        for k, policy in enumerate(np.unique(data[:, 4])):
            values = data[(data[:, 0] == theta_level) & (data[:, 1] == battery_level) & (data[:, 4] == policy), 2]
            #values = [2.86, 9.0, 3.0, 7.18, 13.15, 7.5]
            if (not np.any(values)): continue
            if policy == 0:
                local_data = ["no ev", theta_level, values[0]]
            elif policy == 1:
                local_data = ["unidirectional", theta_level, values[0]]
            elif policy == 2:
                local_data = ["bidirectional", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'Battery (kWh)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['Battery (kWh)'] = df['Battery (kWh)'].astype(np.float32)

plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Battery (kWh)", hue="Theta", hue_order=theta_levels, order=["no ev", "unidirectional", "bidirectional"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 14)
sns_plot.set_xlabel('')

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0.8, 2.86, 3.43, 10542, 0],
    [0.8, 9.0, 5.14, 17974, 1],
    [0.8, 3.0, 3.94, 11986, 2],
    [0.7, 7.18, 5.18, 17237, 0],
    [0.7, 13.15, 6.85, 24494, 1],
    [0.7, 7.5, 5.45, 18129, 2]
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
                local_data = ["no ev", theta_level, values[0]]
            elif policy == 1:
                local_data = ["unidirectional", theta_level, values[0]]
            elif policy == 2:
                local_data = ["bidirectional", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'PV (kW)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['PV (kW)'] = df['PV (kW)'].astype(np.float32)

plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="PV (kW)", hue="Theta", hue_order=theta_levels, order=["no ev", "unidirectional", "bidirectional"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 8)
sns_plot.set_xlabel('')



plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define the data
data = np.array([
    [0.8, 2.86, 10542, 10542, 0],
    [0.8, 9.0, 17974, 17974, 1],
    [0.8, 3.0, 11986, 11986, 2],
    [0.7, 7.18, 17237, 17237, 0],
    [0.7, 13.15, 24494, 24494, 1],
    [0.7, 7.5, 18129, 18129, 2]
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
                local_data = ["no ev", theta_level, values[0]]
            elif policy == 1:
                local_data = ["unidirectional", theta_level, values[0]]
            elif policy == 2:
                local_data = ["bidirectional", theta_level, values[0]]
            data_plot.append(local_data)

import seaborn as sns
import pandas

colors = ["#F8C8DC", "#D8C5E8"]
sns.set_palette(sns.color_palette(colors))

df = pandas.DataFrame(data_plot, columns=['EV operation mode', 'Theta', 'Cost ($)'])
df['EV operation mode'] = df['EV operation mode'].astype('category')
df['Theta'] = df['Theta'].astype('category')
df['Cost ($)'] = df['Cost ($)'].astype(np.float32)

plt.figure(figsize=(6, 5))
sns_plot = sns.barplot(data=df, x="EV operation mode", y="Cost ($)", hue="Theta", hue_order=theta_levels, order=["no ev", "unidirectional", "bidirectional"])

for p in sns_plot.containers:
    sns_plot.bar_label(p, label_type='edge', labels=np.round(p.datavalues,2), padding=2)

sns_plot.set_ylim(0, 26000)
sns_plot.set_xlabel('')



plt.show()









