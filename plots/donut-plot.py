import matplotlib.pyplot as plt
import numpy as np


# Data1 for the donut plot
outcome = ["WGTS reporting", "Failed sample QC (not sequenced)", "Failed WGS QC (sequenced)"]
category = ["Successful", "Unsuccessful", "Unsuccessful"]
values = [45, 5, 4]

# Colors
colors = ["lightgreen", "grey", "lightblue"]

# Create a pie chart
fig, ax = plt.subplots(figsize=(6, 6))
wedges, texts, autotexts = ax.pie(values, colors=colors, labels=outcome, autopct='%1.1f%%', startangle=90,
pctdistance=0.85, labeldistance=1.1, textprops={'fontsize': 10})

# Set aspect ratio to be equal so that the pie is drawn as a circle
ax.axis('equal')

# Create a circle in the center
centre_circle = plt.Circle((0, 0), 0.7, fc='white')
fig.gca().add_artist(centre_circle)
ax.set_title("Samples Overview", fontsize=14)

# Show the chart
# plt.show()

# Data2 for the pie plot
outcome = ["WGTS reporting", "Failed sample QC (not sequenced)", "Failed WGS QC (sequenced)"]
values = [45, 5, 4]

# Define colors
colors = ["forestgreen", "cornflowerblue", "orange"]

# Create a pie chart
fig, ax = plt.subplots(figsize=(6, 6))
ax.pie(values, labels=outcome, autopct='%1.1f%%', startangle=90, colors=colors)
ax.set_title("Samples Overview", fontsize=14)

# Show the chart
plt.show()