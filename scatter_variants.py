## Script used to generate scatter of variant counts across experiments - main figure.
## Usage: python scatter_variants.py

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

# Create the data
data = {
    'Category': [
        'MSA1 - All bulk',
        'MSA1 - Total Bulk single-cell regions',
        'MSA1 - Single cell only',
        'MSA1 - Single cell bulk',
        'MSA1 - Bulk only',
        'MSA2 - All bulk',
        'MSA2 - Total Bulk single-cell regions',
        'MSA2 - Single cell only',
        'MSA2 - Single cell bulk',
        'MSA2 - Bulk only',
        'Control - All bulk',
        'Control - Total Bulk single-cell regions',
        'Control - Single cell only',
        'Control - Single cell bulk',
        'Control - Bulk only'
    ],
    'Values': [
        21665, 4784, 2270, 2946, 1412,
        21831, 10868, 4722, 8466, 3025,
        22040, 5336, 3414, 3571, 1922
    ]
}

# Convert to DataFrame
df = pd.DataFrame(data)

# Create figure and axis with specific figure size
fig, ax = plt.subplots(figsize=(19, 14))

# Define category labels and corresponding positions
categories = ['All bulk', 'Total Bulk single-cell regions', 'Single cell only', 'Single cell bulk', 'Bulk only']
positions = [1, 1.95, 3.05, 3.9, 4.7]

# Define colors for each point in the specified order
colors = ['green', 'red', '#2E4A8B']  # Green, Red, Blue

# Scatter plot with random jitter
np.random.seed(0)  # For reproducibility
for i, category in enumerate(categories):
    # Get values for the corresponding category
    values = df['Values'][df['Category'].str.contains(category)].values
    
    # Loop through each value, assigning colors in order
    for j, value in enumerate(values):
        jittered_x = np.random.normal(loc=positions[i], scale=0.05)  # adding jitter to the x-coordinate
        color = colors[j % len(colors)]  # Alternate colors for each value
        ax.scatter(jittered_x, value, color=color, alpha=0.7, s=600)  # Increased size

# Customize appearance
ax.set_xlabel('Sample group', fontsize=56, fontweight='bold', labelpad=22)
ax.set_ylabel('Variant counts', fontsize=56, fontweight='bold', labelpad=10)
ax.set_title('Variants (INS + DEL) \ndetected in examined brains', fontsize=62, fontweight='bold', pad=50)

# Customize ticks and their appearance
ax.tick_params(axis='x', labelsize=40)
ax.tick_params(axis='y', labelsize=40)  # Increase y-axis tick label size
# Add grid lines behind the boxes
ax.grid(True, which='major', axis='y', linestyle='-', linewidth=0.75, color='lightgray')
ax.set_axisbelow(True)  # Ensures grid is behind the boxes
ax.set_yticklabels([f"{int(tick) // 1000}k" for tick in ax.get_yticks()], fontsize=40)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Setting xticks
ax.set_xticks(positions)
ax.set_xticklabels([
    'All bulk',
    '  Total Bulk\n(single-cell regions)',
    'Single cell +\n bulk',
    'Single cell\n only',
    'Bulk only'
], fontsize=32, fontweight='bold', verticalalignment='top')

# Move the y-axis to the left to reduce overlap
ax.spines['left'].set_position(('outward', 20))

# Adjust the y-axis limits dynamically based on the data
ax.set_ylim(0, 25000)

# Custom legend entries for TE type and group (Insertions and Deletions)
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=36, label='MSA1'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=36, label='MSA2'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2E4A8B', markersize=36, label='Control'),
]

# Add the custom legend to the plot
plt.legend(handles=legend_elements, title='Brain sample', fontsize=32, title_fontsize=36, handletextpad=0.1)

# Adjust layout for better fit
plt.tight_layout()
plt.savefig('scatter_variants_colored.png', dpi=600)

# Show the plot
plt.show()