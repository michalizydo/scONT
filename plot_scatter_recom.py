## Script used to generate scatter of TE recombining - main figure.
## Usage: python plot_scatter_recom.py

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

plt.rcParams.update({'font.size': 12})

# Data for Repeat Recombinants and Insertions into same TE
data = {
    'Category': [
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only",
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only",
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only",
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only"
    ],
    'Value': [
        213, 136, 78, 77, 195, 115, 39, 80,  # Insertions
        92, 56, 73, 36, 111, 76, 62, 35      # Deletions
    ],
    'Group': [
        'Insertions', 'Insertions', 'Insertions', 'Insertions',
        'Insertions', 'Insertions', 'Insertions', 'Insertions',
        'Deletions', 'Deletions', 'Deletions', 'Deletions',
        'Deletions', 'Deletions', 'Deletions', 'Deletions'
    ],
    'Type': [
        'Alu-Alu', 'Alu-Alu', 'Alu-Alu', 'Alu-Alu',
        'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1',
        'Alu-Alu', 'Alu-Alu', 'Alu-Alu', 'Alu-Alu',
        'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1', 'LINE/L1-LINE/L1'
    ]
}

# Convert the data into a DataFrame
df = pd.DataFrame(data)

# Create a new column 'x_offset' to adjust x positions
df['x_offset'] = np.tile([0, 1, 2, 3], 4)

# Adjust x_offset for Deletions and Insertions
df.loc[df['Group'] == 'Deletions', 'x_offset'] = np.tile([0, 1, 2, 3], 2)
df.loc[df['Group'] == 'Insertions', 'x_offset'] = np.tile([4.5, 5.5, 6.5, 7.5], 2)

# Shift Alu points slightly to the left and LINE points to the right
df.loc[df['Type'] == 'Alu-Alu', 'x_offset'] -= 0.1
df.loc[df['Type'] == 'LINE/L1-LINE/L1', 'x_offset'] += 0.1

# Set up color mappings for types
type_colors = {'Alu-Alu': 'red', 'LINE/L1-LINE/L1': 'green'}

# Create the plot
plt.figure(figsize=(9, 7))

# Plot the data using Seaborn scatterplot with adjusted x_offset
sns.scatterplot(
    x='x_offset', y='Value', 
    hue='Type',  
    palette=type_colors, 
    s=150, data=df[df['Group'] == 'Deletions'],
    marker='o', zorder=2
)

sns.scatterplot(
    x='x_offset', y='Value', 
    hue='Type',  
    palette=type_colors, 
    s=150, data=df[df['Group'] == 'Insertions'],
    marker='D', zorder=2
)

# Set x-axis tick positions and labels
plt.xticks(
    ticks=[0, 1, 2, 3, 4.5, 5.5, 6.5, 7.5],
    labels=[
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only",
        "Total bulk\n(single cell regions)", "Single cell+bulk", "Single cell-only", "Bulk-only"
    ],
    ha='center', fontsize=15
)

# Remove the x_offset axis label
plt.xlabel('')

# Rotate the x-axis labels if necessary
plt.xticks(rotation=45, ha='right', fontweight='bold')

# Add Y-axis label and make it bold
plt.ylabel('# of events', fontsize=22, fontweight='bold')

# Add plot title
plt.title('Comparison of TE sequence events\n across brain samples', fontsize=24, fontweight='bold')
plt.grid(True, which='major', axis='y', linestyle='-', linewidth=0.75, color='lightgray', zorder=0)

# Add a vertical dotted line to divide Deletions and Insertions sections
plt.axvline(x=3.7, color='gray', linestyle='--', linewidth=1, zorder=1)

# Custom legend entries for TE type and group (Insertions and Deletions)
legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='red', markersize=12, label='Alu to Alu Insertions'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='green', markersize=12, label='LINE/L1 to LINE/L1 Insertions'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=12, label='Alu-Alu Deletions'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=12, label='LINE/L1-LINE/L1 Deletions')
]

# Add the custom legend to the plot
plt.legend(handles=legend_elements, title='Event Type', fontsize=14, title_fontsize=16, handletextpad=0.1)

# Adjust plot layout
plt.tight_layout(pad=1.0)
plt.yticks(fontsize=17)

# Optionally save the figure
plt.savefig('insertions_deletions_custom_legend_with_divider.png', dpi=600)

# Show the plot
plt.show()
