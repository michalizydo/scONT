## Script used to generate bar chart coverage - main figure.
## Usage: python bar_coverage.py

import matplotlib.pyplot as plt

# Use the 'seaborn-white' style for plots
plt.style.use('seaborn-white')

# Data
categories = [
    "Control - single cell ONT",
    "Control - single cell Illumina",
    "Control - bulk ONT",
    "MSA1 - single cell ONT",
    "MSA1 - single cell Illumina",
    "MSA1 - bulk ONT",
    "MSA2 - single cell ONT",
    "MSA2  - single cell Illumina",
    "MSA2 - bulk ONT"
]
values = [21.73, 47.39, 93.65, 35.44, 62.46, 94.5, 45.72, 61.1, 94.02]

# Identify Illumina and bulk samples
illumina_samples = ["Control - single cell Illumina", "MSA1 - single cell Illumina", "MSA2  - single cell Illumina"]
bulk_samples = ["Control - bulk ONT", "MSA1 - bulk ONT", "MSA2 - bulk ONT"]

# Create figure and axis with specific figure size
fig, ax = plt.subplots(figsize=(13, 6))

# Define colors for bars: green for Illumina samples, dark orange for bulk, blue-gray for others
colors = [
    '#4A70B0' if category not in illumina_samples + bulk_samples else (
        '#A85DA2' if category in illumina_samples else '#D2691E'
    ) for category in categories
]

# Create the bars without borders
bars = ax.bar(categories, values, color=colors)

# Add labels and title with enhanced font properties
ax.set_xlabel('Brain samples', fontsize=22, fontweight='bold', labelpad=65)
ax.set_ylabel('% of genome covered', fontsize=22, fontweight='bold', labelpad=0)
ax.set_title('Breadth of coverage across 3 tested brain samples', fontsize=24, fontweight='bold', pad=15)

# Add custom labels below each bar
custom_labels = ["Single cell\nONT", "Single cell\nIllumina", "Bulk ONT", "Single cell\nONT", "Single cell\nIllumina MDA", "Bulk ONT", "Single cell\nONT", "Single cell\nIllumina MDA", "Bulk ONT"]
for bar, label in zip(bars, custom_labels):
    # Set a higher y-position for blue and orange bar labels based on keywords
    if "Bulk ONT" in label:  # Assuming "ONT" bars are blue and orange
        label_y = -7.5  # Slightly higher position for these bars
    else:
        label_y = -11.5  # Original position for other bars

    ax.text(bar.get_x() + bar.get_width() / 2, label_y, label, ha='center', va='bottom', fontsize=15, fontweight='bold')

ax.ticklabel_format(axis='y')

# Customize ticks and their appearance
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xticks([])
ax.grid(True, which='major', axis='y', linestyle='-', linewidth=0.75, color='lightgray')
ax.set_axisbelow(True)

# Remove top and right spines for a cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add group labels for "Control", "MSA1", and "MSA2"
group_labels = ['Control', 'MSA1', 'MSA2']
group_positions = [1, 4, 7]  # The center of each group of 3 bars

for label, pos in zip(group_labels, group_positions):
    ax.text(pos, -0.15, label, ha='center', va='top', fontsize=18, fontweight='bold', transform=ax.get_xaxis_transform(), color='black')

# Set y-axis limits and adjust layout
ax.set_ylim(0, 100)
plt.tight_layout()

# Save the figure in high resolution
plt.savefig('barchart_illumina_violet_bulk_orange.png', dpi=600)

# Show the plot
plt.show()
