## Script used to generate scatter alu main figure.
## Usage: python plot_scatteralu.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Data
data = {
    'SC': ['MSA1 T7', 'MSA1 RBP', 'Control RBP', 'Control T7', 'MSA2 RBP', 'MSA2 T7'],
    'SC_AluJ': [5, 21, 2, 10, 139, 149],
    'SC_AluS': [25, 40, 6, 32, 283, 311],
    'SC_AluY': [3, 20, 4, 10, 67, 91],
    'SCBulk_AluJ': [1, 31, 20, 38, 23, 76],
    'SCBulk_AluS': [2, 130, 45, 103, 125, 260],
    'SCBulk_AluY': [3, 395, 161, 477, 240, 635],
    'BulkOnly_AluJ': [0, 61, 21, 60, 57, 110],
    'BulkOnly_AluS': [3, 194, 97, 171, 232, 291],
    'BulkOnly_AluY': [1, 156, 92, 167, 153, 236]
}
plt.rcParams.update({'font.size': 17})
# Create a DataFrame
df = pd.DataFrame(data)

# Calculate the average for each Alu type across experiments
averages = {
    'SCBulk_AluJ': df['SCBulk_AluJ'].mean(),
    'SCBulk_AluS': df['SCBulk_AluS'].mean(),
    'SCBulk_AluY': df['SCBulk_AluY'].mean(),
    'SC_AluJ': df['SC_AluJ'].mean(),
    'SC_AluS': df['SC_AluS'].mean(),
    'SC_AluY': df['SC_AluY'].mean(),
    'BulkOnly_AluJ': df['BulkOnly_AluJ'].mean(),
    'BulkOnly_AluS': df['BulkOnly_AluS'].mean(),
    'BulkOnly_AluY': df['BulkOnly_AluY'].mean()
}

# Calculate Total bulk (single cell regions) = SCBulk + BulkOnly for each Alu type
total_bulk_averages = {
    'TotalBulk_AluJ': averages['SCBulk_AluJ'] + averages['BulkOnly_AluJ'],
    'TotalBulk_AluS': averages['SCBulk_AluS'] + averages['BulkOnly_AluS'],
    'TotalBulk_AluY': averages['SCBulk_AluY'] + averages['BulkOnly_AluY']
}

# Create a new DataFrame for the averages, including Total bulk
df_avg = pd.DataFrame({
    'Group': ['TotalBulk', 'TotalBulk', 'TotalBulk', 'SCBulk', 'SCBulk', 'SCBulk',
              'SC', 'SC', 'SC', 'BulkOnly', 'BulkOnly', 'BulkOnly'],
    'AluType': ['AluJ', 'AluS', 'AluY', 'AluJ', 'AluS', 'AluY', 
                'AluJ', 'AluS', 'AluY', 'AluJ', 'AluS', 'AluY'],
    'Average': [total_bulk_averages['TotalBulk_AluJ'], total_bulk_averages['TotalBulk_AluS'], total_bulk_averages['TotalBulk_AluY'],
                averages['SCBulk_AluJ'], averages['SCBulk_AluS'], averages['SCBulk_AluY'],
                averages['SC_AluJ'], averages['SC_AluS'], averages['SC_AluY'],
                averages['BulkOnly_AluJ'], averages['BulkOnly_AluS'], averages['BulkOnly_AluY']]
})

# Mapping for category positions on the x-axis
x_map = {
    'TotalBulk': 0, 'SCBulk': 1, 'SC': 2, 'BulkOnly': 3
}

# Create a new column to assign numerical values to x-axis positions
df_avg['X_numeric'] = df_avg['Group'].map(x_map)

# Set up color and marker mappings for Alu types
alu_colors = {'AluJ': 'blue', 'AluS': 'green', 'AluY': 'red'}
alu_markers = {'AluJ': 'o', 'AluS': 's', 'AluY': 'D'}

# Create the plot
plt.figure(figsize=(8, 6))

# Use seaborn scatterplot to differentiate both by color (hue) and marker (style)
sns.scatterplot(x='X_numeric', y='Average', hue='AluType', style='AluType', 
                palette=alu_colors, markers=alu_markers, data=df_avg, s=150, zorder=10)

# Manually set the x-ticks to correspond to the Group
plt.xticks(ticks=[0, 1, 2, 3], labels=['Total bulk \n(single cell regions)', 'Single cell + \nbulk', 'Single cell-only', 'Bulk-only'], rotation=0,fontweight='bold',fontsize=13)
plt.grid(which='major', axis='y', linestyle='-', linewidth=0.75, color='lightgray')

# Add labels and title with adjusted x label position
plt.title('Average Alu insertion counts \n across experiments', weight='bold',fontsize=24)
plt.xlabel('Sample Group', weight='bold', labelpad=10,fontsize=22)  # Adjust labelpad to lower x-axis label
plt.ylabel('Average Variant Counts', weight='bold',fontsize=22, labelpad=15)

# Customize the legend to be inside the scatter area
plt.legend(title='Alu Type', loc='upper right',fontsize=17)
plt.yticks(fontsize=17) 
# Save the plot as an image file
plt.savefig('alu_insertions_avg_plot_colored_symbols.png', dpi=1200, bbox_inches='tight')

# Show the plot
plt.tight_layout()
plt.show()
