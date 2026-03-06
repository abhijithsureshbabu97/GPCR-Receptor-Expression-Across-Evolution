import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Read the data from Excel file
data = pd.read_excel('/home/asb55/irp/results_file/Droso/Dorso_GPCRfiltered.xlsx')

# Extract the feature names and IDs
feature_names = data['features.plot']
ids = data['id']
OG = data['OG']
# Get unique ortho groups and assign colors to each group
ortho_groups = data['OG'].unique()
num_groups = len(ortho_groups)
colors = sns.color_palette('Set1', num_groups)  # Use a color palette from seaborn

# Create a dictionary to map ortho groups to colors
group_colors = dict(zip(ortho_groups, colors))

# Map the ortho groups to colors for each gene
gene_colors = [group_colors[og] for og in OG]

x_labels = [f'{name} ({og})' for name, og in zip(feature_names, OG)]


# Create the scatter plot with grouped colors
fig=plt.figure(figsize=(50, 20))
plt.scatter(x_labels, ids, c=gene_colors, s=100)
plt.xlabel('Feature Names + OG')
plt.ylabel('ID')
plt.title('Feature Names vs. ID')
plt.xticks(rotation=45, ha='right')
plt.yticks(ids)




plt.grid()
# Create legend
legend_elements = [plt.Line2D([30], [30], marker='o', color='w', label=group, markerfacecolor=color, markersize=10)
                   for group, color in zip(ortho_groups, colors)]
plt.legend(handles=legend_elements, title='Ortho Groups')
plt.savefig('scatter_plot.jpeg')
plt.show()






# Create the joint plot
joint_plot = sns.jointplot(x=x_labels, y=ids, s=100, marginal_kws=dict(bins=67, fill=True))
joint_plot.fig.set_figwidth(90)  # Adjust the width of the figure 
joint_plot.fig.set_figheight(40)  # Adjust the width of the figure 

plt.ylim(-1, 70)
plt.xlim(-1,None)
plt.grid()

plt.yticks(ids)
plt.xticks(rotation=45, ha='right')

plt.savefig('GPCR_DROSO.jpeg')
plt.show()

