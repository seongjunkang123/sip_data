import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler

# remember to change the csv file name
df = pd.read_csv("/Users/user/COPD_peaktable_ver2.csv")
df = df.set_index('pubchem_CID')

scaler = StandardScaler()
scaled_df = pd.DataFrame(scaler.fit_transform(df.T), columns=df.index, index=df.columns)

plt.figure(figsize=(9, 4))

# Adjust colorbar width (cbar_kws) and color scale (cbar_kws)
sns.clustermap(
    scaled_df,
    metric="correlation",
    method="single",
    cmap="Blues",
    standard_scale=1,
    cbar_kws={"shrink": 0.5},  # Adjust the width of the colorbar
    dendrogram_ratio=(0.1, 0.2),  # Adjust the height of the dendrograms
    cbar_pos=(0.2, 0.2, 0.02, 0.4)  # Adjust the position of the colorbar
)

# Adjust subplot parameters to move the plot title away from the top
plt.subplots_adjust(top=0.95)  # You can adjust the top value as needed

# Add a centered, larger, and bold title
plt.suptitle('COPD heat map', fontsize=16, fontweight='bold')

plt.show()
