from matplotlib.lines import lineStyles
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from moepy import lowess
import seaborn as sns 
import matplotlib.font_manager as fm
from datetime import datetime
from sklearn.metrics import roc_curve
from sklearn import metrics
from sklearn.metrics import PrecisionRecallDisplay, precision_recall_curve, average_precision_score
from scipy import stats
import numpy as np
import math
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
import os
from statannot import add_stat_annotation
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from statsmodels.stats.multitest import multipletests




df = pd.read_csv("Result/Coverage_SNV.csv")

# Filter data to exclude the "Total" row
df = df[df['Exon'] != 'Total']

# Convert necessary columns to numeric
df['SNV_Percentage'] = pd.to_numeric(df['SNV_Percentage'])
df['induced'] = pd.to_numeric(df['induced'])
df['total(SNV)'] = pd.to_numeric(df['total(SNV)'])

# Define colors
colors = ['#123753', '#955196']

# Create the bar plot
fig, ax = plt.subplots(figsize=(7, 6))
bars = ax.bar(
    df['Exon'], 
    df['SNV_Percentage'], 
    color=colors[0], 
    edgecolor=None, 
    linewidth=1.2,
    width=0.6
)

# Customize tick parameters
plt.tick_params(labelsize=20, size=10, width=1, colors='black')


# Customize the plot
ax.set_ylabel("SNV Coverage (%)", fontsize=12)
ax.set_xlabel("Exon", fontsize=12)
ax.set_ylim(0, 100)  # Set y-limit slightly above 100 for labels
ax.set_xticks(df['Exon'])
ax.set_xticklabels(df['Exon'], fontsize=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.spines['left'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')


# Show the plot
plt.tight_layout()
plt.show()


fig.savefig('Graph/Barplot/Coverage_SNV.pdf',dpi=300,bbox_inches='tight', transparent = True)



### Scatterplot ###

mutation_types = ['Syn', 'Missense', 'Nonsense', 'Splice', 'Intron', 'Duplication', 'Deletion'] 
color_dic = {'Syn': '#3498DB', 'Missense': '#2ECC71', 'Nonsense': '#E74C3C',  'Splice': '#F39C12', 'Intron' : 'gray', 'Deletion': '#8E44AD', 'Duplication' : 'black'}



# %%
## Drug Scatterplot ##
drug = 'TPX'

fig, axe = plt.subplots(figsize=(8,6))

for mutation_type in ['Missense', 'Nonsense', 'Splice', 'Intron','Deletion','Duplication']:
    eachdf = df[df['Type'] == mutation_type]
    sns.scatterplot(data=eachdf, x='Resistance_Score1', y='Resistance_Score2', label=mutation_type, 
                    s=20, edgecolor='black', linewidth=0.6, alpha=0.7, color=color_dic.get(mutation_type))

eachdf_syn = df[df['Type'] == 'Syn']
sns.scatterplot(data=eachdf_syn, x='Resistance_Score1', y='Resistance_Score2', label='Syn', 
                s=20, edgecolor='black', linewidth=0.6, alpha=0.9, color=color_dic.get('Syn'))


axe.spines['left'].set_color('black')
axe.spines['bottom'].set_color('black')

# Set tick parameters
plt.tick_params(labelsize=20, size=10, width=1)
plt.xticks([-5, 0, 5, 10])
plt.yticks([-5, 0, 5, 10])
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize=14, markerscale=1.5, frameon=False)
axe.grid(False)


r, p = stats.pearsonr(df['Resistance_Score1'], df['Resistance_Score2'])

plt.title('Pearson r = {:.2f}, p = {:.2e} (n={})'.format(r, p, df.shape[0]), fontsize=10, weight='bold')

# Axis labels
plt.xlabel('Resistance_Score1', fontsize=16, weight='bold')
plt.ylabel('Resistance_Score2', fontsize=16, weight='bold')
sns.despine()
plt.tight_layout()
plt.show()

fig.savefig('Graph/Scatter/ALK_{}_Endo_Correlation.pdf'.format(drug),dpi=300,bbox_inches='tight', transparent = True)


## Correlation plot ##


colordic = {'Resistance':"darkred",'Sensitive':'blue','Intermediate':'rosybrown'}


fig,axe=plt.subplots(figsize=(9,6))
for i in ['Resistance','Intermediate','Sensitive']:
    eachdf = df[df['Classification']==i]
    sns.scatterplot(data =eachdf, x='Resistance_Score1', y='Resistance_Score2',label=i, s=20, color= colordic[i])
#axe.spines['left'].set_linewidth(5)
#axe.spines['bottom'].set_linewidth(5)
axe.spines['left'].set_color('black')
axe.spines['bottom'].set_color('black')

plt.tick_params(labelsize = 20,size = 10,width=1)
plt.xticks([-5, 0, 5, 10])
plt.yticks([-5, 0, 5, 10])
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize = 20,markerscale=2.0,frameon=False)
sns.despine()
plt.tight_layout()
r,p = stats.pearsonr(df['Resistance_Score1'],df['Resistance_Score2'])
plt.title('{}_{}_{}'.format(r,p,df.shape[0]), size= 10) 

plt.show()
fig.savefig('Graph/Scatter/ALK_{}_classification.pdf'.format(drug),dpi=300,bbox_inches='tight', transparent = True)

### Drug comparison ###


### Comparison between drugs ##


alec = pd.read_csv('Result/ALK_A_Result1.csv', index_col=0)
lorla = pd.read_csv('Result/ALK_L_Result1.csv', index_col=0)
tpx = pd.read_csv('Result/ALK_T_Result1.csv', index_col=0)


dd = pd.concat([alec,lorla[['Resistance_Score','Classification']], 
                tpx[['Resistance_Score','Classification']]
                ],axis=1)


dd.columns = ['AA_Change', 'Num', 'PE4max-e_score', 'U6', 'sgRNA', 'Scaffold',
       'SynonyRTPBS', 'linker', 'TevopreQ1', 'polyT', 'SortingBC', 'Constant',
       'Buffer', 'RP', 'Exon', 'PBS', 'strand', 'SynonyCodon', 'seq', 'TTTT',
       'TTTTsg', 'template', 'IID', 'Group', 'oligo', 'oligo_len', 'peg',
       'PRIDICT', 'OldID', 'R1', 'IDseq', 'R2', 'REF', 'A', 'L', 'T', 'U',
       'Type', 'R1_RPM', 'R2_RPM', 'A_RPM', 'L_RPM', 'T_RPM', 'UN_RPM', 'LFC',
       'P', 'Odds', 'P_adjusted', 'num', 'syn_LFC', 'stand_LFC',
       'Resistance_Score1', 'Resistance_Score2', 'A_Score',
       'A_Class', 'L_Score', 'L_Class',
       'T_Score', 'T_Class']


dd = dd[['AA_Change','Type','A_Score', 'A_Class',
       'L_Score', 'L_Class',
       'T_Score', 'T_Class']]


alec_res = dd[dd['A_Class']=='Resistance']
lorla_res = dd[dd['L_Class']=='Resistance']
TPX_res = dd[dd['T_Class']=='Resistance']


al_res = alec_res[alec_res['L_Class']=='Resistance']
aoln_res = alec_res[alec_res['L_Class']!='Resistance']
loan_res = lorla_res[lorla_res['A_Class']!='Resistance']


# %%

fig,axe=plt.subplots(figsize=(6.5,6))


axe.spines['left'].set_color('black')
axe.spines['bottom'].set_color('black')

plt.tick_params(labelsize = 20,size = 5,width=1)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize = 20,markerscale=2.0,frameon=False)
sns.despine()
plt.tight_layout()

sns.scatterplot(data = dd, x='A_Score', y='L_Score',color='grey', s=20)
sns.scatterplot(data = al_res, x='A_Score', y='L_Score',color='darkred', label = "Both", s=20)
sns.scatterplot(data = aoln_res, x='A_Score', y='L_Score',color='darkblue', label = "Alec Only", s=20)
sns.scatterplot(data = loan_res, x='A_Score', y='L_Score',color='darkorange', label = "Lorla Only", s=20)

plt.axhline(y=loan_res["L_Score"].min(),linestyle='--',color='grey')
plt.axvline(x=aoln_res["A_Score"].min(),linestyle='--',color='grey')

plt.xticks([-5,0,5,10])
plt.yticks([-5,0,5,10])
r,p = stats.pearsonr(dd['A_Score'],dd['L_Score'])
plt.title('{}_{}_both{}_Alec{}_Lorla{}'.format(r,p,al_res.shape[0],aoln_res.shape[0],loan_res.shape[0]))

plt.show()

fig.savefig('Graph/Scatter/AL_comparison.pdf',dpi=300,bbox_inches='tight', transparent = True)


# %%

at_res = alec_res[alec_res['T_Class']=='Resistance']
aotn_res = alec_res[alec_res['T_Class']!='Resistance']
toan_res = TPX_res[TPX_res['A_Class']!='Resistance']



### A vsT ###

fig,axe=plt.subplots(figsize=(6.5,6))


axe.spines['left'].set_color('black')
axe.spines['bottom'].set_color('black')

plt.tick_params(labelsize = 20,size = 5,width=1)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize = 20,markerscale=2.0,frameon=False)
sns.despine()
plt.tight_layout()

sns.scatterplot(data = dd, x='A_Score', y='T_Score',color='grey', s=20)
sns.scatterplot(data = at_res, x='A_Score', y='T_Score',color='darkred', label = "Both", s=20)
sns.scatterplot(data = aotn_res, x='A_Score', y='T_Score',color='darkblue', label = "Alec Only", s=20)
sns.scatterplot(data = toan_res, x='A_Score', y='T_Score',color='darkorange', label = "TPX Only", s=20)

plt.axhline(y=toan_res["T_Score"].min(),linestyle='--',color='grey')
plt.axvline(x=aotn_res["A_Score"].min(),linestyle='--',color='grey')

plt.xticks([-5,0,5,10])
plt.yticks([-5,0,5,10])
plt.xlim([-8, 10.5])
r,p = stats.pearsonr(dd['A_Score'],dd['T_Score'])
plt.title('{}_{}_both{}_Alec{}_TPX{}'.format(r,p,at_res.shape[0],aotn_res.shape[0],toan_res.shape[0]))

plt.show()

fig.savefig('Graph/Scatter/AT_comparison.pdf',dpi=300,bbox_inches='tight', transparent = True)


### Drug map ###

################## Drugmap ###################

from scipy.stats import gaussian_kde
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


df = pd.read_csv("Result/ALK_T_Result1.csv")

for i in df.index:
    aa = df.loc[i,"AA_Change"]
    cds = aa[1:-1]
    df.loc[i,"CDS"] = cds
    
df['CDS'] = pd.to_numeric(df['CDS'], errors='coerce')

df['Resistance_Score'] = pd.to_numeric(df['Resistance_Score'], errors='coerce')

# Drop rows with invalid values
df = df.dropna(subset=[ 'CDS','Resistance_Score'])


# Sort data by Codon_Number to ensure correct plotting order
df = df.sort_values(by='CDS')



# Function to map classification to colors
def get_marker_color(classification):
    if classification == "Resistance":
        return "red"
    elif classification == "Intermediate":
        return "grey"
    elif classification == "Sensitive":
        return "blue"
    else:
        return "black"  # Default for unknown categories

# Plotting function with classification-based coloring
def plot_resistance_score(data, ax):
    # Get marker colors based on classification
    marker_colors = data['Classification'].apply(get_marker_color)

    # Create scatter plot
    scatter = ax.scatter(
        data['CDS'],
        data['Resistance_Score'],
        c=marker_colors,
        s=10,  # Marker size
        edgecolors="black",
        linewidths=0.4,
        alpha=0.9
    )

    # Add exon boundaries
    exon_boundaries = data.groupby('Exon')['CDS'].agg(['min', 'max']).reset_index()
    for _, exon in exon_boundaries.iterrows():
        ax.axvline(x=exon['max'], color='grey', linestyle='--', linewidth=1)

    # Add horizontal line at y=0
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_title('Resistance Score', fontsize=14)
    ax.set_ylabel('Resistance Score', fontsize=12)
    ax.set_xlabel('Codon Number', fontsize=12)

# Create a figure and subplot
fig, ax = plt.subplots(figsize=(18, 4))

# Plot resistance score
plot_resistance_score(df, ax)

plt.ylim([-5.5, 16])
plt.yticks([-5, 0, 5, 10, 15])
plt.show()




fig.savefig('Graph/Scatter/TPX_Drugmap.pdf',dpi=300,bbox_inches='tight', transparent = True)

### Sankey plot ###


        line=dict(color="black", width=1),
        label=all_nodes,
        color=[node_colors[node] for node in all_nodes]
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=counts['counts'].tolist(),
        color=link_colors  # Use the list of link colors
    )
)])

# Update layout settings (same as before)
fig.update_layout(
    title_text="CAIS to Classification Sankey Diagram",
    font_size=5,
    width=1000,
    height=800
)

fig.show()

fig.write_html("Graph/Sankey/TPX.html")



