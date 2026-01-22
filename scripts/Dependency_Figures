
drug = "UN"

lst = []
for exon in ['E22','E23']:
        
    df = pd.read_csv('Endo_Count/H2228/H2228_{}_R1.csv'.format(exon))
    df.reset_index(drop=True, inplace=True)
    df = rpm(df)

    df['LFC'] = np.log2((df['{}_RPM'.format(drug)]+1)/(df['R1_RPM']+1))
    #df['LFC'] = np.log2((df['UN_RPM']+1)/(df['{}_RPM'.format(rep)]+1))

    df= df[df["R1_RPM"]>0.5]
    lst.append(df)


e2 = pd.concat(lst)


lst = []

for exon in ['E22','E23']:
        
    df = pd.read_csv('Endo_Count/H2228/H2228_{}_R2.csv'.format(exon))
    df.reset_index(drop=True, inplace=True)
    df = rpm(df)

    df['LFC'] = np.log2((df['{}_RPM'.format(drug)]+1)/(df['R1_RPM']+1))
    #df['LFC'] = np.log2((df['UN_RPM']+1)/(df['{}_RPM'.format(rep)]+1))

    df= df[df["R1_RPM"]>0.5]
    lst.append(df)


e3 = pd.concat(lst)

  e2.to_csv("Result/Raw_2228/U_R1.csv")
e3.to_csv("Result/Raw_2228/U_R2.csv")


# %%

df = pd.read_csv("Result/Raw_2228/U_R1.csv", index_col=0)

def calculate_weighted_score(group, drug_column):
    # Calculate total ARV (allele frequency)
    total_drug = group[drug_column].sum()
    # Calculate weighted LFC (Resistance Score)
    weighted_lfc = (group[drug_column] * group['stand_LFC']).sum() / total_drug
    # Assign the score to a new column for all rows in the group
    group['Resistance_Score'] = weighted_lfc
    #group['Fitness_Score'] = weighted_lfc
    return group


# Apply the calculation across groups of IDs
df = df.groupby('ID',as_index=False).apply(calculate_weighted_score, "R1")
df = df.reset_index(drop=True) 


df.to_csv("Result/ALK_U_R1.csv", index=False)




  
dd = pd.read_csv("Result/Processed_2228/U_Result1.csv")
## AR Dependency ROC and PR ##

e2 = dd.copy()

dic  = {'Syn':1,'Nonsense':0}
e2 = e2[e2['Type'].isin(['Nonsense','Syn'])]  
e2['label'] = [dic[i] for i in e2['Type']]

label =e2['label']
value = e2['LFC']

ef, et, thresholds = roc_curve(label, value)
e_auc = metrics.roc_auc_score(label, value)

# Youden's Index calculation
youden_index = et - ef
optimal_idx = np.argmax(youden_index)
optimal_threshold = thresholds[optimal_idx]
optimal_sensitivity = et[optimal_idx]
optimal_specificity = 1 - ef[optimal_idx]

# Plot ROC curve with optimal threshold
fig, axe = plt.subplots(figsize=(5, 5))
plt.plot(ef, et, color = '#2ECC71', label=f'ROC Curve\nAUC: {e_auc:.2f}', linewidth=3)

plt.plot([0, 1], [0, 1], '--', color='gray', lw=2)

plt.tick_params(labelsize = 20,size = 10,width=1)
#plt.tick_params(axis='both', which='both', direction='out', length=6, width=2, colors='black')

plt.xlabel('1 - Specificity (FPR)', fontsize=15)
plt.ylabel('Sensitivity (TPR)', fontsize=15)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize = 20,markerscale=2.0,frameon=False)
plt.grid(False)
## Spine
#axe.spines['left'].set_linewidth(5)
#axe.spines['bottom'].set_linewidth(5)
axe.spines['left'].set_color('black')
axe.spines['bottom'].set_color('black')
sns.despine()


label=f'Optimal Threshold: {optimal_threshold:.2f}\nSensitivity: {optimal_sensitivity:.2f}\nSpecificity: {optimal_specificity:.2f}'


dd = ''
for i in list(set(e2['Type'])):
    eachdf = e2[e2['Type']==i]
    dd += i
    dd += ':'
    dd += str(eachdf.shape[0])
    dd += '_'   

dd += label

plt.title(dd)

plt.show()


fig.savefig('Graph/ROC/H2228_Dependency_ROC.pdf',dpi=300,bbox_inches='tight', transparent = True)



  
### 비교 ROC ###

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn import metrics
import seaborn as sns

fig, ax = plt.subplots(figsize=(5, 5))

dic = {'Syn':1, 'Nonsense':0}

############################
️ Dataset 1 
############################
dd1 = pd.read_csv("Result/Processed_2228/U_Result1.csv")
e1 = dd1[dd1['Type'].isin(['Nonsense', 'Syn'])].copy()
e1['label'] = e1['Type'].map(dic)

fpr1, tpr1, _ = roc_curve(e1['label'], e1['LFC'])
auc1 = metrics.roc_auc_score(e1['label'], e1['LFC'])

ax.plot(
    fpr1, tpr1,
    color='#2ECC71',
    lw=3,
    label=f'H2228 (AUC = {auc1:.2f})',
    zorder=10
)

############################
 Dataset 2 
############################
dd2 = pd.read_csv("Result/ALK_U_Result1.csv")
e2 = dd2[dd2['Type'].isin(['Nonsense', 'Syn'])].copy()
e2['label'] = e2['Type'].map(dic)

fpr2, tpr2, _ = roc_curve(e2['label'], e2['Fitness_Score'])
auc2 = metrics.roc_auc_score(e2['label'], e2['Fitness_Score'])

ax.plot(
    fpr2, tpr2,
    color='#123753',          
    lw=3,
   # linestyle='--',       
    label=f'H3122 (AUC = {auc2:.2f})'
)


  
