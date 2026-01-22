import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, multivariate_logrank_test



sv = pd.read_csv("MSK/SampleBC/EML4_ALK_ALL.csv")
df = pd.read_csv("MSK/msk_chord_2024/data_clinical_patient.txt", sep = '\t')
df2 = pd.read_csv("MSK/msk_chord_2024/data_clinical_sample.txt", sep = '\t')

lst = list(sv["Sample_Id"])

dd = df2[df2["SAMPLE_ID"].isin(lst)]

# %%

mut = pd.read_csv("MSK/msk_chord_2024/ALK_mutation_all.csv")

dic = {}
for i in mut.index:
    ID = mut.loc[i,"Tumor_Sample_Barcode"]
    if ID not in dic:
        dic[ID] = mut.loc[i,"HGVSp_Short"]
    
    else:
        dic[ID] = dic[ID] + ',' + mut.loc[i,"HGVSp_Short"]
        

# %%

for i in dd.index:
    ID = dd.loc[i,"SAMPLE_ID"]
    if ID in dic:
        dd.loc[i,"Mutation"] = dic[ID]
    else:
        dd.loc[i,"Mutation"] = "Intact"

# %%


treatment_df = pd.read_csv("MSK/2024,Nature/msk_chord_2024/data_timeline_treatment.txt", sep = '\t')
sequencing_df = pd.read_csv("MSK/2024,Nature/msk_chord_2024/data_timeline_specimen.txt", sep = '\t')


lst = list(dd["PATIENT_ID"])

treatment_df = treatment_df[treatment_df["PATIENT_ID"].isin(lst)]
sequencing_df = sequencing_df[sequencing_df["PATIENT_ID"].isin(lst)]

## MSK 2024, All 0 timing



# %%

## Alectinib focus ##
survival_df = dd.copy()

treatment_df=treatment_df[treatment_df["PATIENT_ID"].isin(list(survival_df["PATIENT_ID"]))]
treatment_df=treatment_df[treatment_df["AGENT"] =="ALECTINIB"]

dic = {}
for i in treatment_df.index:
    ID = treatment_df.loc[i,"PATIENT_ID"]
    dic[ID]=int(treatment_df.loc[i,"START_DATE"])

for i in survival_df.index:
    ID = survival_df.loc[i,"PATIENT_ID"]
    if ID not in dic:
        survival_df.loc[i,"TIMING"] = "Not_Used"
    else:
            
        num = dic[ID]
        if num>-10: 
            survival_df.loc[i,"TIMING"] = "PRE"
        elif num <-10 :
            survival_df.loc[i,"TIMING"] = "POST"

# %%

## Lorlatinib focus ##
survival_df = dd.copy()

treatment_df=treatment_df[treatment_df["PATIENT_ID"].isin(list(survival_df["PATIENT_ID"]))]
treatment_df=treatment_df[treatment_df["AGENT"] =="LORLATINIB"]

dic = {}
for i in treatment_df.index:
    ID = treatment_df.loc[i,"PATIENT_ID"]
    dic[ID]=int(treatment_df.loc[i,"START_DATE"])

for i in survival_df.index:
    ID = survival_df.loc[i,"PATIENT_ID"]
    if ID not in dic:
        survival_df.loc[i,"TIMING"] = "Not_Used"
    else:
            
        num = dic[ID]
        if num>=-10: 
            survival_df.loc[i,"TIMING"] = "PRE"
        elif num <-10 :
            survival_df.loc[i,"TIMING"] = "POST"




#df = pd.read_csv("MSK/EML4ALK_Alectinib.csv")



# Extract relevant columns and drop missing values
df_survival = df[['OS_MONTHS', 'Status', 'A_Classification']].dropna()
df_survival.columns = ["OS_Months", "Status","Classification"]

# Ensure 'Status' column is coded as 1 = Event (death), 0 = Censored
df_survival['Status'] = df_survival['Status'].astype(int)

# Define classification groups and custom colors
groups = ['Resistant',  'Intact', 'Sensitive']
colors = {'Resistant': '#8B0000', 'Intermediate': '#BC8F8F', 'Sensitive': '#0000FF',"Intact":"lightgrey"}

# Set up plot aesthetics
sns.set_style("white")
fig, axe = plt.subplots(figsize=(8, 6))

# Kaplan-Meier fitting and plotting for each group
kmf_dict = {}  # Store KMF objects for risk table
axe = plt.gca()

for group in groups:
    mask = df_survival['Classification'] == group
    kmf = KaplanMeierFitter()
    kmf.fit(df_survival.loc[mask, 'OS_Months'], df_survival.loc[mask, 'Status'], label=group)
    kmf.plot_survival_function(ax=axe, color=colors[group], lw=2, ci_show=False)  # Remove confidence interval
    kmf_dict[group] = kmf  # Store for risk table



# ---- Step 3: Pairwise Post-Hoc Log-Rank Tests ---- #

# Resistant vs Sensitive
mask_resistant = df_survival['Classification'] == 'Resistant'
mask_intact = df_survival['Classification'] == 'Intact'


results_resistant_vs_intact = logrank_test(
    df_survival.loc[mask_resistant, 'OS_Months'], df_survival.loc[mask_intact, 'OS_Months'],
    df_survival.loc[mask_resistant, 'Status'], df_survival.loc[mask_intact, 'Status']
)
print(f"Posthoc Log-rank test p-value (Resistant vs intact): {results_resistant_vs_intact.p_value:.30f}")





# Customize the plot
plt.xlabel("Overall Survival (Months)")
plt.ylabel("Survival Probability")
# plt.title(f"Posthoc Log-rank test p-value (Resistant vs Sensitive): {results_resistant_vs_sensitive.p_value:.4f}"
#           f"Posthoc Log-rank test p-value (Intermediate vs Sensitive): {results_intermediate_vs_sensitive.p_value:.4f}")

plt.legend(title="Classification", bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize = 20,markerscale=2.0,frameon=False)

plt.grid(False)  # Remove grid
sns.despine()  # Remove axis spines

# Set X-axis ticks every 20 months
max_time = df_survival['OS_Months'].max()
plt.xticks(range(0, int(max_time) + 20, 20))
plt.ylim(0, 1.05)
plt.yticks([0,0.25,0.5,0.75,1.0])
# Add number at risk table (corrected)
add_at_risk_counts(*kmf_dict.values(), ax=axe)

# Show the plot
plt.show()

fig.savefig('Graph/KM/EML4ALK_Alectinib_2.pdf', format='pdf', transparent=True)



