

def lws(df):
    # CDS부분에 num column
    
    for i in df.index:
        AA = df.loc[i,"AA_Change"]
        if df.loc[i,"Type"] == "Splice_S5":
            pass
        elif df.loc[i,"Type"] == "Splice_S3":
            pass
        elif df.loc[i,"Type"] =="Deletion":
            pass
        elif df.loc[i,"Type"] =="Duplication":
            pass
        else:
            df.loc[i,"num"] = int(AA[1:-1])


    maxx = df['num'].max()
    minn = df['num'].min()

    for i in df.index:
        AA = df.loc[i,"AA_Change"]
        if df.loc[i,"Type"] == "Splice_S5":
            df.loc[i,"num"] = minn
        elif df.loc[i,"Type"] == "Splice_S3":
            df.loc[i,"num"] = maxx
        elif df.loc[i,"Type"] == "Deletion":
            de = df[df["Num"]==df.loc[i,"Num"]]
            de.reset_index(drop=True, inplace=True)
            df.loc[i,"num"] = de.loc[0,"num"]
        elif df.loc[i,"Type"] == "Duplication":
            de = df[df["Num"]==df.loc[i,"Num"]]
            de.reset_index(drop=True, inplace=True)
            df.loc[i,"num"] = de.loc[0,"num"]

    df['num'] = df['num'].astype(int)
    df['LFC'] = df['LFC'].astype(float)

    
    # LOWESS regression 실질적인부분
    
    dd = df[df["Type"] =="Syn"]
    # Data generation
    x = np.array(dd["num"])
    y = np.array(dd["LFC"])
    # Model fitting
    lowess_model = lowess.Lowess()
    lowess_model.fit(x, y)   # fitting
    
    # Model prediction
    x_pred = np.linspace(int(maxx), int(minn), int(maxx-minn)+1)  #처음,끝, 안에 표시할 갯수
    y_pred = lowess_model.predict(x_pred)

    x_lst = list(x_pred)
    y_lst = list(y_pred)
    ddic = {}
    for i in range(len(x_lst)):
        ddic[x_lst[i]] = y_lst[i]    
    
    
    # syn_LFC column에 각 number position 별 lowess regression된 syn_LFC값을 넣어줌
    for i in df.index:
        num = df.loc[i,"num"]
        df.loc[i,"syn_LFC"] = ddic[num]
    
    
    # normalization 단계에서 각 (LFC value - syn_LFC)/ df_sd 를 해줬음. 
    
    df_sd = np.std(df[df["Type"] == "Syn"]["LFC"])  #얘네는 regression하지않은 syn 애들의 실제 sd값이므로
    
    df["stand_LFC"] = df["LFC"].sub(df["syn_LFC"]).mul(1/df_sd)

    return df


def d10_fisher(df, rep):
    df = df.copy()
    wt_r1 = df.loc[len(df)-1, rep] + 1
    wt_ref = df.loc[len(df)-1, 'REF'] + 1
    
    p_values = []  # To store raw p-values for adjustment
    
    for idx in df.index:
        r1 = df.loc[idx, rep] + 1
        ref = df.loc[idx, 'REF'] + 1
        
        o, p = stats.fisher_exact([[r1, wt_r1], [ref, wt_ref]])
        df.loc[idx, 'P'] = p
        df.loc[idx, 'Odds'] = o
        p_values.append(p)
    
    # Perform FDR adjustment
    _, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    # Add adjusted p-values to the dataframe
    df['P_adjusted'] = p_adjusted
    
    return df

def classification(df):   #old version
    df = df.copy()
    syn = df[df['Type']=='Syn']
    co1 = np.percentile(syn['Resistance_Score1'],99.7)  #two-tailed 3sigma 
    co2 = np.percentile(syn['Resistance_Score2'],99.7)

    co3 = np.percentile(syn['Resistance_Score1'],95)  #95 
    co4 = np.percentile(syn['Resistance_Score2'],95)  #two-tailed 2sigma


    hit = df[df['Resistance_Score1']>co1]
    hit = hit[hit['Resistance_Score2']>co2]
    hit['Classification']  = 'Resistance'

    sens = df[df['Resistance_Score1']<=co3]
    sens = sens[sens['Resistance_Score2']<=co4]
    sens['Classification'] = 'Sensitive'

    merged = pd.concat([hit,sens])
    inter = np.setdiff1d(df.index,merged.index)
    inter = df.loc[inter,:]
    inter['Classification'] = 'Intermediate'
    final = pd.concat([merged,inter])
    return final


def rpm(df):
    for i in df.index:
        AA = df.loc[i,"AA_Change"]
        if AA[0] == AA[-1]:
            df.loc[i,"Type"] = "Syn"
        elif AA[-1] =="*":
            df.loc[i,"Type"] = "Nonsense"
        elif "S5_" in AA:
            df.loc[i,"Type"] ="Splice_S5"
        elif "S3_" in AA:
            df.loc[i,"Type"] ="Splice_S3"
        elif 'del' in AA:
            df.loc[i,"Type"] ="Deletion"
        elif 'dup' in AA:
            df.loc[i,"Type"] ="Duplication"
        else:
            df.loc[i,"Type"] ="Missense"

    df['R1_RPM'] = df['R1']/df['R1'].sum()*1e6+1
    df['R2_RPM'] = df['R2']/df['R2'].sum()*1e6+1
    df['A_RPM'] = df['A']/df['A'].sum()*1e6+1
    df['L_RPM'] = df['L']/df['L'].sum()*1e6+1
    df['T_RPM'] = df['T']/df['T'].sum()*1e6+1
    df['UN_RPM'] =  df['U']/df['U'].sum()*1e6+1
    #df = df[df['R1']!=0]

    return df


drug = "T"
rep = 'R1'

lst = []

for exon in ['E20','E21','E22','E23','E24','E25','E26','E27','E28']:
    df = pd.read_csv('Endo_Count/ALK_{}_{}_D21.csv'.format(exon,rep))
    df.reset_index(drop=True, inplace=True)
    df = rpm(df)

    #df['LFC'] = np.log2((df['{}_RPM'.format(drug)]+1)/(df['UN_RPM']+1))
    df['LFC'] = np.log2((df['UN_RPM']+1)/(df['{}_RPM'.format(rep)]+1))
    df = d10_fisher(df,rep)
    df= df[df["R1_RPM"]>0.5]
    df = df[df['P_adjusted']<0.05]
    df = df[df['Odds']>3]
    df = lws(df)
    lst.append(df)


e2 = pd.concat(lst)

rep = 'R2'
lst = []


for exon in ['E20','E21','E22','E23','E24','E25','E26','E27','E28']:
    df = pd.read_csv('Endo_Count/ALK_{}_{}_D21.csv'.format(exon,rep))
    df.reset_index(drop=True, inplace=True)
    df = rpm(df)

    #df['LFC'] = np.log2((df['{}_RPM'.format(drug)]+1)/(df['UN_RPM']+1))
    df['LFC'] = np.log2((df['UN_RPM']+1)/(df['{}_RPM'.format(rep)]+1))
    df = d10_fisher(df,rep)
    df= df[df["R1_RPM"]>0.5]
    df = df[df['P_adjusted']<0.05]
    df = df[df['Odds']>3]
    df = lws(df)
    lst.append(df)


e3 = pd.concat(lst)



df = pd.read_csv("Result/Raw/A_R1.csv", index_col=0)


## exmaple 

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


df.to_csv("Result/ALK_A_R1.csv", index=False)


### Example of calculaion ##

df1 = pd.read_csv("Result/Processed/ALK_A_R1.csv")
df2 = pd.read_csv("Result/Processed/ALK_A_R2.csv")

df1 = df1.drop_duplicates("ID")
df2 = df2.drop_duplicates("ID")


df1 = df1.set_index("ID")

df2 = df2.set_index("ID")

dd = pd.concat([df1,df2["Resistance_Score"]], axis= 1, join = 'inner')

dd.columns = ['AA_Change', 'Num', 'PE4max-e_score', 'U6', 'sgRNA', 'Scaffold',
       'SynonyRTPBS', 'linker', 'TevopreQ1', 'polyT', 'SortingBC', 'Constant',
       'Buffer', 'RP', 'Exon', 'PBS', 'strand', 'SynonyCodon', 'seq', 'TTTT',
       'TTTTsg', 'template', 'IID', 'Group', 'oligo', 'oligo_len', 'peg',
       'PRIDICT', 'OldID', 'R1', 'IDseq', 'R2', 'REF', 'A', 'L', 'T', 'U',
       'Type', 'R1_RPM', 'R2_RPM', 'A_RPM', 'L_RPM', 'T_RPM', 'UN_RPM', 'LFC',
       'P', 'Odds', 'P_adjusted', 'num', 'syn_LFC', 'stand_LFC',
       'Resistance_Score1','Resistance_Score2']

dd["Resistance_Score"] = (dd["Resistance_Score1"] + dd["Resistance_Score2"])/2

## Classification ##
dd = classification(dd)



dd.to_csv('Result/ALK_A_Result1.csv')






