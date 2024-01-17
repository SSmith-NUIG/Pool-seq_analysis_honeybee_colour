import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests

cols = ["CHR", "POS"]

df_1_1_2_2 = pd.read_csv(f"/home/stephen/Documents/colour_analysis/CMH/colour_1_1_2_2_3_3.cmh", sep="\t", header=None)
df_1_1_2_2.drop(columns=[2,3,4,5,6,7,8], inplace=True)
df_1_1_2_2.rename({0:"CHR", 1:"POS", 9:"pvalue"}, axis=1, inplace=True)
df_1_1_2_2['CHR_POS'] = df_1_1_2_2[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
p_adjusted = multipletests(df_1_1_2_2["pvalue"], method='bonferroni')
df_1_1_2_2["p_adjusted"] = p_adjusted[1]
df_1_1_2_2['minuslog10padjusted'] = -np.log10(df_1_1_2_2.p_adjusted)


def manhatten_plot(df):

    df.CHR = df.CHR.astype('category')
    df.CHR = df.CHR.cat.set_categories(["LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8",
                                                    "LG9", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15", "LG16"])

    df['ind'] = range(len(df))
    df_grouped = df.groupby('CHR')

    # Create plot for the whole genome
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []

    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10padjusted',color=colors[num % len(colors)], ax=ax, s=1)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, 30])

    # x axis label
    ax.set_xlabel('Chromosome')
    plt.title("CMH test for colour")
    # show the graph
    plt.show()

manhatten_plot(df_1_1_2_2)

df_1_2_2_3 = pd.read_csv(f"/home/stephen/Documents/colour_analysis/CMH/colour_1_2_2_3_3_1.cmh", sep="\t", header=None)
df_1_2_2_3.drop(columns=[2,3,4,5,6,7,8], inplace=True)
df_1_2_2_3.rename({0:"CHR", 1:"POS", 9:"pvalue"}, axis=1, inplace=True)
df_1_2_2_3['CHR_POS'] = df_1_2_2_3[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
p_adjusted = multipletests(df_1_2_2_3["pvalue"], method='bonferroni')
df_1_2_2_3["p_adjusted"] = p_adjusted[1]
df_1_2_2_3['minuslog10padjusted'] = -np.log10(df_1_2_2_3.p_adjusted)
manhatten_plot(df_1_2_2_3)

df_1_3_2_1 = pd.read_csv(f"/home/stephen/Documents/colour_analysis/CMH/colour_1_3_2_1_3_2.cmh", sep="\t", header=None)
df_1_3_2_1.drop(columns=[2,3,4,5,6,7,8], inplace=True)
df_1_3_2_1.rename({0:"CHR", 1:"POS", 9:"pvalue"}, axis=1, inplace=True)
df_1_3_2_1['CHR_POS'] = df_1_3_2_1[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
p_adjusted = multipletests(df_1_3_2_1["pvalue"], method='bonferroni')
df_1_3_2_1["p_adjusted"] = p_adjusted[1]
df_1_3_2_1['minuslog10padjusted'] = -np.log10(df_1_3_2_1.p_adjusted)

manhatten_plot(df_1_3_2_1)

df_merge_1 = pd.merge(df_1_1_2_2, df_1_2_2_3, how='inner', on="CHR_POS")
df_merge_1.drop(columns=["CHR_x", "POS_x", "CHR_y", "POS_y"], axis=1, inplace=True)
df_merge_1.rename({"pvalue_x": "pvalue_1", "pvalue_y":"pvalue_2",
                   "p_adjusted_x": "padjusted_1", "p_adjusted_y": "padjusted_2"}, axis=1, inplace=True)

df_merge_2 = pd.merge(df_merge_1, df_1_3_2_1, how="inner", on="CHR_POS")
df_merge_2.rename({"pvalue": "pvalue_3", "p_adjusted": "padjusted_3"}, axis=1, inplace=True)
df_merge_2.drop(columns=["ind_x", "ind_y", "ind"], inplace=True)
col_order = df_merge_2.columns.to_list()
col_order = [col_order[7]]+[col_order[8]]+[col_order[1]]+[col_order[0]]+[col_order[4]]+[col_order[9]]+[col_order[2]]+[col_order[5]]+[col_order[10]]
df_merge_2 = df_merge_2[col_order]


df_merge_2["pvalue_1"].fillna(1, inplace=True)
df_merge_2["pvalue_2"].fillna(1, inplace=True)
df_merge_2["pvalue_3"].fillna(1, inplace=True)
df_merge_2['average_pvalue'] = df_merge_2[['pvalue_1', 'pvalue_2', 'pvalue_3']].mean(axis=1)
df_merge_2['minuslog10pvalue'] = -np.log10(df_merge_2.average_pvalue)
df_merge_2["padjusted_1"].fillna(1, inplace=True)
df_merge_2["padjusted_2"].fillna(1, inplace=True)
df_merge_2["padjusted_3"].fillna(1, inplace=True)
df_merge_2["padjusted"] = df_merge_2[['padjusted_1', 'padjusted_2', 'padjusted_3']].mean(axis=1)
df_merge_2['minuslog10padjusted'] = -np.log10(df_merge_2.padjusted)
df_merge_2.to_csv("/home/stephen/df_merge2.csv")


# bonferoni padjusted
final_padjusted = df_merge_2[["CHR", "POS", "CHR_POS", "padjusted"]]
final_padjusted['minuslog10padjusted'] = -np.log10(final_padjusted.padjusted)
final_padjusted.CHR = final_padjusted.CHR.astype('category')
final_padjusted.CHR = final_padjusted.CHR.cat.set_categories(["LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8",
                                                "LG9", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15", "LG16"])

final_padjusted['ind'] = range(len(final_padjusted))
final_padjusted_grouped = final_padjusted.groupby('CHR')
final_padjusted.to_csv("/home/stephen/padjustedCMH.csv", index=False)
final_padjusted = pd.read_csv("/home/stephen/padjustedCMH.csv")

def abline(slope, intercept):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')
manhatten_plot(final_padjusted)
abline(0,7)


## WHOLE GENOME STUDY 
all_LG_significant = final_padjusted[final_padjusted["minuslog10padjusted"] > 7]

amGTF = pd.read_csv(f"/home/stephen/Downloads/GO_files/for_cmh.gtf", sep="\t", header=None)
# read in the NCBI honeybee genes file
gene_names = pd.read_csv("/home/stephen/Downloads/honey_bee_genes_NCBI_2022.tsv", sep="\t")

# more accurate filtering step
# get the feature where the genomic location of significance lies
# get GTF info from locations of significant CMH results
all_dfs_list = []
for chr in range(1,17):
    LG = "LG"+str(chr)
    print(LG)
    gtf_loop = amGTF[amGTF[0] == LG]
    #print(gtf_loop.columns)
    gene_info_loop = gtf_loop[8].str.split(';', expand=True)
    gene_info_loop["Start"] = gtf_loop[3]
    gene_info_loop["End"] = gtf_loop[4]
    gene_info_loop[["gene", "Name"]] = gene_info_loop[0].str.split(' ', expand=True)
    gene_info_loop["Name"] = gene_info_loop["Name"].str.replace('"', '')
    gene_names_loop = gene_names[gene_names["Chromosomes"]==LG]
    gene_info_merged_loop = pd.merge(gene_info_loop, gene_names_loop, how="left", left_on="Name", right_on='Symbol')

    significant_loop = all_LG_significant[all_LG_significant["CHR"]==LG]

    df_list = []
    for v in significant_loop['POS']:
        index_of_value = significant_loop.index[significant_loop['POS'] == v].tolist()
        df3 = gene_info_merged_loop[(gene_info_merged_loop['Start'] < v) &
                                    (gene_info_merged_loop['End'] > v)]
        df3["pvalue"] = significant_loop["minuslog10pvalue"][index_of_value[0]]
        df3["POS"] = significant_loop["POS"][index_of_value[0]]

        df_list.append(df3)
    if len(df_list) > 1:
        combined_LG = pd.concat(df_list, axis=0)
        combined_LG = combined_LG.sort_values(['Name']).reset_index(drop=True)
        combined_LG.drop_duplicates(inplace=True)
        combined_LG.to_csv(f"/home/stephen/Documents/colour_analysis/CMH/{LG}_region_gene_info.csv_new")
        all_dfs_list.append(combined_LG)

all_LGS_gene_info = pd.concat(all_dfs_list, axis=0)
all_LGS_gene_info.reset_index(inplace=True)
all_LGS_gene_info.to_csv("/home/stephen/Documents/colour_analysis/CMH/all_LG_gene_info_new.csv",index=False)

gowinda_significant_all_genome = all_LGS_gene_info[["Chromosomes", "POS"]]
gowinda_significant_all_genome.drop_duplicates(inplace=True)
gowinda_significant_all_genome.to_csv("/home/stephen/gowinda/gowinda_significant_all_genome_fixed_new.csv",sep="\t", index=False, header=None)
