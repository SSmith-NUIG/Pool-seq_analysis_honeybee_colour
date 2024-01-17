import pandas as pd

# Apis mellifera Gene Ontology (GO) files downloaded from here:
# https://hymenoptera.elsiklab.missouri.edu/hgd-go-annotation
# gene name column and GO ID column were selected and pasted into a new csv file named gowinda_test.csv
# rows which had their gene name attributed to uniprot were discarded as these were found to be incorrect
# a lot of the time
# This could be done through python by selecting which columns you want to read in with read_csv and then filtering
# for lines which dont contain uniprot in the first column
# This file was then joined with a file downloaded from quickGO
# https://www.ebi.ac.uk/QuickGO/annotations?goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:0009653&taxonId=7460&taxonUsage=descendants
# filtering by taxon ID for apis mellifera 7460 and filtering out the GO:ID and gene name as before

gowinda = pd.read_csv("/home/stephen/Downloads/GO_files/HGD_input_7.csv", sep="\t", header=None)

gowinda.drop_duplicates(inplace=True, keep="first")

gowinda_extra = pd.read_csv("/home/stephen/Downloads/GO_files/gowinda_extra2.csv", sep="\t", header=None)
gowinda_extra.drop_duplicates(inplace=True, keep="first")
gowindas = [gowinda_extra, gowinda]
gowinda = pd.concat(gowindas)
gowinda
# group genes which are associated to the same GO_IDs together into a new column
# and separate them with a space
names = gowinda.groupby([0])[1].apply(' '.join).reset_index()
# Join it to the original dataset
gowinda2 = gowinda.merge(names, on=0)
# Drop the '1_x' column then drop duplicates.
gowinda2 = gowinda2.drop(columns=["1_x"]).drop_duplicates()
gowinda2.rename(columns={0: "GO_ID"}, inplace=True)
# make sure the entire dataframe is in the string format and there are no trailing spaces
# this is required for merging later on
gowinda2 = gowinda2.astype(str)
gowinda2.GO_ID = gowinda2.GO_ID.str.strip()

# Gene ontology annotation file was downloaded from here using wget:
# wget http://purl.obolibrary.org/obo/go/go-basic.obo
# This file was then filtered to only contain the lines which had the GO:IDs and the GO descriptions
# This was done on the command line with the following command:
# grep -w 'id:\|name:' go-basic.obo > filtered_go.obo
go_description = pd.read_csv("/home/stephen/Downloads/GO_files/filtered_go.obo", header=None, sep="\t")

# Above file is in the format:
# id: GO_ID
# name: GO_DESCRIPTION
# id: GO_ID
# name: GO_DESCRIPTION
# etc. etc...
# create new columns GO:ID and DESCRIPTION from every other row respectively
go_description = pd.DataFrame({'GO_ID': go_description[0].iloc[::2].values, 'DESCRIPTION': go_description[0].iloc[1::2].values})

# clean up the cell values in each column by removing the unneeded strings "id:" and "name:"
go_description['GO_ID'] = go_description['GO_ID'].str.replace('id:', '')
go_description['DESCRIPTION'] = go_description['DESCRIPTION'].str.replace('name:', '')

# ensure dataframe is in string format and has no trailing spaces
go_description = go_description.astype(str)
go_description.GO_ID = go_description.GO_ID.str.strip()

# Left join on the gowinda2 dataframe which contains the apis mellifera GO terms and genes
# with the go_description dataframe which contains the GO:IDs and Descriptions for every GO term
go_id_description_genes = pd.merge(gowinda2, go_description, how="left", on="GO_ID")

# reorder the columns in the dataframe to be compatible with gowinda
# The format required is: GO_ID GO_DESCRIPTION GENES_INVOLVED
cols = go_id_description_genes.columns
cols = [cols[0]] + [cols[2]] + [cols[1]]
go_id_description_genes = go_id_description_genes[cols]

# write final dataframe to file with no header and set delimiter as tab
go_id_description_genes.to_csv("/home/stephen/Downloads/GO_files/gowinda_go_input_HGD6.csv", sep="\t", header=None, index=False)
