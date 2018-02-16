import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import collections
import itertools

path = "C:/Users/"

dataset = pd.read_csv(path + "mydata.csv", sep = ",")
labels = pd.read_csv(path + "labels2.csv", sep = ",")

# change the name of the first column to genes
dataset = dataset.rename(columns = {'Unnamed: 0':'genes'})

# row names index is equal to numbers. we need to set the data from the 1st column as rownames
# Set the 1st column as rownames
dataset.index = dataset['genes']


# delete first column
del dataset['genes']



print(dataset)
print(labels)


           ##### find gene expression for every cell

#transpose the dataset (rownames = cell names)

dataset1 = dataset.transpose()
x = dataset1.columns.values

counts = dataset1.loc[:,x].values        #counts give us nd.arrays with the number of transcript for every cell

#####Remove the comments to make a csv file with gene expression for every cell

# with open(path + 'cells.csv','w', newline = '') as f1:
#     writer = csv.writer(f1)
#
#     j=0                                 #counter that gives row index
#     for count in counts:
#         cells = []
#         i = 0
#         for number in count:
#             if number > 0:`
#                 i = i +1
#         #print(j, i)
#         j = j+1
#         cells.append(i)
#         writer.writerow(cells)

gene_exp = pd.read_csv(path + 'cells.csv', header = None)
gene_exp = gene_exp.rename(columns = {0: 'total_number_of_genes'})
print(gene_exp)

#pd.DataFrame of gene expression for every cell along with labels (cell name, cell type)
gene_exp1 = pd.concat([gene_exp,labels[['cell_id']],labels[['cell_type']]],axis = 1)    #add cell name and cell type to the pd.Dataframe
print(gene_exp1)

###### after that I want to plot these values in a histogram to get an overall idea

plt.hist(gene_exp1['total_number_of_genes'], facecolor = 'green', bins = 100)
plt.xlabel('Number of genes')
plt.ylabel('Frequency')
plt.title('Expressed genes across 466 cells')
plt.show()                                         ###How many genes do the cells usually express


#Find average gene expression for each cell type
#Subset the data into cell types (9)

oligo = gene_exp1.query('cell_type == "oligodendrocytes"')
neurons = gene_exp1.query('cell_type == "neurons"')
OPC = gene_exp1.query('cell_type == "OPC"')
microglia = gene_exp1.query('cell_type == "microglia"')
hybrid = gene_exp1.query('cell_type == "hybrid"')
astrocytes = gene_exp1.query('cell_type == "astrocytes"')
endothelial = gene_exp1.query('cell_type == "endothelial"')
fetalq = gene_exp1.query('cell_type == "fetal_quiescent"')
fetalr = gene_exp1.query('cell_type == "fetal_replicating"')



#### Calculate gene expression % (total number of genes = 21627)

k = 21627
def freq (a,b):
    return (a/b)*100 #where a = average number of genes for each type and b = k (total number of genes)

opcs = freq(OPC['total_number_of_genes'].mean(),k)
oligos = freq(oligo['total_number_of_genes'].mean(),k)
neur = freq(neurons['total_number_of_genes'].mean(),k)
glia = freq(microglia['total_number_of_genes'].mean(),k)
astro = freq(astrocytes['total_number_of_genes'].mean(),k)
endo = freq(endothelial['total_number_of_genes'].mean(),k)
fetq = freq(fetalq['total_number_of_genes'].mean(),k)
fetr = freq(fetalr['total_number_of_genes'].mean(),k)
hyb = freq(hybrid['total_number_of_genes'].mean(),k)

#create a list with the cell type, average number of genes and the % of the expressed genes
list = [['OPC', OPC['total_number_of_genes'].mean(), opcs], ['oligodendrocytes',oligo['total_number_of_genes'].mean(),oligos], ['neurons',neurons['total_number_of_genes'].mean(),neur],['microglia',microglia['total_number_of_genes'].mean(),glia],['astrocytes',astrocytes['total_number_of_genes'].mean(),astro], ['endothelial',endothelial['total_number_of_genes'].mean(),endo],['fetal quiescent',fetalq['total_number_of_genes'].mean(),fetq],['fetal_replicating',fetalr['total_number_of_genes'].mean(),fetr],['hybrid', hybrid['total_number_of_genes'].mean(),hyb]]
print(list)

#Make a pandas dataframe from the list and add column names
col = ['Cell_type','Average number of genes','% of expressed genes']
new_list = pd.DataFrame.from_records(list, columns = col)
print(new_list)




#### It gives the number of cells in which every gene is expressed

# This applies for the non transposed dataset
y = dataset.columns.values
genes = dataset.loc[:,y].values

# with open(path + 'genes.csv','w', newline = '') as f2:
#     writer = csv.writer(f2)
#     j=0                                 #counter that gives the cell number
# #print(genes)
#     for gene in genes:
#
#         i = 0
#         for number in gene:
#             if number > 0:
#                 i = i +1
#         #print(j, i)
#         j = j+1
#         alist.append(i)
#         writer.writerow(alist)


cell_gen = pd.read_csv(path + 'genes.csv', header = None)

cell_gen = cell_gen.rename(columns = {0: 'Number of cells'})

x2 = pd.DataFrame(data = x, columns = ['Genes'])         #where X all the gene names

# pd.Dataframe which includes the genes and in how many cells are expressed
cell_gen1 = pd.concat([x2,cell_gen],axis = 1)
print(cell_gen1)

###### Plot my results in a histogram

plt.hist(cell_gen1['Number of cells'], facecolor = 'green', bins = 100)
plt.xlabel('Number of cells')
plt.ylabel('Number of genes')
plt.title('Expressed genes across 466 cells')
plt.show()



######### Find max expression (in terms of highest number of transcripts) for every cell and identify the gene with the highest transcript number

#### use the counts from datataset1

listA = []                       #two empty lists to extract the results
listB = []
for count in counts:                #for every row, i.e. for every cell find the maximum number of transcripts
    position = 0
    max_exp = 0
    i = 0
    for number in count:

        if number > max_exp:
            max_exp = number
            position = i
        i = i + 1
    print(position, max_exp)
    listA.append(position)
    listB.append(max_exp)




### since we got the gene with the most counts for every cell and its position, let's find which one is it
### in the cell_gen1 pd.DataFrame, we have an index next to every gene. We can use this information to identify them


gene_names = []
for position in listA:                    #by using the position of listA, I try to find the index of cell_gen1 which will give me the gene name
    for value in cell_gen1.index.values:
        if position == value:
            gene_names.append(cell_gen1.iloc[position,0])         #gives the gene name



maximum = pd.DataFrame({'Gene':gene_names,'Maximum number of transcripts':listB})
####A pandas dataframe contaning the genes with the maximum number of transcripts for each cell
maximum2 = pd.concat([maximum, labels['cell_id'], labels['cell_type']],axis = 1)
print(maximum2)

##### For every cell print a histogram demonstrating the number of counts

# for count in counts:                                 ### from what I can see, for every cell most of the genes has low numbers of counts (close to 1)
#     plt.hist(count, bins = 100, facecolor = 'green')
#     plt.xlabel('Number of counts')
#     plt.show()


##Find in which cell do every gene has the highest expression (use the genes from dataset)


listC = []                       #two empty lists to extract the results
listD = []
for gene in genes:
    position = 0
    max_exp = 0
    i = 0
    for number in gene:

        if number > max_exp:
            max_exp = number
            position = i
        i = i + 1
    print(position, max_exp)
    listC.append(position)
    listD.append(max_exp)

#find the cells in which every gene shows the highest expression levels (number of counts)

cell_names = []
for position in listC:
    for value in gene_exp1.index.values:
        if position == value:
            cell_names.append(gene_exp1.iloc[position,1])



maximumcell = pd.DataFrame({'Cell':cell_names,'Maximum number of transcripts':listD})

#pandas dataframe, which contains the cells in which every gene shows the highest number of transcripts

maximumcell2 = pd.concat([maximumcell, x2],axis = 1)                          ### where X2 if a dataframe with the gene names
print(maximumcell2)

###tip
### How to see the whole row for a gene: in case you want to print a row and check yourself the number of counts (in case you did something wrong)
#with pd.option_context('display.max_rows', None):
    #print(dataset.loc['PLP1'])




##### Genes with low maximum number of transcripts (low_exp genes = 2877)
# From the genes with the highest number of transcripts there are many that have 1-5 transcripts only as their maximum.

low = maximumcell2[maximumcell2['Maximum number of transcripts'] <= 5]



####Find the frequency of those genes

frequency_low = low['Maximum number of transcripts'].value_counts()
low_exp = frequency_low.sort_index()
print(low_exp)
# Plot the results
plt.bar(low_exp.index, low_exp, facecolor = 'red')
plt.xlabel('Number of transcripts')
plt.ylabel('Gene frequency')
plt.title('Lowest number of transcripts in high expressed genes (gene drop outs)')
plt.show()                      #A high amount of genes (=1600) have one transcript as maximum expression level

##### Keep the highest expressed genes (18750 genes)
high_exp =  maximumcell2 [maximumcell2['Maximum number of transcripts'] > 5]
print(high_exp)


#####Box plot of highest expressed genes, grouped by cell type

gene_exp1.boxplot(column = 'total_number_of_genes', by ='cell_type')
plt.xlabel('Cell types')
plt.ylabel('Number of genes')
plt.title('Number of genes for each cell type')
plt.show()

##### Graph showing the 50 highest expressed genes (based on the maximum number of counts)

top50 = high_exp.nlargest(50, columns = 'Maximum number of transcripts')  #dataframe with the top 50 genes
print(top50)

fig = plt.figure(figsize = (8,8))
top = fig.add_subplot(1,1,1)

top_genes = top50['Genes']
y_pos = np.arange(len(top_genes))

#Horizontal bar plot
top.barh(y_pos, top50['Maximum number of transcripts'], color = 'green')

top.set_yticks(y_pos)
top.set_yticklabels(top_genes)
top.invert_yaxis()  # labels read top-to-bottom
top.set_xlabel('Number of transcripts')
top.set_title('Top 50 genes expression pattern')

plt.show()



#I will follow the process to find the total number of transcripts for each gene
# Total counts calculation

dataset['Total counts'] = dataset.sum(axis = 1)



# Find the total count-based top expressed genes

top50_genes = dataset.nlargest(50, columns = 'Total counts')  #dataframe with the top 50 genes
print(top50_genes)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)

top_genes2 = top50_genes.index.values
y_pos = np.arange(len(top_genes2))

#Horizontal bar plot
ax.barh(y_pos, top50_genes['Total counts'], color = 'green')

ax.set_yticks(y_pos)
ax.set_yticklabels(top_genes2)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Total number of counts')
ax.set_title('Top 50 genes expression pattern')

plt.show()



#remove the genes that have <= 5  transcripts from the original dataset
#use the low dataset as a template to remove the genes (turn the gene names into list first)

dataset2 = dataset.drop(low['Genes'].values.tolist(), axis =0)
print(dataset2)
