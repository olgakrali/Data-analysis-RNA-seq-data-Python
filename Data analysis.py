import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler #for standardization
from sklearn import preprocessing                #for normalization
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE               #for tSNE
from sklearn.cluster import KMeans #for K-means clustering

# working on the genes that are expressed in at least on cell (genes with zero expression have already removed)
path = "C:/Users/"

dataset = pd.read_csv(path + "mydata.csv", sep = ",")
labels = pd.read_csv(path + "labels2.csv", sep = ",")

# change the name of the first column to genes
dataset = dataset.rename(columns = {'Unnamed: 0':'genes'})

# Row names index is equal to numbers. Set the 1st column as rownames

dataset.index = dataset['genes']


# delete first column
del dataset['genes']



print(dataset)
print(labels)
#work on transposed data instead
dataset = dataset.transpose()

           #Principal component analysis
            #data normalization type L2
datanorm = preprocessing.normalize(dataset, norm = 'l2')

#Dataset data
#counts = dataset.columns.values #gives us all columns  (Remove comments if you want to work on standardized data)

#Labels data
label = labels[['cell_type']] #or labels[labels.columns[3]]

# Separating out the genes
# x = dataset.loc[:,counts].values
            # Standardizing the cells
#x = StandardScaler().fit_transform(x)

#PCA on non stardandised data. (Only normalization)
pca = PCA(n_components = 2)
principalComponents = pca.fit_transform(datanorm)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, labels[['cell_type']]], axis = 1)

print(finalDf)

#Create PCA plot

fig = plt.figure(figsize = (8,8))

ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('PCA', fontsize = 20)
cells = ['oligodendrocytes', 'neurons', 'hybrid', 'astrocytes', 'endothelial', 'fetal_quiescent', 'fetal_replicating', 'OPC', 'microglia']
colors = ['r', 'b', 'g', 'teal', 'gray', 'palevioletred', 'm', 'peru', 'orange']
for cell, color in zip(cells,colors):
    indicesToKeep = finalDf['cell_type'] == cell
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(cells)
ax.grid()

plt.show()


print(pca.explained_variance_ratio_)    #1st component explains 14% of the variation and 8.5% the second




#Reduction of multidimensionality through t-distributed stochastic neighbor embedded (tSNE)
               # Use the normalized data datanorm #

                            ### set random seed ###
np.random.seed(123456)

tsne = TSNE(n_components = 2, perplexity = 40, n_iter = 500)
tsnecompon = tsne.fit_transform(datanorm)

tsnedf = pd.DataFrame(data = tsnecompon, columns = ['Dimension 1', 'Dimension 2'])

finaldf = pd.concat([tsnedf, labels[['cell_type']]], axis = 1)


#t-SNE plot


fig = plt.figure(figsize = (8,8))
plot = fig.add_subplot(1,1,1)
plot.set_xlabel('Dimension 1', fontsize = 15)
plot.set_ylabel('Dimension 2', fontsize = 15)
plot.set_title('tSNE plot', fontsize = 20)
celltypes = ['oligodendrocytes', 'neurons', 'hybrid', 'astrocytes', 'endothelial', 'fetal_quiescent', 'fetal_replicating', 'OPC', 'microglia']
colours = ['r', 'b', 'g', 'teal', 'gray', 'palevioletred', 'm', 'peru', 'orange']
for celltype, colour in zip(celltypes, colours):
    indicesToKeep = finaldf['cell_type'] == celltype
    plot.scatter(finaldf.loc[indicesToKeep, 'Dimension 1']
               , finaldf.loc[indicesToKeep, 'Dimension 2']
               , c = colour           #remember to change the name
               , s = 50)
plot.legend(celltypes)
plot.grid()

plt.show()


                            ### Clustering ###

#k means clustering, using datanorm

#Finding the optimum number of clusters for k-means clustering

wcss = []

for i in range(1, 15):
    kmeans = KMeans(n_clusters = i, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
    kmeans.fit(datanorm)
    wcss.append(kmeans.inertia_)

#Plotting the results onto a line graph, allowing us to observe 'The elbow'
plt.plot(range(1, 15), wcss)
plt.title('The elbow method')         #It gave me 8 clusters as the best solution (9 the original ones)
plt.xlabel('Clusters')
plt.ylabel('WCSS') #within cluster sum of squares
plt.show()


                        ### k means clustering with tSNE visualisation



kmeans = KMeans(n_clusters = 9).fit(datanorm)     # we can try the same with 8 clusters
lab = kmeans.labels_

clusters = pd.DataFrame(data = lab, columns = ['clusters'])


np.random.seed(123456)

tsnekmeans = TSNE(n_components = 2, perplexity = 40, n_iter = 500)
tsnek = tsnekmeans.fit_transform(datanorm)
tsnekdf = pd.DataFrame(data = tsnek, columns = ['Dimension 1', 'Dimension 2'])
finalkdf = pd.concat([tsnekdf,labels[['cell_type']],clusters], axis = 1)
print(finalkdf)         #not all the cells of the same cell type are grouped together

###t-SNE plot

fig = plt.figure(figsize = (8,8))
kplot = fig.add_subplot(1,1,1)
kplot.set_xlabel('Dimension 1', fontsize = 15)
kplot.set_ylabel('Dimension 2', fontsize = 15)
kplot.set_title('K-means clustering tSNE plot', fontsize = 20)
clusterings = [0, 1, 2, 3, 4, 5, 6, 7, 8]
colorings = ['r', 'b', 'g', 'teal', 'gray', 'palevioletred', 'm', 'peru', 'orange']
for clustering, coloring in zip(clusterings, colorings):
    indicesToKeep = finalkdf['clusters'] == clustering
    kplot.scatter(finalkdf.loc[indicesToKeep, 'Dimension 1']
               , finalkdf.loc[indicesToKeep, 'Dimension 2']
               , c = coloring           #remember to change the name
               , s = 50)
kplot.legend(clusterings)
kplot.grid()

plt.show()



