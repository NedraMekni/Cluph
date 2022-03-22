import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.ndimage
import math
import matplotlib
import random
import sys
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.cluster import AgglomerativeClustering
from sklearn import datasets, cluster
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.impute import SimpleImputer
from scipy.cluster.hierarchy import dendrogram, linkage
from tqdm import tqdm
from scipy.stats import norm
from kneed import DataGenerator, KneeLocator
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from math import isnan
from datatest import validate

dimension = [2**i for i in range(1,11)]
transf = lambda i: (i+1)*0.5
colors = []

dict_names = {}
title = []
cmp_name = []
raw_value = []

fname = 'HM_molport_MD_test_.csv'

def rand_col(min_v=0.5, max_v=1.0):
  hsv = np.concatenate([np.random.rand(2), np.random.uniform(min_v, max_v, size=1)])
  return matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb(hsv))

def parse_csv(name):
  global dimension,raw_value,cmp_name
  i = 0
  with open(fname,'r') as f:
    title = f.readline().strip().split(',')[1:]
    dimension+=[len(title)]
    for l in f:
      cmp_name += [l.strip().split(',')[0][0:-1]]
      raw_value += [[float(x) if len(x)>0 else -1.7976931348623157e-308 for x in l.strip().split(',')[1:]]]
    '''
    for l in f:
      dict_names[i] = l.strip().split(",")[0] 
      res+=[l.strip().split(",")[1:]]
      i+=1
    '''

def imputing_missing_values(data):
    return SimpleImputer(missing_values=np.NaN, strategy='mean')

def labeling_results(res,plt,labels):
  global dict_names
  for i, txt in enumerate(res):    
    print('{} => {}'.format(i,dict_names[i]))
    plt.annotate(i, (txt[0],txt[1]))
  return plt

def input_lab(crds):
  print("provide labels in CSV format(A for all): ",end='')
  l = input()
  return [x for x in range(len(crds))] if l.upper()=='A' else [int(x.strip()) for x in l.split(',')]

def evaluate_in_components(dsc_list):
  res = []
  #print(dsc_list)
  for n_comp in tqdm(dimension):
    pca = PCA()
    #print('n comp = ',n_comp)
    print(n_comp)
    crds = pca.fit_transform(preprocessing.scale([x[:n_comp] for x in dsc_list]))
    var = np.sum(pca.explained_variance_ratio_)
    print(var)
    res.append([n_comp,var])
  return res   

def evaluate_out_components(dsc_list):
  res = []
  #print(dsc_list)
  for n_comp in tqdm(range(2,50)):
    pca = PCA(n_components=n_comp)
    crds = pca.fit_transform(preprocessing.scale(dsc_list))
    var = np.sum(pca.explained_variance_ratio_)
    res.append([n_comp,var])
  return res   

def inflection_point2(xs, ys):
  kneedle = KneeLocator(xs, ys, S=1.0, curve='concave', direction='increasing')
  return kneedle.elbow

def runPCA(dsc_list, n_comp):
  pca = PCA(n_components=n_comp)
  crds = pca.fit_transform(dsc_list)
  return crds

# clst_type means cluster type e.g. (0)DBSCAN or (1)Hierarchical
def do_clustering(clst_type, data):
  if clst_type == 1:
    # cluster type : DBSCAN
    return DBSCAN(min_samples=2,eps=8.72).fit(data)

  if clst_type == 0:
    # cluster type : Hierarchical
    dend = dendrogram(linkage(data,method='ward'))
    plt.show()
    return AgglomerativeClustering(n_clusters=int(input('Number of cluster obsereved: '))).fit(data)
  print("* clst_type not set properly")
  exit()  

def clusters_evaluation(clustering,data):
  cluster_labels= clustering.fit_predict(data)
  silhouette_avg = silhouette_score(data, cluster_labels)
  return silhouette_avg


if __name__ == "__main__":
  
  title = []
  cmp_name = []
  raw_value = []
  
  parse_csv(fname)
  dim = dimension[-1]
  print(type(raw_value))

  np_raw_value = np.array(raw_value)
  imp = imputing_missing_values(np_raw_value)
  scaler = StandardScaler()
  
  print("Raw data")
  print('Number of compounds',len(np_raw_value))
  print("Number of molecular descriptors = {}".format(dim))
  
  df = pd.DataFrame(np_raw_value)
  #print(df)
  
  #Deal with inf and nan values
  t = df.drop(df.columns[np.isinf(df).any()], axis=1)
  w = df.replace([np.inf, -np.inf], np.nan, inplace =True)
  z = df.dropna(inplace = True)
  
  #Standardize features by removing the mean and scaling to unit variance
  df = pd.DataFrame(scaler.fit_transform(df),columns = df.columns)
  raw_value = df
  corr_matrix =df.corr().abs()
  upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
  plt.figure()
  upper_filt_hmap=sns.heatmap(upper)
  plt.show()

  #Remove correlated features
  to_drop = [column for column in upper.columns if any(upper[column]>0.75)]
  raw_value_clean = df.drop(df[to_drop],axis = 1)
  raw_value = raw_value_clean
  
  print('Features to drop',len(to_drop))
  #print('Cleaned data',len(raw_value))
  ft = pd.DataFrame(raw_value)
  print(ft)
  
  #DataFrame to numpy
  raw_value = df.to_numpy()
  print(type(np_raw_value))
  np_raw_value = scaler.fit_transform(np_raw_value)

  
  #Explained variance and inflection point of eigenvalues
  res_out = evaluate_out_components(raw_value)
  ip_out = [inflection_point2([x[0] for x in res_out],[x[1] for x in res_out])][0]
  print('Inflection point estimation for dimension {}: {}'.format(-1,ip_out))
  print(res_out)
  raw_value=np.array([x[:dim] for x in raw_value])
  
  print(type(raw_value))
  
  raw_value=[x for x in raw_value]
  
  #Perfomrm data reduction with PCA
  #print('len data = {}, len data[0] = {}'.format(len(raw_value),len(raw_value[0])))
  pca = PCA().fit(preprocessing.scale(raw_value))

  i=0
 
  for eigenvalue,eigenvector in zip(pca.explained_variance_,pca.components_):
    print('eigenvector = {}, eigenvalue = {}, ratio = {}'.format(eigenvector,eigenvalue,pca.explained_variance_ratio_[i]))
    i+=1
  i=0
  cov_mat=np.cov(preprocessing.scale(raw_value), rowvar=False)
  #print('len cov_mat = {}, len cov_mat[0]={}'.format(len(cov_mat),len(cov_mat[0])))
  print('components_ ',pca.components_)
  print('len components_ = {}, len components_[0] = {}'.format(len(pca.components_),len(pca.components_[0])))
  plt.rcParams["figure.figsize"] = (12,6)

  fig, ax = plt.subplots()
  xi = np.arange(1, dim+1, step=1)[:140]
  y = np.cumsum(pca.explained_variance_ratio_)[:140]

  print(len(xi))
  print(y)
  print(pca.explained_variance_ratio_)
  plt.ylim(0.0,1.0)
  plt.plot(xi, y, marker='o', linestyle='--', color='b')

  plt.xlabel('Number of Components')
  plt.xticks(np.arange(0, 150, step=4)) #change from 0-based array index to 1-based human-readable label
  plt.ylabel('Cumulative variance (%)')
  plt.title('The number of components needed to explain variance')

  plt.axhline(y=0.85, color='r', linestyle='-')
  plt.axhline(y=0.75, color='b', linestyle='-')
  plt.axhline(y=0.65, color='y', linestyle='-')
  plt.text(0.5, 0.85, '85% cut-off threshold', color = 'red', fontsize=16)
  plt.text(0.4, 0.75, '75% cut-off threshold', color = 'b', fontsize=16)
  plt.text(0.3, 0.65, '65% cut-off threshold', color = 'y', fontsize=16)
  ax.grid(axis='x')
  plt.show()

  fig = plt.figure(figsize=(8,6))
  ax = fig.add_subplot(111)
  print('#####',pca.components_[0], pca.components_)
  
  data=preprocessing.scale(raw_value)
  data_out= runPCA(data,int(input('Number of components: ')))
  
  #Perform Hierachichal clustering 
  
  clst_type = 0
  clustering = do_clustering(clst_type,data_out)
  silhouette_score = clusters_evaluation(clustering,data_out)
  
  col = 1
  
  for label in set(clustering.labels_):
    print(label)
    points = [data_out[i] for i in range(len(data_out)) if clustering.labels_[i] == label]
    indexes = [i for i in range(len(data_out)) if clustering.labels_[i] == label]
    #plt.scatter([el[0] for el in points],[el[1] for el in points],label=label,c=colors[random.randint(0,len(colors)-1)])
    plt.scatter([el[0] for el in points],[el[1] for el in points],label=label,c=rand_col())
    
    plt.legend()
    print('{} => {}'.format(label,[cmp_name[i] for i in indexes ]))
  
 