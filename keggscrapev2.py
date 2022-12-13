#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 20:08:48 2021

@author: hugokleikamp
"""


#%% clear variables
#
#try:
#    from IPython import get_ipython
#    get_ipython().magic('clear')
#    get_ipython().magic('reset -f')
#except:
#    pass

#%% change directory to script directory
import os
from pathlib import Path
import sys
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

#%% Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import random, re, requests
import threading, time, string
import pickle, heapq, math

from itertools import chain, groupby
from collections import Counter
from openpyxl import load_workbook 

import urllib
#%% download KO terms

print("Downloading Kegg KO file")
KEGGurl="https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir="
urllib.request.urlretrieve(KEGGurl, "ko00001.keg")

#read and parse
with open("ko00001.keg") as file:
    print("Parsing Kegg file")
    lines=file.readlines()

unnested=list()
columns=["A","B","C","D"]
for line in lines:
    if line.startswith("A"): 
        cat1=line[1:].strip() 
        unnested.append("/t".join([cat1]))
    if line.startswith("B"): 
        cat2=line[1:].strip()
        unnested.append("/t".join([cat1,cat2]))
    if line.startswith("C"): 
        cat3=line[1:].strip()
        unnested.append("/t".join([cat1,cat2,cat3]))
    if line.startswith("D"): 
        cat4=line[1:].strip()
        unnested.append("/t".join([cat1,cat2,cat3,cat4]))


df=pd.DataFrame(unnested,columns=["val"])
df=df["val"].str.rsplit("/t",expand=True)
df=df.iloc[:,: 4]
df.columns=["cat1","cat2","cat3","cat4"]
df=df[~df["cat4"].isnull()]

#split cat4
df[["KO","Genes"]] = df["cat4"].str.split(" ", 1, expand=True)
df["Description"]=df["Genes"].str.split("; ", 1, expand=True)[1]
df["Genes"]=df["Genes"].str.split("; ", 1, expand=True)[0]
df["ec"]=df["Description"].str.split("EC:", 1, expand=True)[1].str[:-1]
df["Description"]=df["Description"].str.split("EC:", 1, expand=True)[0].str[:-1]

#explode genes and EC nos for exact string matching
df["Genes"]=df["Genes"].str.split(",")
df=df.explode('Genes').reset_index(drop=True)

df["ec"]=df["ec"].str.split(" ")
df=df.explode("ec").reset_index(drop=True)

# strip spaces
for i in df.columns:
    df[i]=df[i].astype(str).str.strip()


df.to_excel("KO_pathways.xlsx")


cats=["09100 Metabolism",
"09120 Genetic Information Processing",
"09130 Environmental Information Processing",
"09140 Cellular Processes"]

metabolicKOs=list(set(df.loc[df["cat1"].isin(cats),"KO"]))

#%% get kegg accessions

def scrape_accs(r,url,i):
    while True:
        try:
            links=requests.get(url,stream=True).text
            s_links=links.split('<a href="/dbget-bin/www_bget?')
            acc=[a.split('">')[0] for a in s_links][1:]
            break
        except:
            print("sleeping")
            time.sleep(2)
      
    r.extend(list(zip(acc,[i]*len(acc)))) # add ko term
    
                        
accs=list()
counter=0
threads=[]
l_url='https://www.genome.jp/dbget-bin/www_bfind_sub?dbkey=genes&keywords='
r_url="&mode=bfind&max_hit=nolimit"


for i in metabolicKOs:
    
    counter+=1
    print(counter)
    url=l_url+i+r_url
    t=threading.Thread(target=scrape_accs, args=[accs,url,i])
    t.start()
    threads.append(t)

for thread in threads:
    thread.join()
    
#%% convert into ncbi acessions

accs = [tuple(l) for l in accs] #convert to tuple for set
accs=list(set(accs))

def chunks(lst,n):
    for i in range(0,len(lst),n):
        yield lst[i:i+n]
batches=[chunk for chunk in chunks(accs,80)] #kegg only gives you a hunderd lines at a time here

def convert_accs(batch,conv_accs):
    
    accs=[i[0] for i in batch]
    kos= [i[1] for i in batch]
    
    while True:
        
        url="http://rest.kegg.jp/conv/ncbi-proteinid/"+("+").join(accs)
        r=requests.get(url,stream=True).text

        if "permission" in r or "too large" in r: 
            print("caught!, sleeping")
            time.sleep(100)
            
        else:
            break
        
        
    #print(url)
    ncbi_accs=[i.split('ncbi-proteinid:')[1] for i in re.split('[\n \t]',r)[1:-1:2]]
    kegg_accs=                            [i for i in re.split('[\n \t]',r)[0:-1:2]]
    
    conv_accs.extend(list(zip(ncbi_accs,kegg_accs,kos)))

conv_accs=list()
counter=0
threads=[]
for batch in batches:
    
    counter+=1
    print(counter)
    time.sleep(random.uniform(0.2,0.25))

    t=threading.Thread(target=convert_accs, args=[batch,conv_accs])    
    t.start()
    threads.append(t)
    
    if counter%200==0:
        print("unwinding")
        for thread in threads:
            thread.join()
    #break
        

for thread in threads:
    thread.join()
        
        
#combine df
convdf=pd.DataFrame(conv_accs,columns=["ncbi_accs","kegg_accs","KOs"])
convdf.to_excel("Kegg2ncbi.xlsx")



#%% scrape ncbi for fastas
from Bio import Entrez
ncbiaccs=list(set(convdf["ncbi_accs"].to_list()))
batches=[chunk for chunk in chunks(ncbiaccs,1000)] 

def ncbi_fastas(batch,f):
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text")
    f.write(handle.read())

with open("kegg_ncbi.fasta","a+") as f:
    counter=0
    threads=[]
    for batch in batches:
        counter+=1
        print(counter)
        time.sleep(random.uniform(0.2,0.25))
        t=threading.Thread(target=ncbi_fastas, args=[batch,f])    
        t.start()
        threads.append(t)
        
        if counter%50==0:
            print("unwinding")
            for thread in threads:
                thread.join()
    
    for thread in threads:
        thread.join()

