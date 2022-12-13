#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:37:32 2022

@author: hugokleikamp
"""




#%% clear variables

# try:
#     from IPython import get_ipython
#     get_ipython().magic('clear')
#     get_ipython().magic('reset -f')
# except:
#     pass

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())

#%%
import urllib
import datetime
import pandas as pd
import time
import requests
import threading
#%% KEGG KO file

#download
print("Downloading Kegg file")
KEGGurl="https://www.genome.jp/kegg-bin/download_htext?htext=ko00002&format=htext&filedir="
urllib.request.urlretrieve(KEGGurl, "ko00002.keg")

#read and parse
with open("ko00002.keg") as file:
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

#%%
df=pd.DataFrame(unnested,columns=["val"])

df=df["val"].str.replace("<b>","").str.replace("</b>","").str.rsplit("/t",expand=True)
df=df.iloc[:,: 4]
df.columns=["cat1","cat2","cat3","cat4"]
df=df[~df["cat4"].isnull()]


df["Pathway ID"]=df["cat4"].str.split(" ",1).apply(lambda x: x[0])



def scrape_accs(res,id):
    
    url="https://www.genome.jp/kegg-bin/module_ko_list?map="+id+"&org=ko"
    
    while True:
        try:
            # links=requests.get(url,stream=True).text
            # s_links=links.split('<a href="/dbget-bin/www_bget?')
            # acc=[a.split('">')[0] for a in s_links][1:]
            
            r=requests.get(url,stream=True).text
            lines=[i for i in r.splitlines() if i.startswith("    <td>")]
            df=pd.DataFrame(list(zip(
                    lines[0::2],
                    lines[1::2])),
                columns=["KO","Unparsed_Description"])
            df["KO"]=df["KO"].str[24:30]
            df=df[df["KO"].str.startswith("K")] #only orthologies, no compounds

            s=df["Unparsed_Description"].str.rsplit(">",expand=True)
            df["ECs"]=s.iloc[:,2::2].apply(lambda x: " ".join(set(x.dropna())),axis=1).str.replace("</a","",regex=False).str.strip()
            df[["Genes","Description"]]=s.iloc[:,1].str.split("<td>").apply(lambda x: x[0]).str.split("[EC:",regex=False).apply(lambda x: x[0]).str.replace("</td","",regex=False).str.strip().str.rsplit("; ",expand=True)
            df=df[["KO","Genes","Description","ECs"]]
            df["Pathway ID"]=id
            res.append(df)
            
            break
        except:
            print("sleeping")
            time.sleep(20)
      

    
                        
res=list()
counter=0
threads=[]

base_thread=threading.active_count()

for i in df["Pathway ID"]:
    
    counter+=1
    print(counter)
    t=threading.Thread(target=scrape_accs, args=[res,i])
    t.start()
    threads.append(t)
    
    #unwind in case of thread overload, manage server traffic
    cur_thread=threading.active_count()
    if (cur_thread-base_thread)>200:
        print("unwinding, query at: "+str(counter/len(df)))
        for thread in threads:
            thread.join()
        threads=[] #this seems to act different on windows?

for thread in threads:
    thread.join()

#%%

c=df.merge(pd.concat(res),on="Pathway ID",how="left")
c.to_csv(str(datetime.datetime.today()).split()[0]+"_modules.tsv",sep="\t")
