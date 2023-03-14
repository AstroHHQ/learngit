import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
hvlist = np.loadtxt('hvlist.txt')
f = open(f'namelist.txt', "r", encoding="utf-8")
str1 = f.read()
namelist = str1.split()
dirlist = ['2_2','2_1']
clist = [['b','lightsteelblue'], \
            ['brown','lightcoral'], \
            ['y','khaki'], \
            ['seagreen','lightgreen'], \
            ['m','violet']]
print(namelist,dirlist)
clist[2][1]

md = 400
plt.xlabel('D')
plt.ylabel('D_wise')

def getratio(df):
    return (df['D']-df['D_wise'])*100/df['D_wise']

for di in range(len(dirlist)): 
    cl = (clist[di][0])
    ecl = (clist[di][1])
    mcmcdir = 'mcmc'+dirlist[di]
    df0 = pd.read_excel(f'./ans_excel/{mcmcdir}.xlsx')
    df = df0[['name','eta','pv_wise','pv','D','Dup','Ddown','D_wise','D_wiseErr']]
    df = df.sort_values('D',ascending = False)
    df1 = df
    for ii in range(len(df)):
        if df.loc[ii,'D_wise']>500 or df.loc[ii,'D_wise']<10 or df.loc[ii,'D']>500:
            df1 = df1.drop(labels=ii)
    #df1 = df1.drop(labels=28)
    df1.loc[:,'delta/WISE'] = df1.apply(getratio,axis=1)
    md = max(md, df1['D'].max())
    x = df1['D']
    y = df1['D_wise']
    #plt.plot(x,y,marker = '*')
    
    Wyerr = df1['D_wiseErr']
    xerr1 = (df1['D']-df1['Ddown'])
    xerr2 = (df1['Dup']-df1['D'])
    plt.errorbar(x,y,xerr=[xerr1,xerr2],yerr=Wyerr,marker = '.',linestyle="none", fmt="o",color=cl,ecolor=ecl, capsize=1.0,label = mcmcdir)
    #plt.errorbar(x,y,xerr=[xerr1,xerr2],linestyle="none")
plt.legend()    
plt.plot(np.linspace(0,md+50),np.linspace(0,md+50),'navajowhite')
ax = plt.gca()
ax.set_aspect(1)
plt.savefig(f'./ans_eps/fig_'+dirlist[0]+'_'+dirlist[1]+'.eps')