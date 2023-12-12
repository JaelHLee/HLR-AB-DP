# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 01:11:47 2022

@author: Thinkpad
"""
#meank 计算指标：平均度值和平均权重
# 运行：
#     0.设置工作路径。
#     1.修改参数。73-79行。（参数说明见“readme2.md”）
# ````对于一个数据文件，每个指标都需要运行4次。````
# ````
#     2.点击run即可出图。
#     3.图片保存，以网络类型命名,存放至各自指标文件夹内。

import matplotlib.pyplot as plt
from sys import argv
import numpy as np
import pandas as pd
from deg_class import AdaptiveDegradation
# from get_input_params import input_params
# import seaborn as sns

def plotFIG4(ins):

    instring = ins + ".npy"
    data = np.load(instring,allow_pickle=True)

    fin_meansk = data[0]
    fin_meansw = data[1]
    yy = data[2]


    
    # sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})

    plt.figure()	
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$y$', fontsize = 20)
    ax1.set_ylabel(r'$\bar{k}(y)$', fontsize = 20)
    col=[(203/256,180/256,123/256),(91/256,183/256,205/256),(71/256,120/256,185/256),(84/256,172/256,117/256),(197/256,86/256,89/256),(153/256,50/256,204/256)]
    for i in range(len(fin_meansk)):
        a=i
        plt.plot(yy,fin_meansk[i],'o',markersize=2.5,label='Scenario=%d'%a,color=col[i]) #平均度值-红色
        plt.plot(yy, fin_meansk[i],'--',linewidth=0.5,color=col[i])
    ax1.set_ylim(bottom=0.)
    ax1.set_xlim(left=0.)
    plt.legend(frameon=True,fontsize=12)
    plt.show()
    plt.close()
    
    fig, ax2 = plt.subplots()
    ax2.set_xlabel(r'$y$', fontsize = 20)
    ax2.set_ylabel(r'$\bar{w}(y)$', fontsize = 20)  
    for i in range(len(fin_meansw)):
        a=i
        plt.plot(yy,fin_meansw[i],'^',markersize=0.5,label='Scenario=%d'% a,color=col[i]) #平均度值-红色
        plt.plot(yy, fin_meansw[i],'--',linewidth=0.5,color=col[i])
    ax2.set_ylim(bottom=0.)
    ax2.set_xlim(left=0.)
    plt.legend(frameon=True,fontsize=9)
    plt.show()
    plt.close()






def main():

            #设置参数
    check=True
    ins='datasets/data/Co_2019.txt'
    dist='no'
    outstring='out_connectivity'
    plot='plot'
    scenario=[0,1,2,3,4,5]
    type='real'
    yy = np.linspace(0.02,3,150)

    if(check): 
        TOLLK = 1e-12 #偏差
         #阈值（50步）
        #创建类实例
        
        adaptive_deg = AdaptiveDegradation(ins,dist,outstring,plot,scenario,type,TOLLK) #创建类示例
        G,n_node,n_edge  = adaptive_deg.initialize_network() #按参数生成路径中的网络，并赋值权重
        if type != "no":
            np.random.seed(1)
            G = adaptive_deg.create_network(G,type, n_node,n_edge) #生成指定type的网络
        weight=adaptive_deg.weight_seq(G) #权重序列
        SAMPLES=len(weight)+1


        fin_meansk=[]
        fin_meansw=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
            fin_meansk.append(MEANSK)
            fin_meansw.append(MEANSW)

            for i in range(1):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 0.xlsx')

            for i in range(2):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 1.xlsx')

            for i in range(3):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 2.xlsx')
            
            for i in range(4):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 3.xlsx')      
            
            for i in range(5):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 4.xlsx')            

            for i in range(6):
                sce=scenario[i]
                MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
                fin_meansk.append(MEANSK)
                fin_meansw.append(MEANSW)
                data=[]
                colunmns = ['sce','yy', 'fin_meansk','fin_meansw']
                data = np.zeros((len(MEANSK), 4))
                for j in range(len(MEANSK)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (MEANSK)[j]
                    data[j, 3] = (MEANSW)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'meansk&w 5.xlsx')


        fin_meansk=[]
        fin_meansw=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            MEANSW,MEANSK = adaptive_deg.caculate_meanskw(yy,G,sce)
            fin_meansk.append(MEANSK)
            fin_meansw.append(MEANSW)
        
        OUT = []
        OUT.append(fin_meansk)#平均权重
        OUT.append(fin_meansw)#平均度值
        OUT.append(yy) #阈值


        
        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)

        if plot == 'plot':
            plotFIG4(outstring)

if __name__ == "__main__":
    main()
    
#run FIG4.py datasets/RG_N5000_k4.txt uniform 5000 out_fig4 plot