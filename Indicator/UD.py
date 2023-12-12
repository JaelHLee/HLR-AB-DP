# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 23:01:33 2022

@author: Thinkpad
"""
#UD 指标：未受影响的贸易量
# 运行：
#     0.设置工作路径。
#     1.修改参数。50-56行。（参数说明见“readme2.md”）
# ````对于一个数据文件，每个指标都需要运行4次。````
# ````
#     2.点击run即可出图。
#     3.图片保存，以网络类型命名,存放至各自指标文件夹内。

import matplotlib.pyplot as plt
from sys import argv
import numpy as np
import pandas as pd
from deg_class import AdaptiveDegradation


def plotFIG2(ins):
        #读取变量文件
        instring = ins + ".npy"
        data = np.load(instring,allow_pickle=True)
        #赋值给新的变量
        yy = data[0]
        fin_UD=data[1]
        
        plt.figure()	
        fig, ax1 = plt.subplots()
        ax1.set_xlabel(r'$y$', fontsize = 20)
        ax1.set_ylabel(r'$UD$', fontsize = 20)
        col=[(203/256,180/256,123/256),(91/256,183/256,205/256),(71/256,120/256,185/256),(84/256,172/256,117/256),(197/256,86/256,89/256),(153/256,50/256,204/256)]

        for i in range(len(fin_UD)):
            a=i
            plt.plot(yy,fin_UD[i],'o',markersize=3,label='Scenario=%d'% a,color=col[i]) #平均度值-红色
            plt.plot(yy,fin_UD[i],'--',color=col[i],linewidth=0.5)
        ax1.set_ylim(bottom=0.)
        ax1.set_xlim(left=0.)
        plt.legend(fontsize=12)
        plt.show()



def main():
    #设置参数
    check=True
    ins='datasets/data/Co_2019.txt'
    dist='no'
    outstring='out_connectivity'
    plot='plot'
    scenario=[0,1,2,3,4,5]
    type='real'
    yy=np.linspace(0.02,3,150)

    if(check): 
        
        TOLLK = 1e-12 #偏差
        

        adaptive_deg = AdaptiveDegradation(ins,dist,outstring,plot,scenario,type,TOLLK) #创建类示例

        G,n_node,n_edge  = adaptive_deg.initialize_network() #按参数生成路径中的网络，并赋值权重
        if type != "no":
            G = adaptive_deg.create_network(G,type, n_node,n_edge) #生成指定type的网络
        weight=adaptive_deg.weight_seq(G) #权重序列
        SAMPLES=len(weight)+1
        
        fin_UD=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
            fin_UD.append(UD)
        
        
        
            for i in range(1):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD0.xlsx')
            
            
            
            for i in range(2):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD1.xlsx')  
            
            
            for i in range(3):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD2.xlsx')


            for i in range(4):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD3.xlsx')     
            
            
            for i in range(5):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD4.xlsx')  
            
            
            for i in range(6):
                sce=scenario[i]
                UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
                fin_UD.append(UD)
                data=[]
                colunmns = ['sce','yy', 'fin_UD']
                data = np.zeros((len(UD), 3))
                for j in range(len(UD)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (UD)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'UD5.xlsx')

        fin_UD=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            UD = adaptive_deg.caculate_ud(yy,G,TOLLK,sce)
            fin_UD.append(UD)             
            
        #将以上结果存储在OUT变量中，并保存在本地文件夹中。若plot，则调用plotFIG2函数，并将OUT变量作为outstring传递给函数。
        OUT = []
        OUT.append(yy)
        OUT.append(fin_UD)

        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)       


        if plot == 'plot' :
            plotFIG2(outstring)


if __name__ == "__main__":
    main()
    
#run in console:
#run FIG2.py datasets/RG_N5000_k4.txt gauss 5000 out_fig2 plot
#run FIG2.py datasets/data/a2021_2610.txt uniform out_fig2 plot 1 random