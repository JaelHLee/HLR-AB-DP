# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 23:50:02 2022

@author: Thinkpad
"""
#connectivity 计算指标：连通率
# 运行：
#     0.设置工作路径。
#     1.修改参数。90-96行。（参数说明见“readme2.md”）
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
import seaborn as sns #可视化的包
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# import matplotlib as mpl
#作图
def plotFIG3(ins):
        
        instring = ins + ".npy"
        data = np.load(instring,allow_pickle=True)


        yy = data[0] #阈值
        fin_con=data[1]

        

        plt.figure()	
        fig, ax1 = plt.subplots()
        ax1.set_xlabel(r'$y$', fontsize = 20)
        ax1.set_ylabel(r'$connectivity$', fontsize = 20)
        N=len(fin_con)
        # colors=[(r,g,b) for (r,g,b) in zip(np.linspace(0,1,N),np.linspace(1,0,N),np.linspace(0,1,N))]
        col=[(203/256,180/256,123/256),(91/256,183/256,205/256),(71/256,120/256,185/256),(84/256,172/256,117/256),(197/256,86/256,89/256),(153/256,50/256,204/256)]

        for i in range(len(fin_con)):
            a=i
            plt.plot(yy,fin_con[i],linestyle=':', linewidth=1,color=col[i])
            plt.plot(yy,fin_con[i],'o',markersize=3,label='Scenario=%d'% a,color=col[i])
        plt.legend(fontsize=12)
        ax1.set_ylim(bottom=0.)
        ax1.set_xlim(left=0.)
#         #子图
#         axins = inset_axes(ax1, width="40%", height="30%", loc='lower left',
#                     bbox_to_anchor=(0.2, 0.4, 1, 1),
#                     bbox_transform=ax1.transAxes)
#         for i in range(len(fin_con)):
#             axins.plot(yy,fin_con[i],'o',markersize=2,color=col[i])
            
#         zone_left = 5
#         zone_right = 20

# # 坐标轴的扩展比例（根据实际数据调整）
#         x_ratio = 0.01 # x轴显示范围的扩展比例
#         y_ratio = 0.01 # y轴显示范围的扩展比例

# # X轴的显示范围
#         xlim0 = yy[zone_left]-(yy[zone_right]-yy[zone_left])*x_ratio
#         xlim1 = yy[zone_right]+(yy[zone_right]-yy[zone_left])*x_ratio

# # Y轴的显示范围
#         y = np.hstack((fin_con[0][zone_left:zone_right], fin_con[1][zone_left:zone_right], fin_con[2][zone_left:zone_right], fin_con[3][zone_left:zone_right], fin_con[4][zone_left:zone_right]))
#         ylim0 = np.min(y)-(np.max(y)-np.min(y))*y_ratio
#         ylim1 = np.max(y)+(np.max(y)-np.min(y))*y_ratio

# # 调整子坐标系的显示范围
#         axins.set_xlim(xlim0, xlim1)
#         axins.set_ylim(ylim0, ylim1)

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
    yy = np.linspace(0.02,3,150) 

    if(check): 

        TOLLK = 1e-12 #误差
        
        M=1

        
        adaptive_deg = AdaptiveDegradation(ins,dist,outstring,plot,scenario,type,TOLLK) #创建类示例

        G,n_node,n_edge  = adaptive_deg.initialize_network() #按参数生成路径中的网络，并赋值权重
        if type != "no":
            G = adaptive_deg.create_network(G,type, n_node,n_edge) #生成指定type的网络
        weight=adaptive_deg.weight_seq(G) #权重序列
        fin_con=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
            fin_con.append(con)
        
            for i in range(1):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
               
                 
            
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con0.xlsx')
            
            
            
            for i in range(2):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con1.xlsx')               
                 
            
            for i in range(3):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con2.xlsx')     
            
            
            for i in range(4):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con3.xlsx')       
            
            
            for i in range(5):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con4.xlsx')                 

            for i in range(6):
                sce=scenario[i]
                con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
                fin_con.append(con)
                data=[]
                colunmns = ['sce','yy', 'fin_con']
                data = np.zeros((len(con), 3))
                for j in range(len(con)):
                    data[j, 0] = scenario[i] 
                    data[j, 1] = (yy)[j]
                    data[j, 2] = (con)[j]
                    data
            result = pd.DataFrame(columns = colunmns, data = data)
            result
            pd.DataFrame.to_excel(result, 'con5.xlsx')     


        fin_con=[]
        for i in range(len(scenario)):
            sce=scenario[i]
            con=adaptive_deg.caculate_connectivity(M,yy,G,sce)
            fin_con.append(con)       
        
        
        
        
        
        OUT = []
        OUT.append(yy)
        OUT.append(fin_con)

        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)

        if plot == 'plot':
            plotFIG3(outstring)

if __name__ == "__main__":
    main()

#run in console:
#run FIG3.py datasets/mouse_kasthuri_directed.txt uniform 5000 out_fig3 plot
