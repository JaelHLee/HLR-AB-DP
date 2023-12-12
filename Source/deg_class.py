from deg_analytic_func import *
from deg_numeric_func import *
import numpy as np
from progress.bar import Bar
import time
from numpy import linalg as LA
from random import sample
import networkx as nx
class AdaptiveDegradation:
    def __init__(self,ins,dist,outstring,plot,scenario,type,TOLLK):
        self.instring = ins
        self.weight_distribution_type = dist
        self.plot_result = plot
        self.tolerance = TOLLK
        self.sce = scenario
        self.type=type
    #生成路径中的网络，并赋值权重
    def initialize_network(self):
        G = nx.read_edgelist(self.instring,create_using=nx.DiGraph()) #创建空有向网络，并将路径中的边赋值给该网络
        assign_weights(G,self.weight_distribution_type,type) #生成权重，并赋值给网络
        n_node=G.number_of_nodes()
        n_edge=G.number_of_edges()

        return G, n_node,n_edge
    #创建权重数组。从0开始，最大入强度-dw结束，步长为dw的数组(横坐标的细分数)
    def initialize_w_ax(self,G,SAMPLES):
        w_max = max([in_strength(G,n) for n in G.nodes()])#获取最大的节点入强度
        dw = w_max / SAMPLES #步长 （最大入强度/间隔数）
        ww = np.arange(0,w_max,dw) #创建从0开始，(最大入强度-dw)结束，步长为dw的数组
        if(self.weight_distribution_type == 'no'):
            ww=np.delete(ww, 0) #删除第一个数0
        ww=ww/max(ww)
        return ww
    #创建出强度序列
    def initialize_ow_ax(self,G,SAMPLES):
        ow_max = max([out_strength(G,n) for n in G.nodes()])#获取最大的节点出强度
        dw = ow_max / SAMPLES #步长 （最大入强度/间隔数）
        ow = np.arange(0,ow_max,dw) #创建从0开始，(最大入强度-dw)结束，步长为dw的数组
        if(self.weight_distribution_type == 'no'):
            ow=np.delete(ow, 0) #删除第一个数0
        return ow
    #创建从0开始，最大入度结束，步长为1的数组
    def initialize_k_ax(self,G):

        k_max = max_in_deg(G) + 1 #网络G最大的入度+1
        kk = np.arange(0,k_max,1)
        
        return kk
    
    #权重序列
    def weight_seq(self,G):
        weight = []
        for e in G.edges():
            weight.append(G[e[0]][e[1]]['weight'])
        return weight

    #网络G的入度概率密度分布（范围：0-最大度，间隔为1）
    def in_degree_distribution(self,G,kk):

        kin_list = list([G.in_degree(n) for n in G.nodes()]) #入度序列
        pk0in = [kin_list.count(i)/len(kin_list) for i in kk] #每个入度值的概率，即入度概率密度分布
        
        return pk0in
    #网络G的出度概率密度分布
    def out_degree_distribution(self,G):
        kout_max = max_out_deg(G) + 1
        kk = np.arange(0,kout_max,1)
        k_list = list([G.out_degree(n) for n in G.nodes()])
        pk0out = [k_list.count(i)/len(k_list) for i in kk]
        
        return pk0out
    #生成指定的权重概率密度分布权重的基础上直接求概率密度函数
    def weight_distribution(self,G,ww,weight):
        dw = ww[1] - ww[0]
        if(self.weight_distribution_type == 'uniform'):
	        pw = np.heaviside(ww,1)*np.heaviside(1-ww,1) #生成均匀分布(权重在0-1之间为1，其余为0).结果全都是1
        if(self.weight_distribution_type == 'gauss'):
	        pw = norm(gaussian(ww,0.5,0.1),dw)#对权重的正态概率密度函数,并标准化
        if(self.weight_distribution_type == 'no'):
            pw=norm(pl_knn(weight),dw)#真实权重的概率密度函数
        return pw
    
        #生成指定的权重概率密度分布权重的基础上直接求概率密度函数
        dw = ow[1] - ow[0]
        if(self.weight_distribution_type == 'uniform'):
	        pow1 = np.heaviside(ow,1)*np.heaviside(1-ow,1) #生成均匀分布(权重在0-1之间为1，其余为0).结果全都是1
        if(self.weight_distribution_type == 'gauss'):
	        pow1 = norm(gaussian(ow,0.5,0.1),dw)#对权重的正态概率密度函数,并标准化
        if(self.weight_distribution_type == 'no'):
            pow1=norm(pl_knn(weight),dw)#真实权重的概率密度函数
        return pow1
    #权重-度联合概率密度分布 fk(x)=lk*wk(x)
    #输入：入度概率分布，权重概率密度分布，入度序列
    def weight_degree_distribution(self,pk,pw,kk):
        pwk = []
        pkex = PFK_EX(pk,kk) #余度分布lk
        for k in range(1,len(kk)):
	        pwk.append(pw*pkex[k-1])
        pwk = np.array(pwk)

        return pwk
    #权重-度联合概率密度分布 fk(x)=lk*wk(x)
    #输入：入度概率分布，权重概率密度分布，入度序列
    def outweight_degree_distribution(self,G,pk,pw):
        powk = []
        kout_max = max_out_deg(G) + 1
        kk = np.arange(0,kout_max,1)
        pkex = PFK_EX(pk,kk) #余度分布lk
        for k in range(1,len(kk)):
	        powk.append(pw*pkex[k-1])
        powk = np.array(powk)

        return powk
    #FIG2上：输入：阈值范围，网络，0-入强度的数组，权重分布.输出：预测的权重，经过无响应渗流后的权重序列
    def degradation_no_adaptive(self,yy,G,ww,pw):
        dw = ww[1] - ww[0]#步长
        G_noad = G.copy()#复制网络G
        PW_noad_T = []
        WW_noad = []
        
        with Bar('NO ADAPTIVE',max=len(yy)) as bar:
            for y in yy:
                pw = pw*np.heaviside(ww-y,1) #权重均匀分布序列，若大于阈值则为原权重，若小于阈值则为0。
                pw = norm(pw,dw) #标准化
                PW_noad_T.append(pw) #更新后的预测权重分布
                
                if self.sce==1:
                    remove_thresh(G_noad,y,0)#无响应渗流，没有重分配
                    print('1')
                elif self.sce == 2:
                    remove_thresh2(G_noad,y,0) #情景2
                    print('2')
                elif self.sce == 3:
                    remove_thresh3(G_noad,y,0) #情景3
                elif self.sce == 4:
                    remove_thresh4(G_noad,y,0) #情景4
                else:
                    remove_thresh5(G_noad,y,0) #情景5

                w_list = [G_noad[e[0]][e[1]]['weight'] for e in G_noad.edges()] #渗流后的网络权重
                WW_noad.append(w_list) #更新后网络的权重

                bar.next()        
        return PW_noad_T,WW_noad
    #FIG3：没有响应的过滤后-计算WGCC的大小
    def degradation_no_adaptive_WGCC(self,M,yy,G,ww,pw,pk0in,TOLLK):
        dw = ww[1] - ww[0] #步长
        G_noad = G.copy() #复制网络
        N = nx.number_of_nodes(G) #节点数
        WGCC_noad = []
        WGCC_noad_T = []
			#利用入度-出度的外积计算WGCC
        with Bar('NO ADAPTIVE - MODEL:', max=len(yy)) as bar:
            t1 = time.time()
            for y in yy:
                pkin_noad = PFK_NOAD(pk0in,y,pw,ww,dw,TOLLK)
                pkout_noad  = pkin_noad
                WGCC_noad_T.append(fixed_point_directed(np.outer(pkout_noad,pkin_noad)))#入度-出度的外积
                bar.next()
            Dt1 = time.time() - t1
			#利用网络计算WGCC大小
        with Bar('NO ADAPTIVE - MC:', max=M) as bar:
            t2 = time.time()
            for i in range(M):
                G_noad2 = G_noad.copy() #复制网络
                assign_weights(G_noad2,self.weight_distribution_type,type) #分配权重（uniform）
                wgcc_n = []
                for y in yy:
                    remove_thresh(G_noad2,y,0)
                    sw_noad = len(max(nx.weakly_connected_components(G_noad2), key=len))/N
                    wgcc_n.append(sw_noad)
                WGCC_noad.append(wgcc_n)
                bar.next()
            Dt2 = time.time() - t2

        print('Model Time, No Adaptive:',Dt1)
        print('MC Time, No Adaptive:',Dt2)        
        
        return WGCC_noad_T,WGCC_noad

    #FIG2下：输入：阈值范围，网络，权重-度联合分布，入度序列，入度分布，权重序列，权重分布，误差。输出：更新后的权重，经过有响应渗流多次重分配后的权重序列
    # adaptive_deg.degradation_adaptive_multi(yy,G,pwk,kk,pk0in,ww,pw,TOLLK,scenario)
    def degradation_adaptive_multi(self,yy,G,pwk,powk,kk,pk,ww,pw,TOLLK,scenario):
        dw = ww[1] - ww[0] #步长
        PW_ad_T = [] 
        WW_ad = []        

        with Bar('ADAPTIVE',max=len(yy)) as bar:
            # G2=G.copy()
            for y in yy:
                wkpred = pwk #权重-度联合分布
                kpred = pk #入度分布
                wpred = pw #权重分布


                pwk_ = PFKEX(wkpred,powk,y,kk,ww,dw,kpred,TOLLK,1,scenario) #联合概率分布f*k平均
                    
                kk = np.arange(0,len(pwk_)+1,1) #更新后的度序列

                pkex = [] #余度分布lk的分子部分
                for k in range(0,len(pwk_)):
                    pkex.append(integ(pwk_[k],dw))
                meank = np.sum(pkex) #平均度

                #更新后的损伤度分布
                pk = [0]*len(kk)
                for k in range(1,len(kk)): 
                    pk[k] = pkex[k-1]/k
                pk[0] = 1 - np.sum(pk)


                pwk = pwk_/meank #联合概率分布f
                pkex = pkex/meank #损伤后的余度分布

                #更新后的权重分布
                pw = 0
                for q in pwk:
                    pw += q

                PW_ad_T.append(pw) #预测的权重分布
                
                # G2=G.copy()

                if self.sce==1:
                    remove_thresh(G,y,1)#无响应渗流，没有重分配
                elif self.sce == 2:
                    remove_thresh2(G,y,1) #情景2
                elif self.sce == 3:
                    remove_thresh3(G,y,1) #情景3
                elif self.sce == 4:
                    remove_thresh4(G,y,1) #情景4
                else:
                    remove_thresh5(G,y,1) #情景5

                w_list = [G[e[0]][e[1]]['weight'] for e in G.edges()]
                WW_ad.append(w_list) #实验模拟的权重分布



                bar.next()

        return PW_ad_T,WW_ad

    #FIG3：多步响应的WGCC大小。输入：次数，阈值，网络，pwk，入度序列，入度分布，出度分布，权重序列，权重分布，偏差
    def degradation_adaptive_multi_WGCC(self,M,yy,G,pwk,kk,pk,pk0out,ww,pw,TOLLK):
        dw = ww[1] - ww[0]
        WGCC_T = []
        WGCC = []

			#使用入度-出度外积计算的WGCC大小
        with Bar('ADAPTIVE - MODEL:',max=len(yy)) as bar:
            t1 = time.time()
            for y in yy:
                wkpred = pwk
                kpred = pk
                wpred = pw

                pwk_ = PFKEX_APPROX(wkpred,y,kk,ww,dw,kpred,TOLLK,1,12) #更新后的联合概率密度函数fk*k平均
                    
                kk = np.arange(0,len(pwk_)+1,1) #更新后的入度序列

                pkex = [] #余度分布的分子：对联合分布f*k平均积分
                for k in range(0,len(pwk_)):
                    pkex.append(integ(pwk_[k],dw)) 
                meank = np.sum(pkex) #平均度

                pk = [0]*len(kk) #入度概率分布
                for k in range(1,len(kk)): 
                    pk[k] = pkex[k-1]/k #余度分布的分子除以度 等于入度概率分布
                pk[0] = 1 - np.sum(pk)

                pkout = PFK_OUT(kpred,wpred,y,ww,dw,kk) #出度概率分布
            
                WGCC_T.append(fixed_point_directed(np.outer(pkout,pk))) #利用入度-出度概率分布计算的WGCC

                pwk = pwk_/meank #联合概率分布f
                pkex = pkex/meank

 					#pw为更新后的网络权重分布，通过累加联合概率分布f得到
                pw = 0
                for q in pwk:
                    pw += q 

                bar.next()
            Dt1 = time.time() - t1
            #使用网络模拟得到有响应的WGCC大小
        with Bar('ADAPTIVE - MC:',max=M) as bar:
            t2 = time.time()
            for i in range(M):
                G2 = G.copy()
                assign_weights(G2,self.weight_distribution_type,type)
                wgcc = []
                for y in yy:
                    # remove_thresh(G2,y,1)
                    if self.sce==1:
                        remove_thresh(G2,y,1)#无响应渗流，没有重分配
                    elif self.sce == 2:
                        remove_thresh2(G2,y,2) #情景2
                    elif self.sce == 3:
                        remove_thresh3(G2,y,3) #情景3
                    elif self.sce == 4:
                        remove_thresh4(G2,y,4) #情景4
                    else:
                        remove_thresh5(G2,y,5) #情景5
                    
                    if not nx.is_empty(G2):	
                        sw = len(max(nx.weakly_connected_components(G2), key=len))/N
                        wgcc.append(sw)
                    else:
                        wgcc.append(0)
                WGCC.append(wgcc)
                bar.next()
            Dt2 = time.time() - t2       

        

        print('Model Time, Adaptive:',Dt1)
        print('MC Time, Adaptive:',Dt2)

        return WGCC_T,WGCC
#FIG4：得到50单步后的平均权重，平均入度值，协方差。一步的损伤平均权重，平均入度，协方差，每步阈值。
#输入：阈值，网络，权重序列，权重概率分布，入度序列，入度概率分布
    def degradation_adaptive_single(self,yy,G,ww,pw,kk,pk,scenario):

        MEANSW = []#平均权重
        MEANSK = []#平均入度值
        RHO = [] #协方差
        yy1 = [] #阈值
        dw = ww[1] - ww[0] #步长
        w_max = ww[-1] #最大权重

        with Bar('MC',max=len(yy)) as bar: #次数
            t1 = time.time()
            for y in yy:
                

                G2 = G.copy()#单步

                if scenario==1:
                    remove_thresh(G2,y,1)#无响应渗流，没有重分配
                elif scenario == 2:
                    remove_thresh2(G2,y,1) #情景2
                elif scenario == 3:
                    remove_thresh3(G2,y,1) #情景3
                elif scenario == 4:
                    remove_thresh4(G2,y,1) #情景4
                else:
                    remove_thresh5(G2,y,1) #情景5
                
                

                
                if not nx.is_empty(G2):
                    w = []
                    wk = []
                    kn = []
                    for e in G2.edges():
                        w.append(G2[e[0]][e[1]]['weight']) #边权
                        wk.append(G2.in_degree(e[1]) * G2[e[0]][e[1]]['weight']) #边的target节点的入度*边权
                        kn.append(G2.in_degree(e[1])) #边的target节点的入度
                    
                    w_mean = np.mean(w) #网络的平均权重
                    MEANSW.append(w_mean)
                    
                    wk_mean = np.mean(wk) #网络的所有边的平均度权
                    kn_mean = np.mean(kn) #网络的所有边的平均入度
                
                    kin_list = [G2.in_degree(n) for n in G2.nodes()] #网络的入度序列
                    k_mean = np.mean(kin_list) #网络的平均入度
                    
                    MEANSK.append(k_mean) #网络的平均入度
                    
                    RHO.append(wk_mean - w_mean*kn_mean) 
                    
                    yy1.append(y)
                bar.next()
            Dt1 = time.time() - t1

        pkex = PFK_EX(pk,kk) #余度分布

        t2 = time.time()
        
        k_0 = mean_k(pk) #平均度
        kn_0 = np.sum([k*pkex[k-1] for k in range(1,max_in_deg(G) + 1)]) #初始平均度

        beta1, beta2 = beta1_2(pw,ww,dw,w_max)
        F = primitive(pw,ww,dw) #累积分布F（y）
        GF = GEN(F,pkex) #生成函数G（F(y)）

        AVGW_T = beta1 + beta2*( (F - GF) / (1 - F) ) 
        AVGK_T = k_0*(1 - F) #损伤后的平均度
        RHO_T = beta2*(kn_0*GF - (F*(1-GF))/(1-F)) 
        
        Dt2 = time.time()-t2

        print('MC Time:',Dt1)
        print('Model Time:',Dt2)
        
        

        return MEANSW,MEANSK,RHO,AVGW_T,AVGK_T,RHO_T,yy1

#FIG5：单步，有无响应的最大lambda值和网络平均强度
    def degradation_single_maxeig(self,M,yy,G,ww,pw,kk,pk):

        L_max_resp = [] #有响应的lambda最大值（50个阈值*25次）
        L_max_noresp = []#无响应的lambda最大值
        dw = ww[1] - ww[0] #步长
        w_max = ww[-1] #最大权重

        with Bar('MC, Single',max=len(yy)) as bar: #阈值50
            t1 = time.time()
            for y in yy: #每个阈值
                lambdas_noresp = []
                lambdas_resp = []

                for i in range(M): #25次
                    G1 = G.copy()
                    assign_weights(G1,self.weight_distribution_type)
                
                    GG1 = G1.copy()
                    remove_thresh(GG1,y,1) 
                    lmax = np.real(max(LA.eigvals(nx.to_numpy_matrix(GG1)))) #得到有响应过滤后的最大特征值
                    lambdas_resp.append(lmax)

                    GG2 = G1.copy()
                    remove_thresh(GG2,y,0)
                    lmax = np.real(max(LA.eigvals(nx.to_numpy_matrix(GG2)))) #得到无响应过滤后的最大特征值
                    lambdas_noresp.append(lmax)

                L_max_noresp.append(lambdas_noresp)
                L_max_resp.append(lambdas_resp)
                bar.next()

            L_max_resp = np.array(L_max_resp)
            L_max_noresp = np.array(L_max_noresp)
            Dt1 = time.time() - t1

        t2 = time.time()

        pkex = PFK_EX(pk,kk)#余度分布
        F = primitive(pw,ww,dw) #累积概率分布F（y）
        beta1, beta2 = beta1_2(pw,ww,dw,w_max) 
        GF = GEN(F,pkex) #生成函数
        k_0 = mean_k(pk) #平均度


        AVGW_T = beta1 + beta2*( (F - GF) / (1 - F) ) #解析结果的平均权重 
        AVGK_T = k_0*(1 - F) #解析结果的损伤后平均度值
        AVG_S_resp = AVGW_T*AVGK_T #有响应的网络的平均强度，
        AVG_S_noresp = beta1*AVGK_T #无响应的网络的平均强度

        Dt2 = time.time() - t2

        print('MC Time, Single:',Dt1)
        print('Model Time, Single:',Dt2)

        return L_max_resp,L_max_noresp,AVG_S_resp,AVG_S_noresp

#FIG5：多步响应
    def degradation_multi_maxeig(self,M,yy,G,pwk,ww,pw,kk,pk,TOLLK):

        dw = ww[1] - ww[0] #步长
        w_max = ww[-1] #最大权重
        MEANSS_T = [] 
        meank = mean_k(pk) #网络平均度

        t1 = time.time()
        
        F = primitive(pw,ww,dw) #累计分布F（y）
        beta1, beta2 = beta1_2(pw,ww,dw,w_max)
        AVGK_T = meank*(1 - F) #损伤后平均度
        AVG_S_noresp = beta1*AVGK_T #无响应的网络平均强度
        
        with Bar('ADAPTIVE - MODEL, Multi',max=len(yy)) as bar: 
            for y in yy:
                wkpred = pwk
                kpred = pk
                wpred = pw

                pwk_ = PFKEX(wkpred,y,kk,ww,dw,kpred,TOLLK,1) ##联合概率分布f*k平均
                        
                kk = np.arange(0,len(pwk_)+1,1) #入度序列
                #pkex为余度分布的分子
                pkex = []
                for k in range(0,len(pwk_)):
                    pkex.append(integ(pwk_[k],dw))
                meank = np.sum(pkex) #平均度
                #更新后度概率分布
                pk = [0]*len(kk)
                for k in range(1,len(kk)): 
                    pk[k] = pkex[k-1]/k
                pk[0] = 1 - np.sum(pk)
                    
                pwk = pwk_/meank #联合概率密度分布f
                pkex = pkex/meank  #余度分布
                #pw为损伤后权重分布
                pw = 0
                for q in pwk:
                    pw += q
                    
                w_m = integ(pw*ww,dw) #平均权重

                MEANSS_T.append(w_m*meank) #网络的平均强度
                bar.next()

        Dt1 = time.time() - t1
        #使用网络模拟得到有/无响应两种情况下最大特征值
        L_resp = []
        L_noresp = []

        t2 = time.time()
        with Bar('MC, Multi',max=len(yy)) as bar:
            for i in range(0,M):
                G_resp = G.copy()
                G_noresp = G.copy()

                assign_weights(G_resp,self.weight_distribution_type)
                assign_weights(G_noresp,self.weight_distribution_type)

                l_resp = []
                l_noresp = []

                for f in yy:     
                    remove_thresh(G_resp,f,1)
                    lmax = np.real(max(LA.eigvals(nx.to_numpy_matrix(G_resp)))) 
                    l_resp.append(lmax)

                    remove_thresh(G_noresp,f,0)
                    lmax = np.real(max(LA.eigvals(nx.to_numpy_matrix(G_noresp)))) 
                    l_noresp.append(lmax)

                L_resp.append(l_resp)
                L_noresp.append(l_noresp)
                bar.next()
    
        L_resp = np.array(L_resp)
        L_noresp = np.array(L_noresp)
        Dt2 = time.time() - t2 

        print('Model Time, Multi:',Dt1) 
        print('MC Time, Multi:',Dt2)
          

        return L_noresp,L_resp,AVG_S_noresp,MEANSS_T

    def create_network(self,G,type,n_node,n_edge):
        G_temp=G
        #ER随机网络。含有n_node个节点、n_edge边的有向ER随机图
        if(type=="random"):
            G=nx.random_graphs.gnm_random_graph(n_node, n_edge,directed=True)
        #BA无标度网络。含有n个节点的有向BA无标度网络
        if(type=="ba"):
            G=nx.scale_free_graph(n_node)
            G=nx.DiGraph(G)
        #小世界网络。含有n个节点、每个节点有k（平均度值）个邻居、以概率p随机化重连边的WS小世界网络。节点数和边数与实际网络相同
        if(type=="sw"):
            kk=self.initialize_k_ax(G)
            k_list = list([G.degree(n) for n in G.nodes()]) #度序列
            pk = [k_list.count(i)/len(k_list) for i in kk] #每个度值的概率，即度概率密度分布
            out = []
            for k in range(1,len(kk)):
                out.append(k*pk[k])
            m=round(np.sum(out))
            G=nx.watts_strogatz_graph(n_node, m, p=0.3)
            G=G.to_directed()
            #随机删除多余边
            diff=len(G.edges)-n_edge
            if diff>0:
                G2=G
                dele=sample(G.edges,diff)
                G2.remove_edges_from(dele)
                #判断网络中有无孤立点（是否有节点的度值为0）
                while (np.array(list([G2.degree(n) for n in G2.nodes()]))==0).any():
                    G2=G
                    dele=sample(G.edges,diff)
                    G2.remove_edges_from(dele)
                G = G2
        if(type=="real"):
            G=G_temp
        #分配边的权重
        assign_weights(G,self.weight_distribution_type,type)
        
        return G
            


#line 200：应该为pkout = PFK_OUT(pk0out,wpred,y,ww,dw,kk)
       
    #指标3
    # adaptive_deg.degradation_adaptive_multi(yy,G,pwk,kk,pk0in,ww,pw,TOLLK,scenario)
    def caculate_ud(self,yy,G,TOLLK,scenario):

        UD=[]

        with Bar('ADAPTIVE',max=len(yy)) as bar:
            G2=G.copy()
            for y in yy:
                #指标3
                A=np.array(nx.adjacency_matrix(G2).todense()) #前一步网络邻接01矩阵
                rows,cols=A.shape#读取A的行列数
                for i in range(rows):
                    for j in range(cols):
                        if A[i][j]:
                            A[i][j]=1#用按照行和列定义X中的每个元素为a
                W=np.array(nx.adjacency_matrix(G2).todense()) #前一步网络邻接权重矩阵
                X=np.zeros((int(rows),int(cols)))#未响应矩阵
                for i in range(cols):
                    if sum(A[:,i])==1:
                        if sum(W[:,i])<=y:
                            for j in range(rows):
                                if A[j][i]==1:
                                    X[j][i]=1
                R=A-X
                W_ini=np.array(nx.adjacency_matrix(G).todense()) #原始网络邻接权重矩阵
                ud=np.dot(R, W_ini.T).trace()/np.sum(W_ini)
                UD.append(ud)

                if scenario == 0:
                    remove_thresh2(G2,y,0)#无响应渗流，没有重分配
                elif scenario == 1:
                    remove_thresh(G2,y,1) #情景1
                elif scenario == 2:
                    remove_thresh2(G2,y,1) #情景2
                elif scenario == 3:
                    remove_thresh3(G2,y,1) #情景3
                elif scenario == 4:
                    remove_thresh4(G2,y,1) #情景4
                else:
                    remove_thresh5(G2,y,1) #情景5

                bar.next()

        return UD

    #计算连通率。输入：次数，阈值，网络，pwk，入度序列，入度分布，出度分布，权重序列，权重分布，偏差
    def caculate_connectivity(self,M,yy,G,scenario):
        N_node1=[]#原始网络节点数
        N_node2=[]#重配后网络节点数
        N_edge1=[]#原始网络边数
        N_edge2=[]#重配后网络边数
        con=[]
        N = nx.number_of_nodes(G)

        with Bar('ADAPTIVE - MC:',max=M) as bar:
            t2 = time.time()
            for i in range(M):
                G2 = G.copy()
                # assign_weights(G2,self.weight_distribution_type,type)
                for y in yy:
                    # remove_thresh(G2,y,1)
                    if scenario == 0:
                        remove_thresh2(G2,y,0)#无响应渗流，没有重分配
                    elif scenario == 1:
                        remove_thresh(G2,y,1) #情景1
                    elif scenario == 2:
                        remove_thresh2(G2,y,1) #情景2
                    elif scenario == 3:
                        remove_thresh3(G2,y,1) #情景3
                    elif scenario == 4:
                        remove_thresh4(G2,y,1) #情景4
                    else:
                        remove_thresh5(G2,y,1) #情景5
                #指标2
                
                    N_node1=G.number_of_nodes()
                    N_edge1=G.number_of_edges()
                    tempt1=1/N_node1/(N_node1-1)*N_edge1
                    N_node2=G2.number_of_nodes()
                    N_edge2=G2.number_of_edges()
                    tempt2=N_edge2/N_node2/(N_node2-1)
                    # tempt1=1
                    con.append(tempt2/tempt1)

                bar.next()


        return con

#平均度值和权重，多步（根据figure4改写）
#输入：阈值，网络，权重序列，权重概率分布，入度序列，入度概率分布
    def caculate_meanskw(self,yy,G,scenario):

        MEANSW = []#平均权重
        MEANSK = []#平均入度值
        yy1 = [] #阈值


        with Bar('MC',max=len(yy)) as bar: #次数
            t1 = time.time()
            G2 = G.copy()#多步
            for y in yy:
                
                

                if scenario == 0:
                    remove_thresh2(G2,y,0)#无响应渗流，没有重分配
                elif scenario == 1:
                    remove_thresh(G2,y,1) #情景1
                elif scenario == 2:
                    remove_thresh2(G2,y,1) #情景2
                elif scenario == 3:
                    remove_thresh3(G2,y,1) #情景3
                elif scenario == 4:
                    remove_thresh4(G2,y,1) #情景4
                else:
                    remove_thresh5(G2,y,1) #情景5
                

                if not nx.is_empty(G2):
                    w = []
                    wk = []
                    kn = []
                    for e in G2.edges():
                        w.append(G2[e[0]][e[1]]['weight']) #边权
                        wk.append(G2.in_degree(e[1]) * G2[e[0]][e[1]]['weight']) #边的target节点的入度*边权
                        kn.append(G2.in_degree(e[1])) #边的target节点的入度
                    
                    w_mean = np.mean(w) #网络的平均权重
                    MEANSW.append(w_mean)
                    
                    wk_mean = np.mean(wk) #网络的所有边的平均度权
                    kn_mean = np.mean(kn) #网络的所有边的平均入度
                
                    kin_list = [G2.in_degree(n) for n in G2.nodes()] #网络的入度序列
                    k_mean = np.mean(kin_list) #网络的平均入度
                    
                    MEANSK.append(k_mean) #网络的平均入度
                if nx.is_empty(G2):
                    MEANSK.append(0)
                    MEANSW.append(0)
                
                bar.next()


        return MEANSW,MEANSK
