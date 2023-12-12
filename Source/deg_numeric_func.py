import networkx as nx
import numpy as np
import random
#计算网络G中节点为n的入强度
def in_strength(G,n):
	eout = G.in_edges(n)#获取节点为n的入边
	s = 0
	for e in eout:
		s += G[e[0]][e[1]]['weight']
	return s
#计算网络G中节点为n的出强度
def out_strength(G,n):
	eout = G.out_edges(n)#应为out_edges(n)
	s = 0
	for e in eout:
		s += G[e[0]][e[1]]['weight']
	return s

#返回网络G节点中最大的出度
def max_out_deg(G):
	deg = [G.out_degree(n) for n in G.nodes]
	return max(deg)
#返回网络G节点中最大的入度
def max_in_deg(G):
	deg = [G.in_degree(n) for n in G.nodes] #入度序列
	return max(deg)
#权重序列
def weight_seq2(G):
    weight1 = []
    for e in G.edges():
        weight1.append(G[e[0]][e[1]]['weight'])
    return weight1


#生成权重函数：若权重分布为uniform，则随机生成0-1的随机数作为权重。若权重分布为gauss，则生成正态分布的随机浮点数，参数为（平均值，标准差）
def assign_weights(G,dist,type):
	if(dist == 'uniform'):
		for e in G.edges():
			G[e[0]][e[1]]['weight'] = random.random()
	if(dist == 'gauss'):
		for e in G.edges():
			G[e[0]][e[1]]['weight'] = random.normalvariate(0.5,0.1)

	if(dist == 'no'):
		weight=weight_seq2(G)
		maxvalue=max(weight)
		minvalue=min(weight)
		for e in G.edges():
			G[e[0]][e[1]]['weight']=(G[e[0]][e[1]]['weight']-minvalue)/(maxvalue-minvalue)


        
        

#原始情景：渗流和重分配。输入：网络，阈值，是否响应（0、1）
def remove_thresh(G,f,a):
	for n in G.nodes():
		nn = [m for m in G.predecessors(n)]#获得节点n的前置节点
		to_remove = []
		Delta = 0
		for m in nn: #渗流过程：统计边权小于阈值的边，并计算总损失权重
			if(G[m][n]['weight'] < f):
				to_remove.append([m,n])
				Delta += G[m][n]['weight']
		for e in to_remove: #从网络G中删除边权小于阈值的边
			G.remove_edge(e[0],e[1])
# 		G.remove_nodes_from(list(nx.isolates(G)))
		nn_new = [m for m in G.predecessors(n)] #获得节点n的新的前置节点
		for m in nn_new:
			G[m][n]['weight'] += (Delta*a)/len(nn_new) #重分配。若a为0，则无响应；若a为1，则有响应
#情景2：根据边权按比例重新负载给存在的边；
def remove_thresh2(G,f,a):
    for n in G.nodes():
        nn=[m for m in G.predecessors(n)]
        to_remove=[]
        Delta=0
        re_we=0 #更新后节点的入边权重和
        for m in nn:
            if(G[m][n]['weight'] < f):
                to_remove.append([m,n])
                Delta += G[m][n]['weight']
            else:
                re_we += G[m][n]['weight']
        for e in to_remove:
            G.remove_edge(e[0],e[1])
        # G.remove_nodes_from(list(nx.isolates(G)))
        nn_new=[m for m in G.predecessors(n)]
        for m in nn_new:
            G[m][n]['weight']+= (Delta*a)*(G[m][n]['weight']/re_we)
#情景3：根据出口国出口总量重新分配
def remove_thresh3(G,f,a):
    for n in G.nodes():
        nn=[m for m in G.predecessors(n)]
        to_remove=[]
        Delta=0
        for m in nn:
            if(G[m][n]['weight'] < f):
                to_remove.append([m,n])
                Delta += G[m][n]['weight']
        for e in to_remove:
            G.remove_edge(e[0],e[1])
        # G.remove_nodes_from(list(nx.isolates(G)))
        nn_new=[m for m in G.predecessors(n)]
        out_sum=0 #更新后网络中节点n所有前置节点的出强度总和
        for m in nn_new:
           out_sum =out_sum+ out_strength(G, m)
        temp=1e-26
        for m in nn_new:
             G[m][n]['weight'] += (Delta*a)*(out_strength(G, m)/(out_sum+temp))
#情景4：选最大的一条边负载
def remove_thresh4(G,f,a):
    for n in G.nodes():
        nn=[m for m in G.predecessors(n)]
        to_remove=[]
        Delta=0
        for m in nn:
            if(G[m][n]['weight'] < f):
                to_remove.append([m,n])
                Delta += G[m][n]['weight']
        for e in to_remove:
            G.remove_edge(e[0],e[1])

        nn_new = [m for m in G.predecessors(n)]
        if nn_new:
            out_seq = [G[m][n]['weight'] for m in nn_new] #更新后网络中节点n的所有入边权重序列
            if out_seq:
                temp=out_seq.index(max(out_seq))
                G[nn_new[temp]][n]['weight'] +=  Delta*a
    # G.remove_nodes_from(list(nx.isolates(G)))

#情景5：最大熵重分配
def remove_thresh5(G,f,a):
    for n in G.nodes():
        nn=[m for m in G.predecessors(n)]
        Delta2=in_strength(G,n) #节点n的初始入强度
        #更新网络
        to_remove=[]
        for m in nn:
            if(G[m][n]['weight'] < f):
                to_remove.append([m,n]) #删除边的节点
        for e in to_remove:
            G.remove_edge(e[0],e[1]) #更新网络，删除边
        if a:
        # G.remove_nodes_from(list(nx.isolates(G)))
        #更新后节点n的前置节点
            nn_new = [m for m in G.predecessors(n)]
        #更新后节点n的前置节点个数
            n_node=len(nn_new)
            if n_node:
                thre=Delta2/n_node
                nn_thre=[]
                weight_thre=0
                for m in nn_new:
                    if(G[m][n]['weight']<=thre):
                        nn_thre.append([m,n]) #满足分配条件的前置节点
                    else:
                        weight_thre += G[m][n]['weight'] #不参与分配的边权和
                        G[m][n]['weight']
                       
                        
                    
        #仅更新满足分配条件的边
                if nn_thre:
                    for m in nn_thre:
                        G[m[0]][m[1]]['weight'] = (Delta2-weight_thre)*a/len(nn_thre)
                        G[m[0]][m[1]]['weight']
                       
                        

            
       
    
