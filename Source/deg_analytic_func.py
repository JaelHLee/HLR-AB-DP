import numpy as np
from scipy.special import binom
from numpy.fft import fft,ifft
import math
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
from scipy.stats import powerlaw 
import powerlaw

def kNN(data_set, x_list,k):
    # assign number
    N = data_set.shape[0]
    # initial knn list
    knn_list = []
    for x in x_list:
        u = np.linalg.norm(data_set.reshape(-1,1)-x, axis=1)
        r = np.sort(u)[k-1]
        v = 2*r
        p = k/(N*v)
        knn_list.append(p)
    return knn_list,np.sum(np.log(knn_list))

def pl_knn(arr):
    # assign independent variable
    arr=np.array(arr)
    x_list = np.linspace(0, 1, arr.shape[0])
    # assign k list
    k=18
    knn, likelihood = kNN(arr, x_list, k)
    return knn
        

#积分，和*步长
def integ(arr,df):
	return np.sum(arr)*df
#得到初始累计权重分布F（y）
def primitive(arr,ww,dw):
	out = []
	for t in ww:
		out.append(integ(arr*np.heaviside(t-ww,1),dw))
	return np.array(out)
#标准化
def norm(arr,dw):
	if(integ(arr,dw) != 0):
		return arr/integ(arr,dw)
	if(integ(arr,dw) == 0):
		return arr
#链接两个数组
def scalex(arr,fac):
	split = arr[::int(fac)] #arr的每第int（fac）个序列，每隔fac间隔，取一个数，一直到结尾end
	add = np.zeros(len(arr)-len(split))#返回零数组,数值为差值个0
	return np.concatenate((split,add))#链接两个数组
#生成正态概率密度分布(省略了系数，但结果不影响)
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) 
#生成幂律分布
#计算卷积：计算a和b的快速傅里叶变换，然后得到二者的逆傅里叶变换的实数部分
def conv_norm(a,b,n,df):	
	fa = fft(a) #快速傅里叶变换
	fb = (fft(b)*df)**n
	return np.real(ifft(fa*fb)) #计算逆傅里叶变换，并返回实数部分
#平均度
def mean_k(pk):
	s = 0
	for k in range(0,len(pk)):
		s += pk[k]*k
	return s
#图2有响应：输入：权重-度联合分布，阈值，入度序列，权重序列，步长，入度分布，偏差，是否多步响应。输出：过滤响应后的余度分布的分子
def PFKEX(wkpred,powk,f,kk,ww,dw,kpred,TOLLK,multi,scenario):
	out = []
	check = 0
	for k in range(1,len(kpred)):#带入每个入度值计算一遍，结果为out
		s = 0
		for j in range(k,len(kk)):
			a_k = norm(wkpred[j-1],dw) #度权联合分布的标准化
			a = integ(a_k*np.heaviside(ww-f,1),dw) #Fn(y):筛选后的权重累积分布函数.比较与阈值的大小，若小于a=0,否则a=a_k*步长，求积分
			wpred_up = wkpred[j-1]*np.heaviside(ww-f,1) #根据阈值筛选联合分布 
			if(scenario == 1):
				wpredk_down = scalex(wkpred[j-1],k)*np.heaviside(f - (k*ww),1) 
			elif(scenario == 2):
				wpredk_down = k*ww*np.heaviside(ww - f,1)/np.sum(ww*np.heaviside(ww-f,1))*wkpred[j-1]*np.heaviside(f-ww,1) 
			elif(scenario == 3):
				wpredk_down = powk[j-1]*np.heaviside(f-ww,1)/np.sum(powk[j-1]*np.heaviside(ww-f,1)+1e-26)*[wkpred[j-1]*np.heaviside(f - ww,1)]
				wpredk_down = wpredk_down.T
				wpredk_down = wpredk_down[:,0]
            #找入边中权重最大的边
			elif(scenario == 4):
				temp=np.argmax(ww*np.heaviside(f - ww,1))
				w_=ww
				w_=[0 for i in w_]
				w_[temp]=1
				wpredk_down = w_*wkpred[j-1]*np.heaviside(f - ww,1) 
            #场景5：最大熵
			elif(scenario == 5):
				m_temp=list(np.heaviside(ww - f,1)).count(1) #留下的边数
				w_sum = np.sum(wkpred[j-1])
				w_temp=wkpred[j-1]*np.heaviside(wkpred[j-1]*np.heaviside(ww - f,1)*m_temp-w_sum,1)
				m=list(np.heaviside(w_temp - w_sum,1)).count(1)
				if m:
					wpredk_down = (w_sum/m-w_temp)*[wkpred[j-1]*np.heaviside(f - ww,1)]
					wpredk_down = wpredk_down.T
					wpredk_down = wpredk_down[:,0]
				else:
					wpredk_down=wkpred[j-1]*np.heaviside(f - ww,1)
					wpredk_down = wpredk_down.T

			conv_ = conv_norm(wpred_up,wpredk_down,(j-k),dw) #卷积运算

			if(np.sum(conv_ > 1e6) != 0):
				conv_[conv_ > 1e6] = 0
			if(np.sum(conv_ < 0) != 0):
				conv_[conv_ < 0] = 0

			conv = norm(conv_,dw)

			s += binom(j,k)*(a**k)*((1-a)**(j-k))*kpred[j]*conv	 

		out.append(s*k)
	if(multi == 1): 
		if(len(out) > 1):
			while True:
				a = np.array(out[-1]) #取最后一项
				p = integ(a,dw)
				if(p > TOLLK):
					break
				else:
					del out[-1] #删除最后一项
			
	return np.array(out)
#图3有响应：得到更新后的联合概率密度函数fk*k平均. 设置了最大度值dkmax=12
def PFKEX_APPROX(wkpred,f,kk,ww,dw,kpred,TOLLK,multi,dkmax):
	out = []
	for k in range(1,len(kpred)):
		s = 0
		km =  len(kk)
		if(k < len(kk) - dkmax):
			km = k + dkmax + 1
		for j in range(k,km):
			a_k = norm(wkpred[j-1],dw) #权重标准化
			a = integ(a_k*np.heaviside(ww-f,1),dw) #权重概率密度函数的累积分布Fk（y）

			wpred_up = wkpred[j-1]*np.heaviside(ww-f,1) #大于阈值的权重分布
			wpredk_down = scalex(wkpred[j-1],k)*np.heaviside(f - (k*ww),1) #小于阈值的权重分布

			conv_ = conv_norm(wpred_up,wpredk_down,(j-k),dw)#卷积运算

			if(np.sum(conv_ > 1e6) != 0):
				conv_[conv_ > 1e6] = 0
			if(np.sum(conv_ < 0) != 0):
				conv_[conv_ < 0] = 0

			conv = norm(conv_,dw)

			s += binom(j,k)*(a**k)*((1-a)**(j-k))*kpred[j]*norm(conv,dw)		

		out.append(s*k)
	if(multi == 1):
		if(len(out) > 1):
			while True:
				a = np.array(out[-1])
				p = integ(a,dw)
				if(p > TOLLK):
					break
				else:
					del out[-1]
			
	return np.array(out)
#图3无响应：输入：入度分布，阈值，权重分布，步长，误差。输出：损伤后的度分布Pk（y）
def PFK_NOAD(pk0,f,pw0,ww,dw,TOLLK):
	out = []
	a = integ(pw0*np.heaviside(ww-f,1),dw) #判断与阈值的大小，更新权重（大于阈值不变，小于阈值为0）。然后积分得F累计分布
	for k in range(0,len(pk0)):
		s = 0
		for j in range(k,len(pk0)):			
			s += binom(j,k)*(a**k)*((1-a)**(j-k))*pk0[j] #二项式累积概率函数
			
		out.append(s)
	if(len(out) > 1):
		while True:
			a = np.array(out[-1]) #取最后一个值
			if(a > TOLLK):
				break
			else:
				del out[-1]
	
	return np.array(out)
#求余度分布：excess degree distribution （lk）
#图4：输入：网络入度概率分布，入度序列。输出：余度分布lk
def PFK_EX(pfk,kk):
	out = []
	for k in range(1,len(kk)): #k=1,2,,,(最大入度值-1)
		out.append(k*pfk[k])
	out = np.array(out)
	return out/np.sum(out)
#计算出度概率分布。输入：入度分布，权重分布，阈值，权重序列，步长，入度序列
def PFK_OUT(kpred,wpred,f,ww,dw,kk):
	out = []
	a = integ(wpred*np.heaviside(ww-f,1),dw) #权重的累积分布Fk（y）
	for k in kk:
		sum = 0
		for j in range(k,len(kk)):
			sum += binom(j,k)*(a**k)*((1-a)**(j-k))*(kpred[j])

		out.append(sum) #损失后的度分布Pk（y）
	
	return(np.array(out))
	
def beta1_2(pw,ww,dw,wmax):
	eps = 1e-5
	b2 = []
	b1 = []
	w0 = integ(ww*pw,dw) 
	for t in ww:
		x1 = integ(ww*pw*np.heaviside(t-ww,1),dw) 
		x2 = integ(pw*np.heaviside(t-ww,1),dw) #F（y）
		b2.append(x1/x2) #beta2
		if(np.abs((w0 - x1)/(1 - x2)) > eps):
			b1.append((w0 - x1)/(1 - x2)) #beta1 
		else:
			b1.append(t)		
	b1 = np.array(b1)
	b2 = np.array(b2)
			
			
	return b1, b2
#生成函数G：余度分布*累积分布的累加
def GEN(F,pkex):
	out = 0
	for k in range(1,len(pkex)+1):
		out += pkex[k-1] * F**k
	return out


powerlaw		
