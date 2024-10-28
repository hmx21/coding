"""
1-site variational MPS method
Spin-1/2 antiferromagnetic Heisenberg model on an open chain

Copyright:
Shuo Yang, shuoyang@tsinghua.edu.cn
Nov 1, 2021, Tsinghua, Beijing, China
"""

import numpy as np
import numpy.linalg as LA
import scipy.sparse.linalg as LAs
import Sub180221 as Sub
import math,copy
#----------------------------------------------------------------
"""
Mi =
S0	Sp	Sm	Sz	0
0	0	0	0	Sm/2
0	0	0	0	Sp/2
0	0	0	0	Sz
0	0	0	0	S0
"""
def GetMpo_Heisenberg_Obc(Dp):
	S0,Sp,Sm,Sz,Sx,Sy = Sub.SpinOper(Dp)
	
	Dmpo = 5
	Mpo = np.zeros((Dmpo,Dp,Dmpo,Dp))
	
	Mpo[0,:,0,:] = S0
	Mpo[0,:,1,:] = Sp
	Mpo[0,:,2,:] = Sm
	Mpo[0,:,3,:] = Sz
	Mpo[1,:,4,:] = Sm/2.0
	Mpo[2,:,4,:] = Sp/2.0
	Mpo[3,:,4,:] = Sz
	Mpo[4,:,4,:] = S0
	
	return Mpo

def InitMps(Ns,Dp,Ds):
	T = [None]*Ns
	for i in range(Ns):
		Dl = min(Dp**i,Dp**(Ns-i),Ds)
		Dr = min(Dp**(i+1),Dp**(Ns-1-i),Ds)
		T[i] = np.random.rand(Dl,Dp,Dr)
	
	U = np.eye(np.shape(T[-1])[-1])
	for i in range(Ns-1,0,-1):
		U,T[i] = Sub.Mps_LQP(T[i],U)
	
	return T

def InitH(Mpo,T):
	Ns = len(T)
	Dmpo = np.shape(Mpo)[0]
	
	HL = [None]*Ns
	HR = [None]*Ns
	
	HL[0] = np.zeros((1,Dmpo,1))
	HL[0][0,0,0] = 1.0
	HR[-1] = np.zeros((1,Dmpo,1))
	HR[-1][0,-1,0] = 1.0
	
	for i in range(Ns-1,0,-1):
		HR[i-1] = Sub.NCon([HR[i],T[i],Mpo,np.conj(T[i])],[[1,3,5],[-1,2,1],[-2,2,3,4],[-3,4,5]])
	
	return HL,HR

def OptTSite(Mpo,HL,HR,T,Method=0):
	DT = np.shape(T)
	Dl = np.prod(DT)
	
	if Method == 0:
		A = Sub.NCon([HL,Mpo,HR],[[-1,1,-4],[1,-5,2,-2],[-6,2,-3]])
		A = Sub.Group(A,[[0,1,2],[3,4,5]])
		Eig,V = LAs.eigsh(A,k=1,which='SA')
		T = np.reshape(V,DT)
	
	if Method == 1:
		def UpdateV(V):
			V = np.reshape(V,DT)
			V = Sub.NCon([HL,V,Mpo,HR],[[-1,3,1],[1,2,4],[3,2,5,-2],[4,5,-3]])
			V = np.reshape(V,[Dl])
			return V
		
		V0 = np.reshape(T,[Dl])
		
		MV = LAs.LinearOperator((Dl,Dl),matvec=UpdateV)
		Eig,V = LAs.eigsh(MV,k=1,which='SA',v0=V0)
		# print(Eig)
		T = np.reshape(V,DT)
		Eig = np.real(Eig)
	
	return T,Eig

def OptT(Mpo,HL,HR,T):
	Ns = len(T)
	Eng0 = np.zeros(Ns)
	Eng1 = np.zeros(Ns)
	
	for r in range(100):
		print(r)
	
		for i in range(Ns-1):
			T[i],Eng1[i] = OptTSite(Mpo,HL[i],HR[i],T[i],Method=1)
			# print(i,Eng1[i])
			T[i],U = Sub.Mps_QR0P(T[i])
			HL[i+1] = Sub.NCon([HL[i],np.conj(T[i]),Mpo,T[i]],[[1,3,5],[1,2,-1],[3,4,-2,2],[5,4,-3]])
			T[i+1] = np.tensordot(U,T[i+1],(1,0))
		
		for i in range(Ns-1,0,-1):
			T[i],Eng1[i] = OptTSite(Mpo,HL[i],HR[i],T[i],Method=1)
			# print(i,Eng1[i])
			U,T[i] = Sub.Mps_LQ0P(T[i])
			HR[i-1] = Sub.NCon([HR[i],T[i],Mpo,np.conj(T[i])],[[1,3,5],[-1,2,1],[-2,2,3,4],[-3,4,5]])
			T[i-1] = np.tensordot(T[i-1],U,(2,0))
		
		print(Eng1)
		if abs(Eng1[1]-Eng0[1]) < 1.0e-7:
			break
		Eng0 = copy.copy(Eng1)
	
	print(Eng1/float(Ns))
	
	return T
#----------------------------------------------------------------
if __name__ == "__main__":
	Ns = 10
	Dp = 2
	Ds = 4
	
	Mpo = GetMpo_Heisenberg_Obc(Dp)
	T = InitMps(Ns,Dp,Ds)
	HL,HR = InitH(Mpo,T)
	T = OptT(Mpo,HL,HR,T)
	