

# Solving Elastsicity 3D problem Deterministic using Wire Basket Grid Method


################################################################
## Sudhi Sharma P V : Modified July-2019
##
################################################################

#from petsc4py import PETSc
import numpy as np
import numpy.matlib
import scipy.io as sio
from scipy import linalg
import matplotlib.pylab as plt
import scipy.sparse as sparse
import time
import os as os

print("==========================================================")
print("Solving 3D deterministic elasticity - Wire Basket Grid...")


def DoPCGM_TwoLevel_WB(nDim,nParts,g_global,nBG,nCG,nRG,tol,tolc):


	u_gamma_n = np.zeros([nDim*nBG])

	q_global = np.zeros([nDim*nBG])

	r_gamma_n = g_global

	z_o = GetPreconditionedResidual_TwolevelWB(nParts,r_gamma_n,nDim,nCG,tolc)

	P_n = z_o

	rho_n = r_gamma_n.T.dot(z_o)

	print(rho_n)

	n = len(g_global)

	# print(n)

	z_n = z_o

	nres = np.zeros([n,1],dtype='float32')

	print("nres is",nres)

	# print(rho_n)


	for i in range(n):


		if i==0:
			print("PCGM Main Iteration Starts")

		# print("Iteration id ",i+1)

		Q = GetMatrixVectorProduct_WB(q_global,nParts,P_n)

		rho_temp = P_n.T.dot(Q)

		alpha = rho_n/rho_temp

		print("alpha is",alpha)

		# u_gamma_o = u_gamma_n

		# r_gamma_o = r_gamma_n

########### Exit Criteria Relative Error #############################

		d1 = float(u_gamma_n.dot(u_gamma_n))

		print("d1 is",d1)

		error = (alpha**2 * P_n.dot(P_n))/d1

		print("Iteration number "+str(i+1)+" with Relative Error of :",error)

		if(error < tol):
			break

#######################################################################

		u_gamma_n = u_gamma_n + alpha * P_n

		print("u update is",u_gamma_n)

		r_gamma_n = r_gamma_n - alpha * Q

		z_n = GetPreconditionedResidual_TwolevelWB(nParts,r_gamma_n,nDim,nCG,tolc)

		rho_o = rho_n

		rho_n = r_gamma_n.T.dot(z_n)

		beta = rho_n/rho_o

		print("Beta is",beta)

		# P_o = P_n

		P_n = z_n + beta * P_n

		print("New P_n",P_n)


########## Exit criteria Norm of Residual ########################
		# nres[i] = linalg.norm(r_gamma_n)/linalg.norm(g_global)
		# print("norm r gamma",linalg.norm(r_gamma_n))
		# print("norm Initial",linalg.norm(g_global))
		# print("Ratio of Norm of Residual is",nres[i])
		# if(nres[i] < tol):
		# 	break
###############################################################



	U_gamma = u_gamma_n
	print("U gamma is ",U_gamma)
	return U_gamma


##############################################################################

def GetMatrixVectorProduct_WB(q_global,nParts,P_n):

	q_global = np.zeros([nDim*nBG])

	for ip in range(nParts):

		npart = ip+1

		Rpath = "../../../data/meshData/R3D_"+str(npart)+".dat"
		R_3D = np.genfromtxt(Rpath)
		# print(np.shape(R_3D))

		BDpath = "../../../data/meshData/D3D_"+str(npart)+".dat"
		D_3D = np.genfromtxt(BDpath)
		# print(np.shape(D_3D))

		Rrpath = "../../../data/meshData/Rr3D_"+str(npart)+".dat"
		Rr_3D = np.genfromtxt(Rrpath)
		# print(np.shape(Rr_3D))

		Rcpath = "../../../data/meshData/Rc3D_"+str(npart)+".dat"
		Rc_3D = np.genfromtxt(Rcpath)
		# print(np.shape(Rc_3D))

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		# print(np.shape(Bc_3D))


		one = 1
		subpath = "../data/Amats/subdom"+str(npart)
		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADggpath = subpath+"/ADgg"+str(one)+".dat"
		ADgg = np.genfromtxt(ADggpath)

		ADigpath = subpath+"/ADig"+str(one)+".dat"
		ADig = np.genfromtxt(ADigpath)

		ADgi = ADig.T

		pb = R_3D.dot(P_n)
		u1 = ADig.dot(pb)
		u2 = linalg.solve(ADii,u1)
		u3 = ADgi.dot(u2)
		u4 = ADgg.dot(pb)
		q_local = u4 - u3
		q_global = q_global + R_3D.T.dot(q_local)

	print("q global is",q_global)
	return q_global


#############################################################################
def DoPCGM_LumpedTwoLevel(Ucor_ini,r_cor_o,nParts,nCG,nDim,tolc):

	print("coarse tolerance is",tolc)

	n = len(Ucor_ini)

	# print("shape of rho is",np.shape(rho))

	Ucor_updt = np.zeros([nCG*nDim])
	# print("shape of Uc update is ",np.shape(Ucor_updt))

	r_cor_n = np.zeros([nCG*nDim])

	P = np.zeros([nCG*nDim])

	z_n = np.zeros([nCG*nDim])

	print("r_cor_o inside Coarse PCGM is ",r_cor_o)

	z_o = GetPreconditionedZLumped(r_cor_o,nParts)

	print("preconditioned Z inside Coarse PCGM is",z_o)

	P = z_o

	rho_n = r_cor_o.T.dot(z_o)

	print("rho_n inside Coarse PCGM is",rho_n)

	Ucor_updt = Ucor_ini

	# print("Uc update is ",Ucor_updt)

	# print("P update is",P)

	z_n = z_o

	r_cor_n = r_cor_o

	for j in range(n):

		if(j==0):
			print("PCGM Lumped Iteration for Coarse Problem Starts")

		print("Iteration ID inside Coarse PCGM:",j+1)


		Q = GetMatrixProductLumpedTwoLevel(nParts,nCG,nDim,P)

		rho_temp = P.T.dot(Q)

		print("rho temp inside Coarse PCGM is",rho_temp)

		alpha = rho_n/rho_temp

########### Exit Criteria Relative Error #############################

		d2 = float(Ucor_updt.dot(Ucor_updt))

		print("d2 inside Coarse PCGM is",d2)

		error = (alpha**2 * P.dot(P))/d2

		print("Iteration number "+str(j+1)+" inside Coarse PCGM with Relative Error of :",error)

		if(error < tolc):
			break

#######################################################################


		Ucor_updt = Ucor_updt + (alpha * P)

		print("U cor update inside Coarse PCGM is",Ucor_updt)

		r_cor_n = r_cor_n - (alpha * Q)
		# print("r cor update is",r_cor_up)

		z_n = GetPreconditionedZLumped(r_cor_n,nParts)

		# print("Z update is",Z)

		rho_o = rho_n

		rho_n = r_cor_n.T.dot(z_n)

		beta = rho_n/rho_o

		print("beta inside Coarse PCGM is",beta)

		P = z_n + (beta * P)

#############   Norm Exit Criteria  ################

		# normval = linalg.norm(r_cor_n)/linalg.norm(r_cor_o)

		# print("Norm value is",normval)

		# if(normval < tolc):
		# 	break

#############################################

	print("U corner inside Coarse PCGM is",Ucor_updt)
	return Ucor_updt





#############################################################################


def GetMatrixProductLumpedTwoLevel(nParts,nCG,nDim,P):

	q_global = np.zeros([nDim*nCG])

	for ip in range(nParts):

		npart = ip+1

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		# print(np.shape(Bc_3D))

		one = 1
		subpath = "../data/Amats/subdom"+str(npart)

		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADirpath = subpath+"/ADir"+str(one)+".dat"
		ADir = np.genfromtxt(ADirpath)

		ADicpath = subpath+"/ADic"+str(one)+".dat"
		ADic = np.genfromtxt(ADicpath)

		ADrrpath = subpath+"/ADrr"+str(one)+".dat"
		ADrr = np.genfromtxt(ADrrpath)

		ADccpath = subpath+"/ADcc"+str(one)+".dat"
		ADcc = np.genfromtxt(ADccpath)

		ADrcpath = subpath+"/ADrc"+str(one)+".dat"
		ADrc = np.genfromtxt(ADrcpath)
		# print("shape of ADrc",np.shape(ADrc))

		# print(np.shape(ADir))
		ADri = ADir.T
		ADcr = ADrc.T
		ADci = ADic.T


		out2 = linalg.solve(ADii,ADir)
		srr = ADrr - (ADri.dot(out2))
		scr = ADcr - (ADci.dot(out2))

		out3 = linalg.solve(ADii,ADic)
		src = ADrc - (ADri.dot(out3))
		scc = ADcc - (ADci.dot(out3))


		p_s = Bc_3D.dot(P)
		u1 = src.dot(p_s)

		if(srr.ndim==0):
			u2 = u1/srr
		else:
			print("dimension of Srr inside MatrixProductLumped is",srr.ndim)
			u2 = linalg.solve(srr,u1)

		u3 = scr.dot(u2)
		u4 = scc.dot(p_s)
		q_local = u4 - u3
		q_global = q_global + Bc_3D.T.dot(q_local)

	# print("shape of q_global",np.shape(q_global))
	# print("q global",q_global)

	return q_global

###############################################################################

def DoDirectSolveCoarse(Ucor_ini,r_cor_o,nParts,nCG,nDim):

	Fcc = np.zeros([nDim*nCG])

	for ip in range(nParts):

		npart = ip+1

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		# print(np.shape(Bc_3D))

		one = 1
		subpath = "../data/Amats/subdom"+str(npart)

		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADirpath = subpath+"/ADir"+str(one)+".dat"
		ADir = np.genfromtxt(ADirpath)

		ADicpath = subpath+"/ADic"+str(one)+".dat"
		ADic = np.genfromtxt(ADicpath)

		ADrrpath = subpath+"/ADrr"+str(one)+".dat"
		ADrr = np.genfromtxt(ADrrpath)

		ADccpath = subpath+"/ADcc"+str(one)+".dat"
		ADcc = np.genfromtxt(ADccpath)

		ADrcpath = subpath+"/ADrc"+str(one)+".dat"
		ADrc = np.genfromtxt(ADrcpath)
		# print("shape of ADrc",np.shape(ADrc))

		# print(np.shape(ADir))
		ADri = ADir.T
		ADcr = ADrc.T
		ADci = ADic.T

		out2 = linalg.solve(ADii,ADir)
		srr = ADrr - (ADri.dot(out2))
		scr = ADcr - (ADci.dot(out2))

		out3 = linalg.solve(ADii,ADic)
		src = ADrc - (ADri.dot(out3))
		scc = ADcc - (ADci.dot(out3))

		temp1 = linalg.solve(srr,src)
		temp2 = scc - (scr.dot(temp1))
		temp3 = temp2.dot(Bc_3D)
		BcSc = Bc_3D.T.dot(temp3)
		Fcc = Fcc + BcSc

	print("Solving Corner Nodes Schur Complement Problem Directly")
	Ucor_ini = linalg.solve(Fcc,r_cor_o)

	return Ucor_ini

###############################################################################

def GetPreconditionedZLumped(r_cor_o,nParts):

	Zc = np.zeros_like(r_cor_o)

	for ip in range(nParts):

		npart = ip+1

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		print(np.shape(Bc_3D))

		one = 1
		subpath = "../data/Amats/subdom"+str(npart)
		ADccpath = subpath+"/ADcc"+str(one)+".dat"
		ADcc = np.genfromtxt(ADccpath)
		bcr = Bc_3D.dot(r_cor_o)
		# cr = ADcc.dot(bcr)
		cr = linalg.solve(ADcc,bcr)
		#print("cr inside preconditioned zlumped is ",cr)
		bcrT = Bc_3D.T.dot(cr)
		# print('shape of Zc is',np.shape(Zc))
		# print('value of Zc before addition is',Zc)
		# print('shape of bcrT and its value is',np.shape(bcrT))
		print(bcrT)
		Zc = Zc + bcrT
		# print('Zc is',Zc)

	# print("Zc is ",Zc)
	return Zc




################################################################################

def GetPreconditionedResidual_TwolevelWB(nParts,r_gamma_i,nDim,nCG,tolc):


	d_global = np.zeros([nDim*nCG])
	# print("shape of d global is",np.shape(d_global))

	for ip in range(nParts):

		npart = ip+1

		## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
		mdpath = "../../../data/meshData/meshdim"+str(npart).zfill(4)+".dat"
		meshdim = np.genfromtxt(mdpath)
		nNodes = meshdim[0].astype(int)     ## For 3D
		nBnodes = meshdim[4].astype(int)    ## For 3D

		mcrdpath = "../../../data/meshData/dimcrn"+str(npart).zfill(4)+".dat"
		meshdim = np.genfromtxt(mcrdpath)
		nCnodes = meshdim[0].astype(int)
		nRnodes = meshdim[1].astype(int)

		Rpath = "../../../data/meshData/R3D_"+str(npart)+".dat"
		R_3D = np.genfromtxt(Rpath)
		# print(np.shape(R_3D))

		BDpath = "../../../data/meshData/D3D_"+str(npart)+".dat"
		D_3D = np.genfromtxt(BDpath)
		# print(np.shape(D_3D))

		Rrpath = "../../../data/meshData/Rr3D_"+str(npart)+".dat"
		Rr_3D = np.genfromtxt(Rrpath)
		# print(np.shape(Rr_3D))

		Rcpath = "../../../data/meshData/Rc3D_"+str(npart)+".dat"
		Rc_3D = np.genfromtxt(Rcpath)
		# print(np.shape(Rc_3D))

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		# print(np.shape(Bc_3D))


		one = 1
		subpath = "../data/Amats/subdom"+str(npart)
		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADirpath = subpath+"/ADir"+str(one)+".dat"
		ADir = np.genfromtxt(ADirpath)

		ADicpath = subpath+"/ADic"+str(one)+".dat"
		ADic = np.genfromtxt(ADicpath)

		ADrrpath = subpath+"/ADrr"+str(one)+".dat"
		ADrr = np.genfromtxt(ADrrpath)

		ADccpath = subpath+"/ADcc"+str(one)+".dat"
		ADcc = np.genfromtxt(ADccpath)

		ADrcpath = subpath+"/ADrc"+str(one)+".dat"
		ADrc = np.genfromtxt(ADrcpath)

		ADri = ADir.T
		ADcr = ADrc.T
		ADci = ADic.T

####### Initial Residue for Whole Interface Problem ##########
		out2 = linalg.solve(ADii,ADir)
		print("shape of out2 is",np.shape(out2))
		srr = ADrr - (ADri.dot(out2))
		print("srr is ",srr)
		scr = ADcr - (ADci.dot(out2))
		print("dim of srr is",srr.ndim)

		# print("r_gamma_i is ",r_gamma_i)
		rg_local = D_3D.dot(R_3D.dot(r_gamma_i))

		f_r = Rr_3D.dot(rg_local)
		f_c = Rc_3D.dot(rg_local)
		# print("fc is",f_c)

		if(srr.ndim==0):
			v1 = f_r/srr
		else:
			v1 = linalg.solve(srr,f_r)
		# print("v1 is",v1)
		d_local = f_c - scr.dot(v1)
		# print("d_local is ",d_local)
		# plt.spy(Bc_3D)
		# plt.show()
		d_temp = Bc_3D.T.dot(d_local)
		# print("shape of d temp is",np.shape(d_temp))
		# print("d temp is",d_temp)
		d_global = d_global + Bc_3D.T.dot(d_local)
		# print("dc global is",d_global)

		# dcglobalpath = "../../../data/meshData/dcglobal_"+str(npart)+".dat"
		# np.savetxt(dcglobalpath, d_global)

######## Solve "Fcc Zc = dc " using Lumped Preconditioner #######

	Ucor_ini = np.zeros_like(d_global)
	r_cor_o = d_global
	Zc = DoPCGM_LumpedTwoLevel(Ucor_ini,r_cor_o,nParts,nCG,nDim,tolc)
	# Zc = DoDirectSolveCoarse(Ucor_ini,r_cor_o,nParts,nCG,nDim)
	print("Zc is",Zc)

##############Coarse Problem Solved ##########################
#############Remianing/Face Nodes of WB #####################

	Z_global = np.zeros_like(r_gamma_i)
	print("shape of Z global", np.shape(Z_global))

	for ip in range(nParts):

		npart = ip+1

		Rpath = "../../../data/meshData/R3D_"+str(npart)+".dat"
		R_3D = np.genfromtxt(Rpath)
		# print(np.shape(R_3D))

		BDpath = "../../../data/meshData/D3D_"+str(npart)+".dat"
		D_3D = np.genfromtxt(BDpath)
		# print(np.shape(D_3D))

		Rrpath = "../../../data/meshData/Rr3D_"+str(npart)+".dat"
		Rr_3D = np.genfromtxt(Rrpath)
		# print(np.shape(Rr_3D))

		Rcpath = "../../../data/meshData/Rc3D_"+str(npart)+".dat"
		Rc_3D = np.genfromtxt(Rcpath)
		# print(np.shape(Rc_3D))

		Bcpath = "../../../data/meshData/Bc3D_"+str(npart)+".dat"
		Bc_3D = np.genfromtxt(Bcpath)
		# print(np.shape(Bc_3D))

		one = 1
		subpath = "../data/Amats/subdom"+str(npart)
		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADirpath = subpath+"/ADir"+str(one)+".dat"
		ADir = np.genfromtxt(ADirpath)

		ADicpath = subpath+"/ADic"+str(one)+".dat"
		ADic = np.genfromtxt(ADicpath)

		ADrrpath = subpath+"/ADrr"+str(one)+".dat"
		ADrr = np.genfromtxt(ADrrpath)

		ADccpath = subpath+"/ADcc"+str(one)+".dat"
		ADcc = np.genfromtxt(ADccpath)

		ADrcpath = subpath+"/ADrc"+str(one)+".dat"
		ADrc = np.genfromtxt(ADrcpath)

		# print(np.shape(ADir))
		ADri = ADir.T
		ADcr = ADrc.T
		ADci = ADic.T


		out2 = linalg.solve(ADii,ADir)
		srr = ADrr - (ADri.dot(out2))

		out3 = linalg.solve(ADii,ADic)
		src = ADrc - (ADri.dot(out3))

		zc_local = Bc_3D.dot(Zc)
		rg_local = D_3D.dot(R_3D.dot(r_gamma_i))

		f_r = Rr_3D.dot(rg_local)
		v2 = f_r - (src.dot(zc_local))

		if(srr.ndim==0):
			zr_local = v2/srr
		else:
			zr_local = linalg.solve(srr,v2)

		z_local = Rr_3D.T.dot(zr_local) + Rc_3D.T.dot(zc_local)
		Z_global = Z_global + R_3D.T.dot(D_3D.dot(z_local))


	print("Z_global is",Z_global)
	return Z_global






############################################################################





###############################################################################

def GetInitialResidue(rbg,nParts):


	for ip in range(nParts):

		npart = ip+1

		## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
		mdpath = "../../../data/meshData/meshdim"+str(npart).zfill(4)+".dat"
		meshdim = np.genfromtxt(mdpath)
		nNodes = meshdim[0].astype(int)     ## For 3D
		nBnodes = meshdim[4].astype(int)    ## For 3D

		mcrdpath = "../../../data/meshData/dimcrn"+str(npart).zfill(4)+".dat"
		meshdim = np.genfromtxt(mcrdpath)
		nCnodes = meshdim[0].astype(int)
		nRnodes = meshdim[1].astype(int)

		Rpath = "../../../data/meshData/R3D_"+str(npart)+".dat"
		R_3D = np.genfromtxt(Rpath)
		print(np.shape(R_3D))


		one = 1
		subpath = "../data/Amats/subdom"+str(npart)
		ADiipath = subpath+"/ADii"+str(one)+".dat"
		ADii = np.genfromtxt(ADiipath)

		ADggpath = subpath+"/ADgg"+str(one)+".dat"
		ADgg = np.genfromtxt(ADggpath)

		ADigpath = subpath+"/ADig"+str(one)+".dat"
		ADig = np.genfromtxt(ADigpath)

		bipath = subpath+"/bi"+str(one)+".dat"
		bi = np.genfromtxt(bipath)

		bgpath = subpath+"/bg"+str(one)+".dat"
		bg = np.genfromtxt(bgpath)

		########## One Level Schur Complement and Residual ##############

		ADgi = ADig.T

		out2 = linalg.solve(ADii,bi)
		Gg = ADgi.dot(out2)
		Gg = bg - Gg
		RGi = R_3D.T.dot(Gg)
		rbg = rbg + RGi
		# print("\n")
		# print(ip)
		# print(rbg)

	return rbg


# def GetPreconditionedResidual_TwolevelNN(n_part,r_gamma_o):







##############################################################################

       ##             MAIN CODE  ###

###############################################################################
t1 = time.time()

tol = 1e-8
tolc = 1e-5

## Path to the global meshdims.dat
meshdimpath = '../../../data/meshData/meshdim.dat'
meshdim = np.genfromtxt(meshdimpath)
nPG = meshdim[0].astype(int)     ## Total Number of Points
nParts = meshdim[5].astype(int)   ## Total Number of Partitions

meshdimpath = '../../../data/meshData/boundary_nodes.dat'
Gbfile = np.genfromtxt(meshdimpath)
nodesBG = Gbfile.astype(int)
nBG = len(nodesBG)

meshdimpath = '../../../data/meshData/corner_nodes.dat'
CGfile = np.genfromtxt(meshdimpath)
nodesCG = CGfile.astype(int)
nCG = len(nodesCG)

meshdimpath = '../../../data/meshData/remaining_nodes.dat'
RGfile = np.genfromtxt(meshdimpath)
nodesRG = RGfile.astype(int)
nRG = len(nodesRG)

pdeDimPath = '../../../data/meshData/PDE_Dim.txt'
pdedim = np.genfromtxt(pdeDimPath)
nDim = pdedim.astype(int)


### Getting the Initial residue / One Level Schur Complement ######
res_in = np.zeros([nDim*nBG])
g_global = GetInitialResidue(res_in,nParts)

####### PCGM Iteration Two Level Neumann Neumann Preconditioner #######

U_gamma = DoPCGM_TwoLevel_WB(nDim,nParts,g_global,nBG,nCG,nRG,tol,tolc)

subpath = '../../../data/solution/'
Ugamma_path = subpath+'interface_solu.dat'
np.savetxt(Ugamma_path,U_gamma)

############## Interior Nodes Solution ######################

U = np.zeros([nPG*nDim,1])
print("Shape of U",np.shape(U))

for ip in range(nParts):

	npart = ip+1

	mdpath = "../../../data/meshData/meshdim"+str(npart).zfill(4)+".dat"
	meshdim = np.genfromtxt(mdpath)
	nP = meshdim[0].astype(int)     ## For 3D
	nB = meshdim[4].astype(int)    ## For 3D

	nodespath = "../../../data/meshData/nodes"+str(npart).zfill(4)+".dat"
	nodei = np.genfromtxt(nodespath) # Nodes in the subdomain interior at top
	nodei = nodei.astype(int)

	nI = nP - nB # Number of interior nodes in the sub domain

	one = 1
	subpath = "../data/Amats/subdom"+str(npart)
	ADiipath = subpath+"/ADii"+str(one)+".dat"
	ADii = np.genfromtxt(ADiipath)

	ADigpath = subpath+"/ADig"+str(one)+".dat"
	ADig = np.genfromtxt(ADigpath)

	bipath = subpath+"/bi"+str(one)+".dat"
	bi = np.genfromtxt(bipath)

	Rpath = "../../../data/meshData/R3D_"+str(npart)+".dat"
	R_3D = np.genfromtxt(Rpath)

	Ub = R_3D.dot(U_gamma)
	out1 = bi - ADig.dot(Ub)
	Ui = linalg.solve(ADii,out1)
	# print("Interior Solution",Ui)
	npi = len(Ui)

	######### Re arranging the nodes to form Global Solution #########

	for k in range(npi):

		if(k < nI):

			nni = nodei[k]
			# print("nni",nni)
			U[nni-1,:] = Ui[k]

		elif(k < 2*nI):

			nni = nodei[k-nI]
			# print("nni",nni)
			U[(nPG-1)+nni,:] = Ui[k]

		else:
			nni = nodei[k-(2*nI)]
			# print("nni",nni)
			U[(2*nPG-1)+nni,:] = Ui[k]


subpath = '../../../data/solution/'
Uinterior_path = subpath+'interior_solu.dat'
np.savetxt(Uinterior_path,U)
# print("U interior is",U)


npB = len(U_gamma)
# print("U_gamma is",U_gamma)
# print("U before putting gamma nodes is",U)

for j in range(npB):

	if(j < nBG):

		tempG = nodesBG[j]
		# print("tempG",tempG)
		U[tempG-1,:] = U_gamma[j]

	elif(j < 2*nBG):

		tempG = nPG + nodesBG[j-nBG]
		# print("tempG",tempG)
		U[tempG-1,:] = U_gamma[j]

	else:
		tempG = 2*nPG + nodesBG[j-(2*nBG)]
		# print("tempG",tempG)
		U[tempG-1,:] = U_gamma[j]

subpath = '../../../data/solution/'
U_path = subpath+'final_solu.dat'
np.savetxt(U_path,U)
# print("U final is",U)

### Timing ####

t2 = time.time()

print("Execution Time for processor is",t2-t1)

















































