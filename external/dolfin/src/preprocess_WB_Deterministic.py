


######### FUNCTIONS to GET RESTRICTION MATRICES ###########


#from petsc4py import PETSc
import numpy as np
import scipy.io as sio
import os as os
import matplotlib.pylab as plt
import scipy.sparse as sparse

   ## For 3D


# for ip in range(nParts):
#   ## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data


def GetR(nBG,nodesBG,nB,nodesB,nPart):

	# Taking the boundary nodes
	Rmat = np.zeros([nB,nBG])

	for i in range(nB):


		Ri = nodesB[i]

		for j in range(nBG):

			Rj = nodesBG[j]

			if Ri == Rj:
				Rmat[i,j] = 1
			else:
				Rmat[i,j] = 0

	return (Rmat.astype(int))

################################

def GetRr(nB,nodesB,nR,nodesR,nPart):
	# Taking the boundary nodes

	Rrmat = np.zeros([nR,nB])

	for i in range(nR):

		if nR==1:
			Ri = nodesR
		else:
			Ri = nodesR[i]


		for j in range(nB):

			Rj = nodesB[j]

			if Ri == Rj:
				Rrmat[i,j] = 1
			else:
				Rrmat[i,j] = 0

	return Rrmat.astype(int)

################################

def GetRc(nB,nodesB,nC,nodesC,nPart):
	# Taking the boundary nodes

	Rcmat = np.zeros([nC,nB])

	for i in range(nC):

		if nC == 1:
			Ri = nodesC
		else:
			Ri = nodesC[i]

		for j in range(nB):

			Rj = nodesB[j]

			if Ri == Rj:
				Rcmat[i,j] = 1
			else:
				Rcmat[i,j] = 0

	return Rcmat.astype(int)

################################

def GetBc(nCG,nodesCG,nC,nodesC,nPart):
	# Taking the boundary nodes

	Bcmat = np.zeros([nC,nCG])

	for i in range(nC):

		if nC == 1:
			Ri = nodesC
		else:
			Ri = nodesC[i]

		for j in range(nCG):

			Rj = nodesCG[j]

			if Ri == Rj:
				Bcmat[i,j] = 1
			else:
				Bcmat[i,j] = 0

	return Bcmat.astype(int)

################################

def GetBD(nBG,nodesBG,nB,nodesB,nParts):

	count = np.zeros([nBG,1])

	print("Global nodes inside BD",nodesBG)
	for k in range(1,nParts+1):

		mcrdpath = "../../../data/meshData/nbnodes"+str(k).zfill(4)+".dat"
		nbfile = np.genfromtxt(mcrdpath)
		nodesB_local = nbfile.astype(int)
		nB_local = len(nodesB_local)
		# print("local nodes inside BD from",nodesB_local,mcrdpath)

		for i in range(nBG):

			Ri = nodesBG[i]

			for j in range(nB_local):

				Rj = nodesB_local[j]

				if Ri == Rj:
					count[i] = count[i]+1

	# print(count)

	Dmat = np.zeros([nB,nB])

	for i in range(nB):

		Ri_2 = nodesB[i]

		for j in range(nBG):

			Rj_2 = nodesBG[j]

			if Ri_2 == Rj_2:

				Dmat[i,i] = 1/count[j]

	print("d is",Dmat)
	return Dmat


##################################################








#########################################################
#............MAIN CODE ......................
########################################################

pdeDimPath = '../../../data/meshData/PDE_Dim.txt'
pdedim = np.genfromtxt(pdeDimPath)
nDim = pdedim.astype(int)
print("nDim is",nDim)

## Reading Number of subdomains
meshdimpath = '../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[5].astype(int)
print("Number of Subdomains",nParts)

#Reading Global Boundary nodes
meshdimpath = '../../../data/meshData/boundary_nodes.dat'
Gbfile = np.genfromtxt(meshdimpath)
nodesBG = Gbfile.astype(int)
nBG = len(nodesBG)
print("Number of Global Boundary Nodes",nBG)

# mcrdpath = "../../../data/meshData/dimcrn"+str(nPart).zfill(4)+".dat"
# meshdim = np.genfromtxt(mcrdpath)
# nCnodes = meshdim[0].astype(int)
# nRnodes = meshdim[1].astype(int)

####  Restriction Matrix 3D

for ip in range(nParts):

	nPart = ip+1

	mdpath = "../../../data/meshData/meshdim"+str(nPart).zfill(4)+".dat"
	meshdim = np.genfromtxt(mdpath)
	nB = meshdim[4].astype(int)    ## Number of Boundary Nodes
	print("Number of Subdomain "+str(nPart)+" Boundary Nodes",nB)


	mcrdpath = "../../../data/meshData/nbnodes"+str(nPart).zfill(4)+".dat"
	print("Taking boundary nodes from ",mcrdpath)
	nbfile = np.genfromtxt(mcrdpath)
	nodesB = nbfile.astype(int)

	R_part = GetR(nBG,nodesBG,nB,nodesB,nPart)

	Rpath = "../../../data/meshData/R_"+str(nPart)+".dat"
	print("Saving Restriction Matrix R in ", Rpath)
	np.savetxt(Rpath, R_part)
	# plt.figure(nPart)
	# plt.spy(R_part)
	# print(R_part[1:nB,1:nBG])
	# print(R_part)

	##### Creating 3D restriction Matrix for R
	R_3D = np.zeros([nB*nDim,nBG*nDim])
	print(np.shape(R_3D))
	for j in range(nDim):
		R_3D[j*nB:(j+1)*nB,j*nBG:(j+1)*nBG] = R_part[0:nB,0:nBG]
		print(np.shape(R_3D))

	# plt.figure(nPart*2)
	# plt.spy(R_3D)
	R3Dpath = "../../../data/meshData/R3D_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, R_3D in ", R3Dpath)
	np.savetxt(R3Dpath, R_3D)


##### Restriction Matrix Rr -Scatter operator from Subdomain level Boundary to Subdomain Level Remaining

	mcrdpath = "../../../data/meshData/rnodes"+str(nPart).zfill(4)+".dat"
	print("Taking remaining nodes from ",mcrdpath)
	nbfile = np.genfromtxt(mcrdpath)
	nodesR = nbfile.astype(int)
	nR = nodesR.size
	print("Number of Subdomain "+str(nPart)+" Remaining Nodes",nR)

	Rr_part = GetRr(nB,nodesB,nR,nodesR,nPart)

	Rrpath = "../../../data/meshData/Rr_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Rr in ", Rrpath)
	np.savetxt(Rpath, Rr_part)

	##### Creating 3D restriction Matrix for R

	Rr_3D = np.zeros([nR*nDim,nB*nDim])

	for j in range(nDim):
		Rr_3D[j*nR:(j+1)*nR,j*nB:(j+1)*nB] = Rr_part[0:nR,0:nB]


	Rr3Dpath = "../../../data/meshData/Rr3D_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Rr_3D in ", Rr3Dpath)
	np.savetxt(Rr3Dpath, Rr_3D)



##### Restriction Matrix Rc -Scatter operator from Subdomain level Boundary to Subdomain Level Corner

	mcrdpath = "../../../data/meshData/cnodes"+str(nPart).zfill(4)+".dat"
	print("Taking corner nodes from ",mcrdpath)
	nbfile = np.genfromtxt(mcrdpath)
	nodesC = nbfile.astype(int)
	nC = nodesC.size
	print("Number of Subdomain "+str(nPart)+" Corner Nodes",nB)
	# print("Corner Nodes",nodesC)

	Rc_part = GetRc(nB,nodesB,nC,nodesC,nPart)

	Rcpath = "../../../data/meshData/Rc_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Rc in ", Rcpath)
	np.savetxt(Rpath, Rc_part)

	##### Creating 3D restriction Matrix for Rc
	Rc_3D = np.zeros([nC*nDim,nB*nDim])

	for j in range(nDim):
		Rc_3D[j*nC:(j+1)*nC,j*nB:(j+1)*nB] = Rc_part[0:nC,0:nB]

	Rc3Dpath = "../../../data/meshData/Rc3D_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Rc_3D in ", Rc3Dpath)
	np.savetxt(Rc3Dpath, Rc_3D)
	# plt.spy(Rc_3D)
	# plt.show()

##### Restriction Matrix Bc -Scatter operator from Global Corner node to Subdomain Level Corner


#Reading Global Corner nodes
	meshdimpath = '../../../data/meshData/corner_nodes.dat'
	Gbfile = np.genfromtxt(meshdimpath)
	nodesCG = Gbfile.astype(int)
	nCG = len(nodesCG)
	print(nCG)

	Bc_part = GetBc(nCG,nodesCG,nC,nodesC,nPart)

	Bcpath = "../../../data/meshData/Bc_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Bc in ", Bcpath)
	np.savetxt(Bcpath, Bc_part)

	##### Creating 3D restriction Matrix for Rc
	Bc_3D = np.zeros([nC*nDim,nCG*nDim])

	for j in range(nDim):
		Bc_3D[j*nC:(j+1)*nC,j*nCG:(j+1)*nCG] = Bc_part[0:nC,0:nCG]

	Bc3Dpath = "../../../data/meshData/Bc3D_"+str(nPart)+".dat"
	print("Saving Restriction Matrix, Bc_3D in ", Bc3Dpath)
	np.savetxt(Bc3Dpath, Bc_3D)

##### Diagonal Scaling Matrix D - Global Boundary to Subdomain Boundary

	BD_part = GetBD(nBG,nodesBG,nB,nodesB,nParts)
	print("BD_part is",BD_part)

	BDpath = "../../../data/meshData/D_"+str(nPart)+".dat"
	print("Saving Scaling Matrix, BD_part in ", BDpath)
	np.savetxt(BDpath, BD_part)

##### Creating 3D Scaling Matrix for BD

	BD_3D = np.zeros([nB*nDim,nB*nDim])

	for j in range(nDim):
		BD_3D[j*nB:(j+1)*nB,j*nB:(j+1)*nB] = BD_part[0:nB,0:nB]

	BD3Dpath = "../../../data/meshData/D3D_"+str(nPart)+".dat"
	print("Saving Sacling Matrix, BD_3D in ", BD3Dpath)
	np.savetxt(BD3Dpath, BD_3D)
	# plt.spy(BD_3D)
	# plt.show()
























