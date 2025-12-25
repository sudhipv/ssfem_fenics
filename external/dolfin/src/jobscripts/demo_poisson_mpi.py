## Modified from FEniCS demo:
#
# Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
#
# Code to test MPI with FEniCS

import dolfin as d
from mpi4py import MPI

print("************* Running FEniCS with MPI... ***************")

comm = MPI.COMM_WORLD
ip = comm.Get_rank()

print('My rank is ',ip)

mesh = d.UnitSquareMesh(d.mpi_comm_self(), 2, 2)

V = d.VectorFunctionSpace(mesh, "Lagrange", 1)
gfun = d.Function(V)

print("========================================================")
print("======================== Success =======================")
print("========================================================")
