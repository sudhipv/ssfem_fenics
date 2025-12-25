## PREPROCESSING
To change mesh density
$ vi ../foo3D.geo
Use GMSH to partition and create 3D decomposed mesh@../inputs/gmsh3D.geo


### Development: ONE-LEVEL DDM Data (interior, interface)
Go to 3D_meshDecomposer/ then from Matlab:
$ globalDecomposer3D.m   :: Generate global (p,e,t,h,b) data
$ localDecomposer3D.m     :: Generate local (for each subdomain) (p,e,t,h,b) data
$ oneLevelRenum3D.m       :: Re-arrange points->rpoints, nodes->Nnodes(i,b), etc.


### NoteStarted: TWO-LEVEL DDM Data (interior, interface[remaining, corner])
Go to 3D_meshDecomposer/ then from Matlab:
$ globalDecomposer3D.m   :: Generate global (p,e,t,h,b) data
$ localDecomposer3D.m     :: Generate local (for each subdomain) (p,e,t,h,b) data
$ cornerRemaingNodes.m   :: Generate local corner/remaining nodes
$ twoLevelRenum3D.m       :: Re-arrange points->rpoints, nodes->Nnodes(i,r,c(v)) etc.


### NotStarted: Wire-Basket DDM Data (interior, interface[remaining, corner])
Go to 3D_meshDecomposer/ then from Matlab:
$ globalDecomposer3D.m   :: Generate global (p,e,t,h,b) data
$ localDecomposer3D.m     :: Generate local (for each subdomain) (p,e,t,h,b) data
$ wireBasketNodes.m         :: Generate local corner/remaining nodes
$ twoLevelRenum3D.m       :: Re-arrange points->rpoints, nodes->Nnodes(i,r(f),c(v+e),) etc.
