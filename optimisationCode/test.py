from __future__ import print_function
import numpy as np

# modules required for fitting 
#from __future__ import print_function
#import numpy as np
from opencmiss.iron import iron
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.context import Context


from FitGaussToNode import Fitting
	
kinematicRates = np.load('rates1.00.pkl.npy')

kinematicRatesDiag = np.zeros((12,9,3,3))
kinematicRatesDiag [:,:,:,0] = kinematicRates [:,:,:,0,0]
kinematicRatesDiag [:,:,:,1] = kinematicRates [:,:,:,1,1]
kinematicRatesDiag [:,:,:,2] = kinematicRates [:,:,:,2,2]

numberOfElements = 12 
growthRateDivided =  np.zeros((16,3))        #  x.reshape(16,3)            #    np.zeros((16,3)       #16*3
growthElems = np.zeros((numberOfElements,3,3,3,3))                        #12*3*3*3*3
elemId = np.zeros((numberOfElements,4)) 
elemId[0,:] =  0,1,4,5
elemId[1,:] =  1,2,5,6
elemId[2,:] =  2,3,6,7
elemId[3,:] =  3,0,7,4
elemId[4,:] =  4,5,8,9
elemId[5,:] =  5,6,9,10
elemId[6,:] =  6,7,10,11
elemId[7,:] =  7,4,11,8
elemId[8,:] =  8,9,12,13
elemId[9,:] =  9,10,13,14
elemId[10,:] =  10,11,14,15
elemId[11,:] =  11,8,15,12
for i in range (numberOfElements):
    growthElems [i,0,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.21)		
    growthElems [i,1,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.21)		
    growthElems [i,2,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.21)		

    growthElems [i,0,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.50)		
    growthElems [i,1,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.50)		
    growthElems [i,2,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.50)		

    growthElems [i,0,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.79)		
    growthElems [i,1,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.79)		
    growthElems [i,2,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.79)						
                    
print ('finished up to here ... ')
realRates = np.zeros ((12,9,3,3))
for circumElem in range (4):
    for longElem in range (3):
        for xi1Index in range (3):
            for xi2Index in range (3):
                for xi3Index in range (3):
                    realRates [circumElem*3+xi1Index,longElem*3+xi2Index,xi3Index,:]  = growthElems [numberOfElements-1,xi1Index,xi2Index,xi3Index,:]  
                
f2  = np.sum(np.linalg.norm(rts1-rts2,axis=1)) 

                
obj1 = Fitting(4,3,1)
obj1.setupProblem()
rts1 = obj.solveGrowthRates(realRates)
obj2 = Fitting(4,3,1)
obj2.setupProblem()
rts2 = obj.solveGrowthRates(kinematicRatesDiag)				
				
