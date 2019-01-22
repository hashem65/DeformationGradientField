'''
Created on 23/11/2018

@author: rjag008
'''
from __future__ import print_function
import numpy as np
from TubeGenerator import Tube3d
from collections import OrderedDict
from scipy.linalg import logm


class GaussQuadarture(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.gaussPointLocations = np.zeros((4,4))
        self.gaussPointWeights = np.zeros((4,4))
        
        self.gaussPointLocations[0,0] = 0.5
        self.gaussPointWeights[0,0] = 1.0
        
        self.gaussPointLocations[1,0] = (-1.0/np.sqrt(3.0)+1.0)/2.0 
        self.gaussPointLocations[1,1] = (1.0/np.sqrt(3.0)+1.0)/2.0
        self.gaussPointWeights[1,0] = 0.5
        self.gaussPointWeights[1,1] = 0.5
        
        self.gaussPointLocations[2,0] = (-np.sqrt(0.6)+1.0)/2.0
        self.gaussPointLocations[2,1] = 0.5
        self.gaussPointLocations[2,2] = (np.sqrt(0.6)+1.0)/2.0
        self.gaussPointWeights[2,0] = 5.0/18.0        
        self.gaussPointWeights[2,1] = 4.0/9.0
        self.gaussPointWeights[2,2] = 5.0/18.0      
        
        self.gaussPointLocations[3,0] = (-np.sqrt((3.0+2.0*np.sqrt(6.0/5.0))/7.0)+1.0)/2.0
        self.gaussPointLocations[3,1] = (-np.sqrt((3.0-2.0*np.sqrt(6.0/5.0))/7.0)+1.0)/2.0
        self.gaussPointLocations[3,2] = (+np.sqrt((3.0-2.0*np.sqrt(6.0/5.0))/7.0)+1.0)/2.0
        self.gaussPointLocations[3,3] = (+np.sqrt((3.0+2.0*np.sqrt(6.0/5.0))/7.0)+1.0)/2.0
        self.gaussPointWeights[3,0] = (18.0-np.sqrt(30.0))/72.0              
        self.gaussPointWeights[3,1] = (18.0+np.sqrt(30.0))/72.0       
        self.gaussPointWeights[3,2] = (18.0+np.sqrt(30.0))/72.0
        self.gaussPointWeights[3,3] = (18.0-np.sqrt(30.0))/72.0    
        
    def getQuadratureFeaturesFor(self,num):
        if num < 1 or num > 4:
            raise ValueError('Not implemented for xi %d'%num)
        return [self.gaussPointLocations[num-1,:num],self.gaussPointWeights[num-1,:num]]

class GenerateVolumeData(object):
    '''
    Load the exregion files at specified times
    Compute the deformation tensor at the gauss points for each element
    '''


    def __init__(self, exfiles,elementsAlongX1=16,elementsAlongX2=16,elementsAlongX3=1,deltax=0.001,outputMeshes=False,numGauss=3):
        '''
        exfiles - dict with timexexfile
        elementsAlong* Specify the number of elements along that XI - used to calculate the normalized coordinate
        deltax - neighbor radius in xi units
        '''
        #Load coordinate data
        #Calculate the number of time values
        meshFiles = dict(exfiles)
        if 'reference' in meshFiles:
            del meshFiles['reference']
        tubeGenerator = Tube3d(elementsAlongX1,elementsAlongX2,elementsAlongX3,len(meshFiles))
        region = tubeGenerator.getRegion()
        tubeGenerator.generateMesh(region)
        self.context = tubeGenerator.context
        for t,exnode in meshFiles.items():
            tubeGenerator.loadStateFromFile(exnode[0], t, coordinateFieldName=exnode[1],nodeOffset=int(exnode[2]))
        self.times = np.array(sorted(meshFiles.keys()))
        #Load the reference coordinates
        if 'reference' in exfiles:
            exnode = exfiles['reference']
        else:
            exnode = exfiles[self.times[0]]
        tubeGenerator.loadReferenceCoordinatesFromFile(exnode[0], coordinateFieldName=exnode[1],nodeOffset=int(exnode[2]))
        if outputMeshes:
            for t in exfiles:
                tubeGenerator.writeToFileAtTime('tube%d.ex2'%t, t)
        
        self.deltax = deltax
        
        self.elements = OrderedDict()
        fieldModule = region.getFieldmodule()
        fieldModule.beginChange()        
        coordinatesField = fieldModule.findFieldByName('coordinates').castFiniteElement()
        
        self.coordinatesField = coordinatesField
        mesh = fieldModule.findMeshByDimension(3)
        ei   = mesh.createElementiterator()
        elem = ei.next()
        while elem.isValid():
            self.elements[elem.getIdentifier()-1] = elem
            elem = ei.next()
        self.numberOfElements = len(self.elements)
        qgen = GaussQuadarture()
        self.gaussPointLocations = qgen.getQuadratureFeaturesFor(numGauss)[0]
        self.numGaussPoints = self.gaussPointLocations.shape[0]**3
        self.numPointsPerGauss = 1 # (xi + number of neighbors(6) + self)
        #Create the derivative field
        xifield = tubeGenerator.referenceCoordinates
        self.dx = fieldModule.createFieldGradient(coordinatesField,xifield)
        self.volume = fieldModule.createFieldDeterminant(self.dx)

        self.xifield = xifield        
        self.ex1 = elementsAlongX1
        self.ex2 = elementsAlongX2
        self.ex3 = elementsAlongX3
        nx = self.ex1*self.gaussPointLocations.shape[0]
        ny = self.ex2*self.gaussPointLocations.shape[0]
        nz = self.ex3*self.gaussPointLocations.shape[0]
        self.volumes = np.zeros((nx,ny,nz))
        self.tubeGenerator = tubeGenerator

    def saveMeshAtTime(self,filename,t):
        self.tubeGenerator.writeToFileAtTime('%s.ex2'%filename,t)
        self.tubeGenerator.saveAsNumpy(filename, t)

    def computeRateForTime(self,t):
        if t > self.times.max() or t < self.times.min():
            raise ValueError('Requested time %g lies out of bounds %g - %g'% (t,self.times.min(),self.times.max()))

        region = self.context.getDefaultRegion()
        fieldModule = region.getFieldmodule()
        fieldCache = fieldModule.createFieldcache()
        fieldCache.setTime(t)
        
        grid = self.volumes.shape
        rates = np.zeros((grid[0],grid[1],grid[2],3,3))
        #Assumes elements start at 0
        ngauss = self.gaussPointLocations.shape[0]
        xyn = self.ex1*self.ex2
        for eid in self.elements:
            elem = self.elements[eid]
            XV = (eid%self.ex1)*ngauss
            YV = (int(eid/self.ex1)%(self.ex2))*ngauss
            ZV = int(eid/xyn)*ngauss            
            
            for k,xi3 in enumerate(self.gaussPointLocations):
                for j,xi2 in enumerate(self.gaussPointLocations):
                    for i,xi1 in enumerate(self.gaussPointLocations):
                        fieldCache.setMeshLocation(elem,[xi1,xi2,xi3])
                        _,vol = self.dx.evaluateReal(fieldCache,9)                        
                        ev = logm(np.array(vol).reshape((3,3)))
                        rates[XV+i,YV+j,ZV+k] = ev
    
        return rates

    def computeDefGradForTime(self,t):
        if t > self.times.max() or t < self.times.min():
            raise ValueError('Requested time %g lies out of bounds %g - %g'% (t,self.times.min(),self.times.max()))

        region = self.context.getDefaultRegion()
        fieldModule = region.getFieldmodule()
        fieldCache = fieldModule.createFieldcache()
        fieldCache.setTime(t)
        
        grid = self.volumes.shape
        dgrad = np.zeros((grid[0],grid[1],grid[2],3,3))
        #Assumes elements start at 0
        ngauss = self.gaussPointLocations.shape[0]
        xyn = self.ex1*self.ex2
        for eid in self.elements:
            elem = self.elements[eid]
            XV = (eid%self.ex1)*ngauss
            YV = (int(eid/self.ex1)%(self.ex2))*ngauss
            ZV = int(eid/xyn)*ngauss            
            
            for k,xi3 in enumerate(self.gaussPointLocations):
                for j,xi2 in enumerate(self.gaussPointLocations):
                    for i,xi1 in enumerate(self.gaussPointLocations):
                        fieldCache.setMeshLocation(elem,[xi1,xi2,xi3])
                        _,vol = self.dx.evaluateReal(fieldCache,9)                        
                        dgrad[XV+i,YV+j,ZV+k] = np.array(vol).reshape((3,3))
    
        return dgrad

        
    def computeVolForTime(self,t):
        if t > self.times.max() or t < self.times.min():
            raise ValueError('Requested time %g lies out of bounds %g - %g'% (t,self.times.min(),self.times.max()))
        
        region = self.context.getDefaultRegion()
        fieldModule = region.getFieldmodule()
        fieldCache = fieldModule.createFieldcache()
        fieldCache.setTime(t)
        #Load and solve volumes        
        #Assumes elements start at 0
        ngauss = self.gaussPointLocations.shape[0]
        xyn = self.ex1*self.ex2
        for eid in self.elements:
            elem = self.elements[eid]
            XV = (eid%self.ex1)*ngauss
            YV = (int(eid/self.ex1)%(self.ex2))*ngauss
            ZV = int(eid/xyn)*ngauss            
            
            for k,xi3 in enumerate(self.gaussPointLocations):
                for j,xi2 in enumerate(self.gaussPointLocations):
                    for i,xi1 in enumerate(self.gaussPointLocations):
                        fieldCache.setMeshLocation(elem,[xi1,xi2,xi3])
                        _,vol = self.volume.evaluateReal(fieldCache,1)
                        self.volumes[XV+i,YV+j,ZV+k] = vol


def generateStageData(files,storeloc,ex1,ex2,ex3,discret=250):
    obj = GenerateVolumeData(files,ex1,ex2,ex3,numGauss=3)
    fks = list(sorted(files.keys()))
    times = np.linspace(fks[0], fks[-1], discret)
    obj.computeVolForTime(0.0)
    
    vol = obj.volumes.flatten()
        
    for i in times[1:]:
        obj.computeVolForTime(i)
        nv = obj.volumes.flatten()
        vol = np.c_[vol,nv]
    
    nvol = (vol[:,1:]-vol[:,:-1])/vol[:,:-1]
    vx = np.max(np.fabs(nvol),axis=0)
    
    #Find time slots by breaking when derivative jumps
    dvx = np.fabs(np.diff(vx))
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(dvx, height=0)
    
    import matplotlib.pyplot as plt
    plt.plot(dvx)
    plt.plot(peaks, dvx[peaks], "x")
    #plt.scatter(vx[1:],dvx)
    plt.show()
    
    print(times[peaks])
    
    #Peaks correspond points that we be intermediate stages
    stagetimes = [times[0]]
    stagetimes.extend(times[peaks].tolist())
    stagetimes.append(times[-1])
    #Find midpoints between peaks
    stagetimes = np.array(stagetimes)
    dist = (stagetimes[1:]-stagetimes[:-1])*0.9
    mstages = stagetimes[:-1]+dist
    ilace = np.c_[stagetimes[:-1],mstages]
    stagetimes = ilace.flatten()
    stagetimes = np.concatenate([stagetimes,[times[-1]]])
    np.savetxt('%s/stagetimes.csv'%storeloc,stagetimes,delimiter=',')
    
    with open('%s/vis.cmgui'%storeloc,'w') as ser:
        for t in stagetimes:
            print("Saving file %s/stage%0.2f"%(storeloc,t))
            obj.saveMeshAtTime('%s/stage%0.2f'%(storeloc,t), t)
            print("gfx read node stage%0.2f.ex2 time %f"%(t,t),file=ser)

def generateStageRates(ex1,ex2,ex3,storeloc):
    times = np.loadtxt('%s/stagetimes.csv'%storeloc, dtype=np.float, delimiter=',')
    for i in range(1,times.shape[0]):
        files = {0:['%s/stage%0.2f.ex2'%(storeloc,times[i-1]),'coordinates',0],
                 1:['%s/stage%0.2f.ex2'%(storeloc,times[i]),'coordinates',0]}
        obj = GenerateVolumeData(files,ex1,ex2,ex3,numGauss=3)
        rates = obj.computeRateForTime(1.0)
        print ('############################')
        print ('############################')
        print ('############################')
        print (rates)
		
        np.save('%s/rates%0.2f.pkl'%(storeloc,times[i]),rates)
        print("Saving rates %s/rates%0.2f"%(storeloc,times[i]))

def generateStageDefGrad(ex1,ex2,ex3,storeloc):
    times = np.loadtxt('%s/stagetimes.csv'%storeloc, dtype=np.float, delimiter=',')
    for i in range(1,times.shape[0]):
        files = {0:['%s/stage%0.2f.ex2'%(storeloc,times[i-1]),'coordinates',0],
                 1:['%s/stage%0.2f.ex2'%(storeloc,times[i]),'coordinates',0]}
        obj = GenerateVolumeData(files,ex1,ex2,ex3,numGauss=3)
        rates = obj.computeDefGradForTime(1.0)
        np.save('%s/defgrad%0.2f.pkl'%(storeloc,times[i]),rates)
        print("Saving defgrad %s/defgrad%0.2f"%(storeloc,times[i]))    
    
if __name__ == '__main__':

    files = {0:[r'../meshfiles04/mesh6-8x8.part0.exnode','Geometry',0],\
             1:[r'../meshfiles04/mesh8-8x8.part0.exnode','Geometry',0]}
    
    generateStageData(files,'../meshfiles04',8,8,1,discret=2)
    #Load the stage data and compute rates
    generateStageRates(8,8,1,'../meshfiles04')
    generateStageDefGrad(8,8,1,'../meshfiles04')
    
    