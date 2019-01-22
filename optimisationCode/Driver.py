from __future__ import print_function
import numpy as np
import json,os,sys
import traceback
import datetime 
import math 


from FitGaussToNode import Fitting
from io import StringIO
from subprocess import Popen,PIPE
from sqlitedict import SqliteDict  
size = 1
myrank = 0
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    size = comm.Get_size()
except:
    raise ImportError('mpi4py is required for parallelization')

my_env = os.environ.copy()
'''
Even when launching iron through a subprocess,
MPI may find out that it was a part of a distributed process through environment and filesystem based flags
Fortunately, in the case of mpich, it has a few environmental flagsm which if removed the subprocess is loaded
with its own communicator!!

https://stackoverflow.com/questions/21090085/calling-mpi-binary-in-serial-as-subprocess-of-mpi-application

Here are the flags 
'''
try:
    del my_env['PMI_RANK']
except:
    pass
try:    
    del my_env['PMI_SIZE']
except:
    pass
try:    
    del my_env['PMI_FD']
except:
    pass

obj1 = Fitting(4,3,1)
obj1.setupProblem()
kinematicRates = np.load('rates1.00.pkl.npy')
#print (kinematicRates[0,0,0,0,0], kinematicRates[0,0,0,0,1], kinematicRates[0,0,0,0,2])
#print (kinematicRates[0,0,0,1,0],kinematicRates[0,0,0,1,1], kinematicRates[0,0,0,1,2])
#print (kinematicRates[0,0,0,2,0],kinematicRates[0,0,0,2,1], kinematicRates[0,0,0,2,2])
#print (kinematicRates.shape)

def findDeformedCoordinates(rates):
    try:
        p = Popen(['python','parallelOptimization.py'],stdin=PIPE,stderr=PIPE,env=my_env)
        inputData = dict()
        inputData['optmin'] = rates.tolist() 
        subin = json.dumps(inputData)
        #print('##Initializing problem for parameters %s ' % (' ,'.join(map(str,rates.tolist()))))
        _,perr = p.communicate(str.encode(subin))
        p.wait()
        perr = perr.decode()
        #Ensure that we are starting at a json prefix
        ps = perr.index('POSResSTART') + 11
        pe = perr.index('POSResEND')
        result = json.loads(perr[ps:pe])
        if u'error' in result:
            raise ValueError(result[u'error'])
        return np.array(result[u'nodePositions'])
    except:
        traceback.print_exc(file=sys.stdout)

dictionaryLoc = 'solutions431.sql'
solutionsDictionary = SqliteDict(dictionaryLoc, autocommit=True)

def getSolutionsDictionary():
    return solutionsDictionary

import random
import pygmo as pg

class GrowthOptimization(object):
    precision = 1e3
    bestAnswer = 1e16
    def __init__(self):
        growthRate = np.zeros((16,6))
        self.grshape = growthRate.shape
        growthRate[0,:3] =  0.05,0.08,0.14
        growthRate[1,:3] = 0.08,0.06,0.08
        growthRate[2,:3] = -0.13,-0.03,0.04
        growthRate[3,:3] = -0.02,-0.12,0.06	
		
        growthRate[4,:3] = 0.18,0.23,0.06
        growthRate[5,:3] = 0.19,0.14,0.05
        growthRate[6,:3] = -0.13,-0.03,0.04
        growthRate[7,:3] = -0.02,-0.12,0.06	

        growthRate[8,:3] = 0.18,0.23,0.06
        growthRate[9,:3] = 0.19,0.14,0.05
        growthRate[10,:3] =-0.13,-0.03,0.04
        growthRate[11,:3] =-0.02,-0.12,0.06	

        growthRate[12,:3] =  0.18,0.23,0.06
        growthRate[13,:3] =  0.19,0.14,0.05
        growthRate[14,:3] = -0.13,-0.03,0.04
        growthRate[15,:3] = -0.02,-0.12,0.06

        growthRate[0:4,3:]=0.0
        growthRate[4:8,3:]=0.0
        growthRate[8:12,3:]=0.0
        growthRate[12:,3:]=0.0
        self.targetCoordinates = findDeformedCoordinates(growthRate)
    
    def checkSolution(self,key):
        solutions = getSolutionsDictionary()
        mkey = (np.array(key)*self.precision).astype('int')
        with StringIO() as sfile:
            print(mkey, file=sfile)
            skey = sfile.getvalue()
            if skey in solutions:
                print("Found solution for ",key," ",solutions[skey])
                return solutions[skey]
        return None
    
    def addSolution(self,key,value):
        solutions = getSolutionsDictionary()
        mkey = (np.array(key)*self.precision).astype('int')
        with StringIO() as sfile:
            print(mkey, file=sfile)
            skey = sfile.getvalue()
            solutions[skey]=value
            
    def printSolutions(self):
        solutions = getSolutionsDictionary()
        for k,v in solutions.items():
            print(k," ",v)
                       
    def fitness(self,x):
        #f = self.checkSolution(x)
        #if f is None:
        #    try: 	



        kinematicRatesDiag = np.zeros((12,9,3,3))
        kinematicRatesDiag [:,:,:,0] = kinematicRates [:,:,:,0,0]
        kinematicRatesDiag [:,:,:,1] = kinematicRates [:,:,:,1,1]
        kinematicRatesDiag [:,:,:,2] = kinematicRates [:,:,:,2,2]

        growthRateDivided = x.reshape(16,3)            #    np.zeros((16,3)       #16*3
        growthElems = np.zeros((12,3,3,3,3))                        #12*3*3*3*3
        elemId = np.zeros((12,4)) 
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
        for i in range (12):
            growthElems [i,0,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.21)		
            growthElems [i,1,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.21)		
            growthElems [i,2,0,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.21)		

            growthElems [i,0,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.50)		
            growthElems [i,1,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.50)		
            growthElems [i,2,1,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.50)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.50)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.50)		

            growthElems [i,0,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.79)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.21)*(0.79)		
            growthElems [i,1,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.50)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.50)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.50)*(0.79)		
            growthElems [i,2,2,:,:] = growthRateDivided[int(elemId[i,0]),:]*(0.21)*(0.21)+ growthRateDivided[int(elemId[i,1]),:]*(0.79)*(0.21)+ growthRateDivided[int(elemId[i,2]),:]*(0.21)*(0.79)+ growthRateDivided[int(elemId[i,3]),:]*(0.79)*(0.79)						
                    
        numberOfElements = 12 
        realRates = np.zeros ((12,9,3,3))
        for circumElem in range (4):
            for longElem in range (3):
                for xi1Index in range (3):
                    for xi2Index in range (3):
                        for xi3Index in range (3):
                            realRates [circumElem*3+xi1Index,longElem*3+xi2Index,xi3Index,:]  = growthElems [numberOfElements-1,xi1Index,xi2Index,xi3Index,:]  
                
        realRatesOffZero = np.zeros ((12,9,3,3,3))
        realRatesOffZero [:,:,:,0,0] = realRates [:,:,:,0]
        realRatesOffZero [:,:,:,1,1] = realRates [:,:,:,1]
        realRatesOffZero [:,:,:,2,2] = realRates [:,:,:,2]

                
        #obj1 = Fitting(4,3,1)
        #obj1.setupProblem()
        rts1 = obj1.solveGrowthRates(kinematicRates)
        obj1.saveResults('kinematic')
        gradient1 = obj1.computeGradient('kinematic')
        kinematicRatesGrads = np.load('kinematicrates.npy')
        #print ('kinematicRatesGrads shape =', kinematicRatesGrads.shape)

        #obj2 = Fitting(4,3,1)
        #obj2.setupProblem()
        rts2 = obj1.solveGrowthRates(realRatesOffZero)				
        obj1.saveResults('iteration')
        gradient2 = obj1.computeGradient('iteration')
        iterationRatesGrads = np.load('iterationrates.npy')		
        #print ('iterationRatesGrads shape =', iterationRatesGrads.shape)


        f2  = np.sum(np.linalg.norm(kinematicRatesGrads-iterationRatesGrads)) 
        #print ('===============================================')
        #print ('===============================================')
        #print ('===============================================')
        #print ('===============================================')
        #print ('f2 = ', f2)




			
        growthFull = np.zeros((16,6))
        growthRatesMain = np.zeros((16,3))
        growthRatesCross = np.zeros((16,3))
        growthRatesMain[:,:] = x.reshape(16,3)
        growthFull[:,0:3] = growthRatesMain[:,:] 
        coordinates = findDeformedCoordinates(growthFull)
        f1 = np.sum(np.linalg.norm(coordinates-self.targetCoordinates,axis=1)) 
        f3  =  2e-9*f1
        f= 2e-9*f1+f2
        print ('f1 = ',2e-9*f1,'f2 = ',f2)
        self.addSolution(x, f)
        currentAnswer = f
        if (currentAnswer<0.97*self.bestAnswer):
            with open("answers.txt","a") as fileAnswer:
                fileAnswer.write ('  f1=    %f  ' %f1)
                fileAnswer.write (' weighted  f1=    %f  ' %f3)
                fileAnswer.write ('   f2=    %f  ' %f2)
                fileAnswer.write ('  totalObjctive=    %f \r\n' %currentAnswer)
                now = datetime.datetime.now()
                datetime.time(now.hour, now.minute, now.second)
                fileAnswer.write ('  hour %f  ' %now.hour)
                fileAnswer.write ('minute %f  ' %now.minute)
                np.savetxt(fileAnswer, x, newline ='\n')
                print ('Rate   ',x,'    objective',f)
                fileAnswer.write ('    *************************	\n \n ')
                fileAnswer.close()
            self.bestAnswer = currentAnswer
            #except:
                #f = 1e12	
        if (f < 1e-2):
            print ('Rate ',x,' objective ',f)
            sys.exit()
        return [f]

    def get_bounds(self):
        return ([-0.2]*48,[0.25]*48)
    
runParallel = False
import random
class MonteCarlo(pg.algorithm):
    """
     Monte-Carlo (random sampling) algorithm implemented purely in Python.
    """

    def __init__(self,iterations = 10):
        """
            Constructs a Monte-Carlo (random sampling) algorithm

             USAGE: algorithm.my_algorithm(iter = 10)

             NOTE: At the end of each iteration, the randomly generated
                     point substitutes the worst individual in the population if better

             * iter: number of random samples
        """
        #We start calling the base constructor
        super(MonteCarlo,self).__init__()
        #We then define the algorithm 'private' data members
        self.__iter = iterations

    #This is the 'juice' of the algorithm, the method where the actual optimzation is coded.
    def evolve(self,pop):
        #If the population is empty (i.e. no individuals) nothing happens
        if len(pop) == 0:
            return pop

        #Here we rename some variables, in particular the problem
        prob = pop.problem
        #Its dimensions (total and continuous)
        dim, cont_dim = prob.get_nx(), prob.get_ncx()
        #And the lower/upper bounds for the chromosome
        lb, ub = prob.get_bounds()
        
        #The algorithm now starts manipulating the population
        for _ in range(self.__iter):
            #We create a random vector within the bounds ... first the continuous part
            tmp_cont = [random.uniform(lb[i],ub[i]) for i in range(cont_dim)]
            #then the integer part
            tmp_int = [float(random.randint(lb[i],ub[i])) for i in range(cont_dim,dim)]
            #and we assemble them into one decision vector
            tmp_x = tmp_cont + tmp_int
            #which we push back in the population
            pop.push_back(tmp_x)
            #to then remove the worst individual
            #pop.erase(pop.get_worst_idx())
            #at the end of it all we return the 'evolved' population
            return pop

    def get_name(self):
        return "Monte Carlo (Python)"

    def human_readable_extra(self):
        return "n_iter=" + str(self.__n_iter)

if __name__ == '__main__':
    gp = GrowthOptimization()
    prob = pg.problem(gp)
    algo = pg.algorithm(pg.pso(gen = 300))
    #algo = MonteCarlo(100)
    try:
        if not runParallel:
            pop = pg.population(prob,40)
            print(dir(pop))
            pop = algo.evolve(pop)
            
            print(pop.champion_f) 
        else:
            archi = pg.archipelago(n = 16, algo = algo, prob = prob, pop_size = 20, seed = 32)
            archi.evolve()
            print(archi)
            archi.wait()
            res = archi.get_champions_f()
            print(res) 
    except:
        traceback.print_exc(file=sys.stdout)
    