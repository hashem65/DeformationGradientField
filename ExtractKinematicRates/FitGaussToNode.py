'''
Created on 2/12/2018

@author: rjag008
'''
from __future__ import print_function
import numpy as np
from opencmiss.iron import iron
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.context import Context

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

class Fitting(object):
    '''
    classdocs
    '''
    def __init__(self, circumferentialElements,axialElements,wallElements,discret=10,length=10.0,innerRadius=1.25,outerRadius=2.0):
        '''
        Constructor
        '''
        self.circumferentialElements = circumferentialElements
        self.axialElements = axialElements
        self.wallElements = wallElements
        self.length=length
        self.innerRadius=innerRadius
        self.outerRadius=outerRadius
        self.numPoints = discret**3

    def setupGeometry(self,geometricField):
        numberOfCircumfrentialElements = self.circumferentialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumfrentialNodes = numberOfCircumfrentialElements
        numberOfWallNodes = self.wallElements+1
            
        # Create the geometric field
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                        (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    radius = self.innerRadius + (self.outerRadius - self.innerRadius)*float(wallNodeIdx-1)/float(numberOfWallNodes)
                    theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*np.pi
                    x = radius*np.cos(theta)
                    y = radius*np.sin(theta)
                    xtangent = -np.sin(theta)
                    ytangent = np.cos(theta)
                    xnormal = np.cos(theta)
                    ynormal = np.sin(theta)
                    z = float(lengthNodeIdx-1)/float(numberOfLengthNodes)*self.length
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)

    def setupProblem(self,showProgress=False):
        # Number of Gauss points used
        numberOfGaussXi = 3
        numberOfCircumfrentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumfrentialNodes = numberOfCircumfrentialElements
        numberOfWallNodes = self.wallElements+1
        numberOfWallElements = self.wallElements
        
        coordinateSystemUserNumber = 1
        regionUserNumber = 1
        tricubicHermiteBasisUserNumber = 1
        meshUserNumber = 1
        decompositionUserNumber = 1
        geometricFieldUserNumber = 1
        tau=0.1
        kappa=0.05
        lambdaFieldUserNumber = 12
        fittingEquationsSetUserNumber = 13
        fittingEquationsSetFieldUserNumber = 14
        fittingDependentFieldUserNumber = 15
        fittingIndependentFieldUserNumber = 16
        fittingMaterialsFieldUserNumber = 17
        fittingProblemUserNumber = 18
        
        # Get the number of computational nodes and this computational node number
        numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
        
        # Create a 3D rectangular cartesian coordinate system
        coordinateSystem = iron.CoordinateSystem()
        coordinateSystem.CreateStart(coordinateSystemUserNumber)
        # Set the number of dimensions to 3
        coordinateSystem.DimensionSet(3)
        # Finish the creation of the coordinate system
        coordinateSystem.CreateFinish()
        
        # Create a region and assign the coordinate system to the region
        region = iron.Region()
        region.CreateStart(regionUserNumber,iron.WorldRegion)
        region.LabelSet("HeartTubeRegion")
        # Set the regions coordinate system to the 3D RC coordinate system that we have created
        region.coordinateSystem = coordinateSystem
        # Finish the creation of the region
        region.CreateFinish()
        self.region = region
        # Define basis
        # Start the creation of a tricubic Hermite basis function
        tricubicHermiteBasis = iron.Basis()
        tricubicHermiteBasis.CreateStart(tricubicHermiteBasisUserNumber)
        tricubicHermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        tricubicHermiteBasis.numberOfXi = 3
        tricubicHermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
        tricubicHermiteBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
        tricubicHermiteBasis.CreateFinish()
        
        # Start the creation of a manually generated mesh in the region
        numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements
        
        
        # Define nodes for the mesh
        nodes = iron.Nodes()
        nodes.CreateStart(region,numberOfNodes)
        nodes.CreateFinish()
        mesh = iron.Mesh()
        
        # Create the mesh. The mesh will have two components - 1. tricubic Hermite elements; 2. trilinear Lagrange elements
        mesh.CreateStart(meshUserNumber,region,3)
        mesh.NumberOfComponentsSet(1)
        mesh.NumberOfElementsSet(numberOfElements)
        
        tricubicHermiteElements = iron.MeshElements()
        tricubicHermiteElements.CreateStart(mesh,1,tricubicHermiteBasis)
        
        elementNumber = 0
        for wallElementIdx in range(1,numberOfWallElements+1):
            for lengthElementIdx in range(1,numberOfLengthElements+1):
                for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                    elementNumber = elementNumber + 1
                    localNode1 = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    if circumfrentialElementIdx == numberOfCircumfrentialElements:
                        localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                            (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    else:
                        localNode2 = localNode1 + 1
                    localNode3 = localNode1 + numberOfCircumfrentialNodes
                    localNode4 = localNode2 + numberOfCircumfrentialNodes
                    localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
                    localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
                    localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
                    localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
                    localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
                    tricubicHermiteElements.NodesSet(elementNumber,localNodes)
        
        tricubicHermiteElements.CreateFinish()
        
        # Finish the mesh creation
        mesh.CreateFinish() 
        
        # Create a decomposition for the mesh
        decomposition = iron.Decomposition()
        decomposition.CreateStart(decompositionUserNumber,mesh)
        # Set the decomposition to be a general decomposition with the specified number of domains
        decomposition.type = iron.DecompositionTypes.CALCULATED
        decomposition.numberOfDomains = numberOfComputationalNodes
        # Finish the decomposition
        decomposition.CreateFinish()
        
        # Create a field for the geometry
        geometricField = iron.Field()
        geometricField.CreateStart(geometricFieldUserNumber,region)
        # Set the decomposition to use
        geometricField.MeshDecompositionSet(decomposition)
        geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        # Set the field label
        geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
        # Set the domain to be used by the field components to be tricubic Hermite
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        # Set the scaling type
        geometricField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        # Finish creating the field
        geometricField.CreateFinish()
            
        self.setupGeometry(geometricField)        
        # Update the geometric field
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        
        
        lambdaField = iron.Field()
        lambdaField.CreateStart(lambdaFieldUserNumber,region)
        lambdaField.TypeSet(iron.FieldTypes.GENERAL)
        # Set the decomposition
        lambdaField.MeshDecompositionSet(decomposition)
        # Set the geometric field
        lambdaField.GeometricFieldSet(geometricField)
        lambdaField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
        # Set the field variables
        lambdaField.NumberOfVariablesSet(1)
        lambdaField.VariableTypesSet([iron.FieldVariableTypes.U])
        # Set the variable label
        lambdaField.VariableLabelSet(iron.FieldVariableTypes.U,"NodeLambda")
        # Set the components to be tricubic-hermite
        lambdaField.NumberOfComponentsSet(iron.FieldVariableTypes.U,9)
        for comp in range(1,10):
            lambdaField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,comp,1)
            # Set the interpolation types
            lambdaField.ComponentInterpolationSet(iron.FieldVariableTypes.U,comp,iron.FieldInterpolationTypes.NODE_BASED)
        
        lambdaField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        
        lambdaField.CreateFinish()
        # Initialise the lambda field
        for comp in range(1,10):
            lambdaField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,0.0)

        
        # Create Gauss point fitting equations set
        fittingEquationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                                     iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                     iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                                     iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
        fittingEquationsSetField = iron.Field()
        fittingEquationsSet = iron.EquationsSet()
        fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber,region,geometricField,
                fittingEquationsSetSpecification,fittingEquationsSetFieldUserNumber,fittingEquationsSetField)
        fittingEquationsSet.CreateFinish()
        
        # Create the fitting dependent field
        fittingDependentField = iron.Field()
        fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber,fittingDependentField)
        fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"FittingU")
        fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"FittingDelUdelN")
        # Set the number of components to 9
        fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,9)
        fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,9)
        # Set the field variables to be tricubic hermite
        for comp in range(1,10):
            fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,comp,1)
            fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,comp,1)
        
        # Finish creating the fitting dependent field
        fittingEquationsSet.DependentCreateFinish()
        
        # Create the fitting independent field
        fittingIndependentField = iron.Field()
        fittingEquationsSet.IndependentCreateStart(fittingIndependentFieldUserNumber,fittingIndependentField)
        fittingIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,"GaussLambda")
        fittingIndependentField.VariableLabelSet(iron.FieldVariableTypes.V,"LambdaWeight")
        # Set the number of components to 9
        fittingIndependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,9)
        fittingIndependentField.NumberOfComponentsSet(iron.FieldVariableTypes.V,9)
        # Finish creating the fitting independent field
        fittingEquationsSet.IndependentCreateFinish()
        # Initialise data point vector field to 0.0
        for comp in range(1,10):
            fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,0.0)
        # Initialise data point weight field to 1.0
            fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,comp,1.0)
        
        # Create material field (Sobolev parameters)
        fittingMaterialField = iron.Field()
        fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber,fittingMaterialField)
        fittingMaterialField.VariableLabelSet(iron.FieldVariableTypes.U,"SmoothingParameters")
        fittingEquationsSet.MaterialsCreateFinish()
        # Set kappa and tau - Sobolev smoothing parameters
        fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
        fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)
        
        # Create the fitting equations
        fittingEquations = iron.Equations()
        fittingEquationsSet.EquationsCreateStart(fittingEquations)
        # Set the fitting equations sparsity type
        fittingEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        # Set the fitting equations output type to none
        fittingEquations.outputType = iron.EquationsOutputTypes.NONE
        # Finish creating the fitting equations
        fittingEquationsSet.EquationsCreateFinish()
        
        # Create fitting problem
        fittingProblemSpecification = [iron.ProblemClasses.FITTING,
                                iron.ProblemTypes.DATA_FITTING,
                                iron.ProblemSubtypes.STATIC_FITTING]
        fittingProblem = iron.Problem()
        fittingProblem.CreateStart(fittingProblemUserNumber,fittingProblemSpecification)
        fittingProblem.CreateFinish()
        
        # Create control loops
        fittingProblem.ControlLoopCreateStart()
        fittingProblem.ControlLoopCreateFinish()
        
        # Create problem solver
        fittingSolver = iron.Solver()
        fittingProblem.SolversCreateStart()
        fittingProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,fittingSolver)
        fittingSolver.outputType = iron.SolverOutputTypes.NONE
        fittingProblem.SolversCreateFinish()

        # Create fitting solver equations and add fitting equations set to solver equations
        fittingSolverEquations = iron.SolverEquations()
        fittingProblem.SolverEquationsCreateStart()
        # Get the solver equations
        fittingSolver.SolverEquationsGet(fittingSolverEquations)
        fittingSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        fittingEquationsSetIndex = fittingSolverEquations.EquationsSetAdd(fittingEquationsSet)
        fittingProblem.SolverEquationsCreateFinish()
        
        # Prescribe boundary conditions for the fitting problem
        
        fittingBoundaryConditions = iron.BoundaryConditions()
        fittingSolverEquations.BoundaryConditionsCreateStart(fittingBoundaryConditions)
        fittingBoundaryConditions.AddNode(fittingDependentField,iron.FieldVariableTypes.U,
                                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,8,
                                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        
        fittingSolverEquations.BoundaryConditionsCreateFinish()
        
        self.fittingIndependentField = fittingIndependentField
        self.fittingDependentField = fittingDependentField
        self.lambdaField = lambdaField
        self.fittingProblem = fittingProblem            

    def solveGrowthRates(self,growthRates):
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfWallElements = self.wallElements
        numberOfGaussXi = 3 
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements
        fittingIndependentField = self.fittingIndependentField
        #The fitting dependent field will follow the 3x3 matrix order rather than iron's expectations
        xyn = numberOfCircumferentialElements*numberOfLengthElements
        for elementNumber in range (1, numberOfElements+1):
            eid = elementNumber -1
            XV = (eid%numberOfCircumferentialElements)*numberOfGaussXi
            YV = (int(eid/numberOfCircumferentialElements)%(numberOfLengthElements))*numberOfGaussXi
            ZV = int(eid/xyn)*numberOfGaussXi                        
            for xiIdx3 in range(numberOfGaussXi):
                for xiIdx2 in range(numberOfGaussXi):
                    for xiIdx1 in range(numberOfGaussXi):
                        gaussPointNumber = xiIdx1 + (xiIdx2)*numberOfGaussXi + (xiIdx3)*numberOfGaussXi*numberOfGaussXi + 1
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,0,0])
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,0,1])
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,0,2])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,4,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,1,0])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,5,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,1,1])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,6,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,1,2])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,7,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,2,0])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,8,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,2,1])                        
                        fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,9,growthRates[XV+xiIdx1,YV+xiIdx2,ZV+xiIdx3,2,2])
                        
        fittingIndependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        fittingIndependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)    
   

        try:
            self.fittingProblem.Solve()
        except Exception as e:
            print("Error during fitting solve.\n")
            raise e
        #Copy fitted dependent field to lambda field
        for comp in range(1,10):
            iron.Field.ParametersToFieldParametersComponentCopy(
                self.fittingDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                self.lambdaField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)

        #Get node values
        numberOfCircumfrentialElements = self.circumferentialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumfrentialNodes = numberOfCircumfrentialElements
        numberOfWallNodes = self.wallElements+1
        
        fieldValues = []
        x=np.zeros((8,9))
        # Create the geometric field
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                        (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes

                    for comp in range(1,10):
                        x[0,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,comp)

                        x[1,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,comp)

                        x[2,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,comp)

                        x[3,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,comp)

                        x[4,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,comp)

                        x[5,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,nodeNumber,comp)

                        x[6,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,nodeNumber,comp)

                        x[7,comp-1] = self.lambdaField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,nodeNumber,comp)
                    fieldValues.append(x.flatten().tolist())
                    
        return np.array(fieldValues)
        
    def saveResults(self,filename):
        fields = iron.Fields()
        fields.CreateRegion(self.region)
        fields.NodesExport(filename,"FORTRAN")
        fields.ElementsExport(filename,"FORTRAN")            
        fields.Finalise()
        
    def computeGradient(self,filename,numGauss=3):
        ctx = Context('Grad')
        region = ctx.getDefaultRegion()
        sir = region.createStreaminformationRegion()
        sir.createStreamresourceFile('%s.part0.exnode'%filename)
        sir.createStreamresourceFile('%s.part0.exelem'%filename)
        region.read(sir)
        
        fieldModule = region.getFieldmodule()        
        coordinatesField = fieldModule.findFieldByName('Geometry').castFiniteElement()
        ratesField = fieldModule.findFieldByName('NodeLambda')
        gradient = fieldModule.createFieldGradient(ratesField,coordinatesField)
        qgen = GaussQuadarture()
        self.gaussPointLocations = qgen.getQuadratureFeaturesFor(numGauss)[0]
        self.numGaussPoints = self.gaussPointLocations.shape[0]**3        
        self.numPointsPerGauss = 1 # (xi + number of neighbors(6) + self)
        #Create the derivative field
        fieldCache = fieldModule.createFieldcache()
        nx = self.circumferentialElements*self.gaussPointLocations.shape[0]
        ny = self.axialElements*self.gaussPointLocations.shape[0]
        nz = self.wallElements*self.gaussPointLocations.shape[0]
        dgrad = np.zeros((nx,ny,nz,3,3,3))
        ngauss = self.gaussPointLocations.shape[0]
        xyn = self.circumferentialElements*self.axialElements
        elements = dict()
        mesh = fieldModule.findMeshByDimension(3)
        ei   = mesh.createElementiterator()
        elem = ei.next()
        while elem.isValid():
            elements[elem.getIdentifier()-1] = elem
            elem = ei.next()
        
        for eid in elements:
            elem = elements[eid]
            XV = (eid%self.circumferentialElements)*ngauss
            YV = (int(eid/self.circumferentialElements)%(self.axialElements))*ngauss
            ZV = int(eid/xyn)*ngauss            
            
            for k,xi3 in enumerate(self.gaussPointLocations):
                for j,xi2 in enumerate(self.gaussPointLocations):
                    for i,xi1 in enumerate(self.gaussPointLocations):
                        fieldCache.setMeshLocation(elem,[xi1,xi2,xi3])
                        _,vol = gradient.evaluateReal(fieldCache,27)                       
                        dgrad[XV+i,YV+j,ZV+k] = np.array(vol).reshape((3,3,3))
        
        np.save('%srates'%filename,dgrad)
            
if __name__ == '__main__':
    obj = Fitting(4,4,1)
    obj.setupProblem()
    rates = np.zeros((12,12,3,3,3))
    rates[:,:,:,0,0] = 1.0
    rates[:,:,:,1,1] = 2.0
    rates[:,:,:,2,2] = 3.0
    rts = obj.solveGrowthRates(rates)
    obj.saveResults('ftest')
    obj.computeGradient('ftest')