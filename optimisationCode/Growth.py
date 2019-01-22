from __future__ import print_function
import numpy as np
from opencmiss.iron import iron
class TubeGrower(object):
    maxSolvebleRate = 0.1
    def __init__(self, circumferentialElements,axialElements,wallElements,discret=10,length=6,innerRadius=1.0,outerRadius=2,fixBottom=True,fixTop=False,DMBC=True,humphrey=False,neoHookean=False):
        self.circumferentialElements = circumferentialElements
        self.axialElements = axialElements
        self.wallElements = wallElements
        self.length=length
        self.innerRadius=innerRadius
        self.outerRadius=outerRadius
        self.fixBottom=fixBottom
        self.fixTop=fixTop
        self.DMBC = DMBC
        self.humphrey = humphrey
        self.neoHookean = neoHookean
        xi = np.linspace(0,1.0,discret)
        xi1,xi2,xi3 = np.meshgrid(xi,xi,xi)
        self.xi = np.c_[xi1.flatten(),xi2.flatten(),xi3.flatten()].T.tolist()
        self.numPoints = discret**3
                
    def createBoundaryConditions(self,nonlinearEquations):
        if (self.DMBC):
            dependentField = self.dependentField
            boundaryConditions = iron.BoundaryConditions()
            nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)
            numberOfCircumferentialElementsPerQuarter = int(self.circumferentialElements/4)
            numberOfCircumferentialElements = self.circumferentialElements
            numberOfLengthNodes = self.axialElements+1
            numberOfCircumferentialNodes = numberOfCircumferentialElements       
            #nodelistx = [91,99,107,115,123,131,139]
            #nodelisty = [73,77,81,85,139]
            #nodelistz = [65,66,67,68,69,70,71,72,137,138,139,140,141,142,143,144]
            nodelistx = [20,28] 
            nodelisty = [20,29,31]
            nodelistz = [20,18]
            for nodeNumber in nodelistz:        
            # Fix S3 and Z direction
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0) 
            for nodeNumber in nodelisty:
            # Fix S2 and Y direction
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0) 
            for nodeNumber in nodelistx:
            # Fix S1 and X direction
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)       
            # Changing the X, Y, Z of the points... 
			# nodeNumbers show where needs to be fixed 
			
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,24,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,24,2,iron.BoundaryConditionsTypes.FIXED,0.35)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,24,3,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,24,2,iron.BoundaryConditionsTypes.FIXED,-0.2)

            nonlinearEquations.BoundaryConditionsCreateFinish()
        else:
            dependentField = self.dependentField
            boundaryConditions = iron.BoundaryConditions()
            nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)
            numberOfCircumfrentialElementsPerQuarter = int(self.circumferentialElements/4)
            numberOfCircumfrentialElements = self.circumferentialElements
            numberOfLengthNodes = self.axialElements+1
            numberOfCircumfrentialNodes = numberOfCircumfrentialElements
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                if (lengthNodeIdx == 1 and self.fixBottom) or (lengthNodeIdx == numberOfLengthNodes and self.fixTop):
                    for wallNodeIdx in range(1,self.wallElements+1):
                        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialElements+1):
                            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
                            # Fix z direction
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                            # Fix S1 (circumfrential) direction derivatives
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                            # Fix S2 (length) direction derivatives
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                            # Fix S3 (wall) direction derivatives
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                    #Set symmetry conditions on the ring to prevent rotation                                      
                    nodeNumber = 1 + (lengthNodeIdx-1)*numberOfCircumfrentialNodes 
                    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                    nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
                    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                    nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
                    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                    nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
                    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            nonlinearEquations.BoundaryConditionsCreateFinish()
		
    def setupGrowthRates(self,growthElementRate):                     #12x3x3x3x6
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfWallElements = self.wallElements
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements
        growthCellMLParametersField = self.growthCellMLParametersField
        for elementNumber in range (1, numberOfElements+1):
            for xi1Index in range (1,3+1):
                for xi2Index in range (1,3+1):
                    for xi3Index in range (1,3+1):
                        gaussPointNumber = (xi3Index - 1)*3*3 + (xi2Index - 1)*3 + xi1Index                   # 1-27
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,0])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,1])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,2])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,4,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,3])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,5,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,4])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,6,growthElementRate[elementNumber-1,xi1Index-1,xi2Index-1,xi3Index-1,5])
        growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES) 
				
    def setupGeometry(self,geometricField):
        numberOfCircumfrentialElements = self.circumferentialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumfrentialNodes = numberOfCircumfrentialElements
        numberOfWallNodes = self.wallElements+1
        # Create the geometric field
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    radius = self.innerRadius + (self.outerRadius - self.innerRadius)*float(wallNodeIdx-1)/float(numberOfWallNodes)
                    theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*np.pi
                    x = radius*np.cos(theta)
                    y = radius*np.sin(theta)
                    xtangent = -np.sin(theta)
                    ytangent = np.cos(theta)
                    xnormal = np.cos(theta)
                    ynormal = np.sin(theta)
                    z = float(lengthNodeIdx-1)/float(numberOfLengthNodes)*self.length
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)
            
    def setupProblem(self,showProgress=True):

        pInit = -8.0 
        fibreAngle = 0.0 
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
        trilinearLagrangeBasisUserNumber = 2
        meshUserNumber = 1
        decompositionUserNumber = 1
        geometricFieldUserNumber = 1
        originalGeometricFieldUserNumber = 20
        fibreFieldUserNumber = 2
        dependentFieldUserNumber = 3
        equationsSetUserNumber = 1
        equationsSetFieldUserNumber = 5
        growthCellMLUserNumber = 1
        growthCellMLModelsFieldUserNumber = 6
        growthCellMLStateFieldUserNumber = 7
        growthCellMLParametersFieldUserNumber = 8
        constitutiveCellMLUserNumber = 2
        constitutiveCellMLModelsFieldUserNumber = 9
        constitutiveCellMLParametersFieldUserNumber = 10
        constitutiveCellMLIntermediateFieldUserNumber = 11
        problemUserNumber = 1
        
        numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()

        coordinateSystem = iron.CoordinateSystem()
        coordinateSystem.CreateStart(coordinateSystemUserNumber)
        coordinateSystem.DimensionSet(3)
        coordinateSystem.CreateFinish()
        region = iron.Region()
        region.CreateStart(regionUserNumber,iron.WorldRegion)
        region.LabelSet("HeartTubeRegion")
        region.coordinateSystem = coordinateSystem
        region.CreateFinish()
        
        tricubicHermiteBasis = iron.Basis()
        tricubicHermiteBasis.CreateStart(tricubicHermiteBasisUserNumber)
        tricubicHermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        tricubicHermiteBasis.numberOfXi = 3
        tricubicHermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
        tricubicHermiteBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
        tricubicHermiteBasis.CreateFinish()
        trilinearLagrangeBasis = iron.Basis()
        trilinearLagrangeBasis.CreateStart(trilinearLagrangeBasisUserNumber)
        trilinearLagrangeBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        trilinearLagrangeBasis.numberOfXi = 3
        trilinearLagrangeBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
        trilinearLagrangeBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
        trilinearLagrangeBasis.CreateFinish()
        
        numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements
        nodes = iron.Nodes()
        nodes.CreateStart(region,numberOfNodes)
        nodes.CreateFinish()
        mesh = iron.Mesh()
        mesh.CreateStart(meshUserNumber,region,3)
        mesh.NumberOfComponentsSet(2)
        mesh.NumberOfElementsSet(numberOfElements)
        tricubicHermiteElements = iron.MeshElements()
        tricubicHermiteElements.CreateStart(mesh,1,tricubicHermiteBasis)
        trilinearLagrangeElements = iron.MeshElements()
        trilinearLagrangeElements.CreateStart(mesh,2,trilinearLagrangeBasis)
        
        elementNumber = 0
        for wallElementIdx in range(1,numberOfWallElements+1):
            for lengthElementIdx in range(1,numberOfLengthElements+1):
                for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                    elementNumber = elementNumber + 1
                    localNode1 = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialNodes + (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    if circumfrentialElementIdx == numberOfCircumfrentialElements:
                        localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
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
                    trilinearLagrangeElements.NodesSet(elementNumber,localNodes)
        tricubicHermiteElements.CreateFinish()
        trilinearLagrangeElements.CreateFinish()
        mesh.CreateFinish() 
        
        decomposition = iron.Decomposition()
        decomposition.CreateStart(decompositionUserNumber,mesh)
        decomposition.type = iron.DecompositionTypes.CALCULATED
        decomposition.numberOfDomains = numberOfComputationalNodes
        decomposition.CreateFinish()
        
        geometricField = iron.Field()
        geometricField.CreateStart(geometricFieldUserNumber,region)
        geometricField.MeshDecompositionSet(decomposition)
        geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        geometricField.CreateFinish()

        originalGeometricField= iron.Field()
        originalGeometricField.CreateStart(originalGeometricFieldUserNumber,region)
        originalGeometricField.MeshDecompositionSet(decomposition)
        originalGeometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        originalGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,"OriginalGeometry")
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        originalGeometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        originalGeometricField.CreateFinish()
        
        self.setupGeometry(geometricField)        
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        iron.Field.ParametersToFieldParametersComponentCopy(geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
                                                            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
        iron.Field.ParametersToFieldParametersComponentCopy(geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
                                                            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
        iron.Field.ParametersToFieldParametersComponentCopy(geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
                                                            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
        originalGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        originalGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        
        fibreField = iron.Field()
        fibreField.CreateStart(fibreFieldUserNumber,region)
        fibreField.TypeSet(iron.FieldTypes.FIBRE)
        fibreField.MeshDecompositionSet(decomposition)
        fibreField.GeometricFieldSet(geometricField)
        fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
        fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U,3)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)
        fibreField.CreateFinish()
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,fibreAngle)
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,0.0)
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
        fibreField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        fibreField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        dependentField = iron.Field()
        dependentField.CreateStart(dependentFieldUserNumber,region)
        dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
        dependentField.MeshDecompositionSet(decomposition)
        dependentField.GeometricFieldSet(geometricField) 
        dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
        dependentField.NumberOfVariablesSet(5)
        dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2,iron.FieldVariableTypes.U3])

        dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"del U/del n")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"Strain")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"Stress")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U3,"Growth")
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U3,6)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
        for comp in range(1,7):
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,comp,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,comp,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,comp,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
        dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        dependentField.CreateFinish()
        for comp in range(1,4):
            iron.Field.ParametersToFieldParametersComponentCopy(geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                                                                dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
        iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,pInit)
        
        equationsSetField = iron.Field()
        equationsSet = iron.EquationsSet()
        equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,iron.EquationsSetTypes.FINITE_ELASTICITY,iron.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
        equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
        equationsSet.CreateFinish()
        
        equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
        equationsSet.DependentCreateFinish()
        
        equations = iron.Equations()
        equationsSet.EquationsCreateStart(equations)
        equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        equations.outputType = iron.EquationsOutputTypes.NONE
        equationsSet.EquationsCreateFinish()
        
        growthCellML = iron.CellML()
        growthCellML.CreateStart(growthCellMLUserNumber,region)
        
        growthCellMLIdx = growthCellML.ModelImport("simplefullgrowth.cellml")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibresheetrate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrenormalrate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetnormalrate")        
        growthCellML.CreateFinish()
        
        growthCellML.FieldMapsCreateStart()
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,1,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,2,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda12",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,4,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda13",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,5,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda23",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U3,6,iron.FieldParameterSetTypes.VALUES)
        growthCellML.FieldMapsCreateFinish()
        
        # Create the CELL models field
        growthCellMLModelsField = iron.Field()
        growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
        growthCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthModelMap")
        growthCellML.ModelsFieldCreateFinish()
        
        # Create the CELL parameters field
        growthCellMLParametersField = iron.Field()
        growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
        growthCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthParameters")
        growthCellML.ParametersFieldCreateFinish()
        
        # Create the CELL state field
        growthCellMLStateField = iron.Field()
        growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
        growthCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthState")
        growthCellML.StateFieldCreateFinish()
        
        # Create the CellML environment for the consitutative law
        constitutiveCellML = iron.CellML()
        constitutiveCellML.CreateStart(constitutiveCellMLUserNumber,region)
        if (self.humphrey):
            constitutiveCellMLIdx = constitutiveCellML.ModelImport("Humphrey.cellml")
        elif (self.neoHookean):
            constitutiveCellMLIdx = constitutiveCellML.ModelImport("neoHookean.cellml")
        else: 
            constitutiveCellMLIdx = constitutiveCellML.ModelImport("mooneyrivlin.cellml")
        # Flag the CellML variables that OpenCMISS will supply
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C11")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C12")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C13")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C22")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C23")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C33")
        #constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c1")
        #constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c2")
        # Flag the CellML variables that OpenCMISS will obtain
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev11")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev12")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev13")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev22")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev23")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev33")
        constitutiveCellML.CreateFinish()
        
        # Create CellML <--> OpenCMISS field maps
        constitutiveCellML.FieldMapsCreateStart()
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C11",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C12",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C13",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C22",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C23",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6,iron.FieldParameterSetTypes.VALUES,constitutiveCellMLIdx,"equations/C33",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev11",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev12",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev13",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev22",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev23",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev33",iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.FieldMapsCreateFinish()
        
        # Create the CELL models field
        constitutiveCellMLModelsField = iron.Field()
        constitutiveCellML.ModelsFieldCreateStart(constitutiveCellMLModelsFieldUserNumber,constitutiveCellMLModelsField)
        constitutiveCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveModelMap")
        constitutiveCellML.ModelsFieldCreateFinish()
        
        # Create the CELL parameters field
        constitutiveCellMLParametersField = iron.Field()
        constitutiveCellML.ParametersFieldCreateStart(constitutiveCellMLParametersFieldUserNumber,constitutiveCellMLParametersField)
        constitutiveCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveParameters")
        constitutiveCellML.ParametersFieldCreateFinish()
        
        # Create the CELL intermediate field
        constitutiveCellMLIntermediateField = iron.Field()
        constitutiveCellML.IntermediateFieldCreateStart(constitutiveCellMLIntermediateFieldUserNumber,constitutiveCellMLIntermediateField)
        constitutiveCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveIntermediate")
        constitutiveCellML.IntermediateFieldCreateFinish()
        
        # Define the problem
        problem = iron.Problem()
        problemSpecification = [iron.ProblemClasses.ELASTICITY,iron.ProblemTypes.FINITE_ELASTICITY,iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
        problem.CreateStart(problemUserNumber,problemSpecification)
        problem.CreateFinish()
        
        # Create control loops
        timeLoop = iron.ControlLoop()
        problem.ControlLoopCreateStart()
        problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],timeLoop)
        problem.ControlLoopCreateFinish()
        
        # Create problem solvers
        odeIntegrationSolver = iron.Solver()
        nonlinearSolver = iron.Solver()
        linearSolver = iron.Solver()
        cellMLEvaluationSolver = iron.Solver()
        problem.SolversCreateStart()
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
        nonlinearSolver.outputType = iron.SolverOutputTypes.NONE
        #if showProgress:
            # nonlinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
        #    nonlinearSolver.outputType = iron.SolverOutputTypes.MONITOR
        nonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
        nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
        nonlinearSolver.NewtonLinearSolverGet(linearSolver)
        linearSolver.linearType = iron.LinearSolverTypes.DIRECT
        problem.SolversCreateFinish()
        
        # Create nonlinear equations and add equations set to solver equations
        nonlinearEquations = iron.SolverEquations()
        problem.SolverEquationsCreateStart()
        nonlinearSolver.SolverEquationsGet(nonlinearEquations)
        nonlinearEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
        problem.SolverEquationsCreateFinish()
        
        # Create CellML equations and add growth and constitutive equations to the solvers
        growthEquations = iron.CellMLEquations()
        constitutiveEquations = iron.CellMLEquations()
        problem.CellMLEquationsCreateStart()
        odeIntegrationSolver.CellMLEquationsGet(growthEquations)
        growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
        cellMLEvaluationSolver.CellMLEquationsGet(constitutiveEquations)
        constitutiveEquationsIndex = constitutiveEquations.CellMLAdd(constitutiveCellML)
        problem.CellMLEquationsCreateFinish()
        
        # Prescribe boundary conditions (absolute nodal parameters)
        self.dependentField = dependentField
        self.elementNumber = elementNumber
        self.growthCellMLParametersField = growthCellMLParametersField
        self.constitutiveCellMLIntermediateField = constitutiveCellMLIntermediateField
        self.constitutiveCellMLParametersField = constitutiveCellMLParametersField
        self.geometricField = geometricField
        self.originalGeometricField = originalGeometricField
        self.timeLoop = timeLoop
        self.problem = problem
        self.growthCellMLStateField = growthCellMLStateField
        self.createBoundaryConditions(nonlinearEquations)
        self.region = region
        self.equationsSet = equationsSet
        
    def resetFieldsForNewSolution(self):
        #Reset to initial geometry 
        for comp in range(1,4):
            iron.Field.ParametersToFieldParametersComponentCopy(
                self.originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
            iron.Field.ParametersToFieldParametersComponentCopy(
                self.originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)        
        self.geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)        
        #Reset the dependent field
        for comp in range(1,7):
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.U1,iron.FieldParameterSetTypes.VALUES,comp,0.0)            
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.U2,iron.FieldParameterSetTypes.VALUES,comp,0.0)
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.U3,iron.FieldParameterSetTypes.VALUES,comp,0.0)
        for comp in range(1,4):
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,comp,0.0)
        iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,4,0.0)
        for comp in range(1,7):
            self.constitutiveCellMLIntermediateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,0.0)
            self.constitutiveCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,1.0)
        self.constitutiveCellMLIntermediateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLIntermediateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)                   
        #Reset growth field
        for comp in range(1,4):
            self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,1.0)
        for comp in range(4,7):
            self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,0.0)            
        self.growthCellMLStateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.growthCellMLStateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

    def solveAndGetSurfaceDescriptors(self,growthRate):
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumferentialNodes = numberOfCircumferentialElements
        numberOfWallNodes = self.wallElements+1
        numberOfWallElements = self.wallElements
        numberOfNodes = numberOfCircumferentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements		
        growthRateDivided = np.zeros((16,6))                        #16*6
        maxRate = growthRate.max()
        exponent = maxRate/self.maxSolvebleRate
        if exponent>1.0:
            divi = np.power(10,int(np.log10(exponent))+1)
            nsteps = divi
        else:
            nsteps = 1
            divi = 1.0
        #print("Using ",nsteps," timesteps(s) maxgrowth rate is ",self.maxSolvebleRate," new maxrate is ",maxRate/divi," rates ",growthRate," new rates ",growthRate/divi)
        growthRateDivided = growthRate/divi                                      #16*6
        #growthElementRate = self.growthRatesInterpolation(growthRate)
        growthElems = np.zeros((numberOfElements,3,3,3,6))                        #12*3*3*3*6
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
        # Initialise the parameters field
        self.setupGrowthRates(growthElems)
        self.resetFieldsForNewSolution()
        time = 0.0
        timeIncrement = 1.0
        for _ in range(nsteps):
            self.timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
            # Solve the problem
            self.problem.Solve()
            # Set geometric field to current deformed geometry
            for comp in range(1,4):
                iron.Field.ParametersToFieldParametersComponentCopy(self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                                                                    self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
            # Reset growth state to 1.0
            for comp in range(1,4):
                self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,1.0)
            for comp in range(4,7):
                self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,0.0)
                
            self.geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            self.geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                
            self.growthCellMLStateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            self.growthCellMLStateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        return self.getSampleCoordinates()
		
    def getVolume(self):
        vol = 0.0
        for eno in range(1,self.elementNumber+1):
            vol += self.geometricField.GeometricParametersElementVolumeGet(eno)
        return vol
        
    def getSampleCoordinates(self):        
        #Get coordinate values at discretized points
        coords = None
        for eno in range(1,self.elementNumber+1):
            gcoord = self.geometricField.ParameterSetInterpolateMultipleXiDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,eno,self.xi,(3,self.numPoints))
            if not coords is None:
                coords = np.r_[coords,np.array(gcoord).T]
            else:
                coords = np.array(gcoord).T
        rcoords = (coords*1e5).astype('int')
        #print("Volume ",self.getVolume())
        return rcoords
    
    def saveResults(self,filename):
        fields = iron.Fields()
        fields.CreateRegion(self.region)
        fields.NodesExport(filename,"FORTRAN")
        fields.ElementsExport(filename,"FORTRAN")            
        fields.Finalise()
    
