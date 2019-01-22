"""
Generates a 3-D unit tube mesh with variable numbers of elements around, along and
through wall, plus variable wall thickness for unit diameter.
A reference coordinate field is created 
"""

import numpy as np
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.context import Context

class Tube3d(object):
    '''
    Creates a tube with specified number of wall, circumferential and lengthwise elements
    In addition to coordinates field a referenceCoordinates field is defined
    '''
    def __init__(self, xi1Elements=4,xi2Elements=1,xi3Elements=1,timepoints=1):
        '''
        Constructor
        '''
        self.context = Context('Tube')    
        self.lengthElements = xi2Elements
        self.circumferentialElements = xi1Elements
        self.wallElements = xi3Elements
        self.numberOfTimePoints = timepoints
        
    def generateMesh(self,region):
        lengthElements = self.lengthElements
        circumferentialElements = self.circumferentialElements
        wallElements = self.wallElements
        wallThickness = 1

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = fm.createFieldFiniteElement(3)
        coordinates.setName('coordinates')
        coordinates.setManaged(True)
        coordinates.setTypeCoordinate(True)
        coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        coordinates.setComponentName(1, 'x')
        coordinates.setComponentName(2, 'y')
        coordinates.setComponentName(3, 'z')

        self.coordinates = coordinates

        referenceCoordinates = fm.createFieldFiniteElement(3)
        referenceCoordinates.setName('referenceCoordinates')
        referenceCoordinates.setManaged(True)
        referenceCoordinates.setTypeCoordinate(True)
        referenceCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        referenceCoordinates.setComponentName(1, 'x')
        referenceCoordinates.setComponentName(2, 'y')
        referenceCoordinates.setComponentName(3, 'z')

        self.referenceCoordinates = referenceCoordinates
                
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        self.tvals = np.array([0.0])
        if self.numberOfTimePoints > 1:
            self.tvals = np.arange(self.numberOfTimePoints,dtype=np.float)
            timeSequence = fm.getMatchingTimesequence(self.tvals.tolist())        
            nodetemplate.setTimesequence(coordinates,timeSequence)
            
        nodetemplate.defineField(referenceCoordinates)
        for field in [coordinates,referenceCoordinates]:
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(field, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(tricubicHermiteBasis)

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)
        result = elementtemplate.defineField(referenceCoordinates, -1, eft)
        
        cache = fm.createFieldcache()
        self.nodes = dict()
        # create nodes
        radiansPerElementAround = 2.0*np.pi/circumferentialElements
        wallThicknessPerElement = wallThickness/wallElements
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / lengthElements ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        numberOfWallNodes = wallElements + 1
        numberOfCircumfrentialNodes = circumferentialElements
        numberOfLengthNodes = lengthElements + 1        
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            radius = 0.5 + wallThickness*((wallNodeIdx-1)/(numberOfWallNodes - 1.0))
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                x[2] = float(lengthNodeIdx-1)/ lengthElements
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
                    radiansAround = circumfrentialNodeIdx*radiansPerElementAround
                    cosRadiansAround = np.cos(radiansAround)
                    sinRadiansAround = np.sin(radiansAround)
                    x[0] = radius*cosRadiansAround
                    x[1] = radius*sinRadiansAround
                    dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                    dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                    dx_ds3[0] = wallThicknessPerElement*cosRadiansAround
                    dx_ds3[1] = wallThicknessPerElement*sinRadiansAround
                    node = nodes.createNode(nodeNumber, nodetemplate)
                    self.nodes[nodeNumber] = node
                    cache.setNode(node)
                    
                    for t in self.tvals:
                        cache.setTime(t)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    referenceCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)                        

        # create elements
        elementNumber = 0
        elems = []
        for wallElementIdx in range(1,wallElements+1):
            for lengthElementIdx in range(1,lengthElements+1):
                for circumfrentialElementIdx in range(1,circumferentialElements+1):
                    elementNumber = elementNumber + 1
                    localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*circumferentialElements + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                    if circumfrentialElementIdx == circumferentialElements:
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
                    element = mesh.createElement(elementNumber, elementtemplate)
                    result = element.setNodesByIdentifier(eft, localNodes)
                    elems.append(element)
        
        fm.defineAllFaces()
        fm.endChange()

    def loadStateFromFile(self,filename,time=0,coordinateFieldName='coordinates',nodeOffset=272):
        '''
        Load a ex file and gets the node's geometry data from the coordinateFieldName field
        if nodeMap is supplied then the mesh nodes values are assigned based on this mapping, else by node number
        nodeOffset - some exfile nodes will not start at 1, if so specify offset such that the nodeMapping can be done
        '''
        context = Context('Load')
        region = context.getDefaultRegion()
        region.readFile(filename)
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        fm.beginChange()
        coordinateValues = dict()
        coordinates = fm.findFieldByName(coordinateFieldName).castFiniteElement()
        if not coordinates.isValid():
            raise ValueError('Unable to find coordinate field %s, failed to load coordinates from file %s' % (coordinateFieldName,filename))

        nodeset = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        niter = nodeset.createNodeiterator()
        node = niter.next()
        while node.isValid():
            cache.setNode(node)
            _,coord = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,3)
            _,dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            _,dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            _,dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            _,dx_ds1ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
            _,dx_ds1ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
            _,dx_ds2ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)
            _,dx_ds1ds2ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)
                
            coordinateValues[node.getIdentifier()-nodeOffset] = [coord,dx_ds1,dx_ds2,dx_ds3,dx_ds1ds2,dx_ds1ds3,dx_ds2ds3,dx_ds1ds2ds3]
            node = niter.next()
            
        fm.endChange()
        
        #Load it into current mesh
        region = self.context.getDefaultRegion()
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        fm.beginChange()
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()
        cache.setTime(time)
        for nd,node in self.nodes.items():
            cache.setNode(node)
            coord,dx_ds1,dx_ds2,dx_ds3,dx_ds1ds2,dx_ds1ds3,dx_ds2ds3,dx_ds1ds2ds3 = coordinateValues[nd]
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, coord)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, dx_ds1ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, dx_ds1ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, dx_ds2ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, dx_ds1ds2ds3)
        
        smoothing = fm.createFieldsmoothing()
        coordinates.smooth(smoothing)                            
        fm.endChange()

    def loadReferenceCoordinatesFromFile(self,filename,coordinateFieldName='coordinates',nodeOffset=272):
        '''
        Load a ex file and gets the node's geometry data from the coordinateFieldName field
        if nodeMap is supplied then the mesh nodes values are assigned based on this mapping, else by node number
        nodeOffset - some exfile nodes will not start at 1, if so specify offset such that the nodeMapping can be done
        '''
        context = Context('Load')
        region = context.getDefaultRegion()
        region.readFile(filename)
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        fm.beginChange()
        coordinateValues = dict()
        coordinates = fm.findFieldByName(coordinateFieldName).castFiniteElement()
        if not coordinates.isValid():
            raise ValueError('Unable to find coordinate field %s, failed to load coordinates from file %s' % (coordinateFieldName,filename))

        nodeset = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        niter = nodeset.createNodeiterator()
        node = niter.next()
        while node.isValid():
            cache.setNode(node)
            _,coord = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,3)
            _,dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            _,dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            _,dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            _,dx_ds1ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
            _,dx_ds1ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
            _,dx_ds2ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)
            _,dx_ds1ds2ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)
                
            coordinateValues[node.getIdentifier()-nodeOffset] = [coord,dx_ds1,dx_ds2,dx_ds3,dx_ds1ds2,dx_ds1ds3,dx_ds2ds3,dx_ds1ds2ds3]
            node = niter.next()
            
        fm.endChange()
        
        #Load it into current mesh
        region = self.context.getDefaultRegion()
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        fm.beginChange()
        coordinates = self.referenceCoordinates
        for nd,node in self.nodes.items():
            cache.setNode(node)
            coord,dx_ds1,dx_ds2,dx_ds3,dx_ds1ds2,dx_ds1ds3,dx_ds2ds3,dx_ds1ds2ds3 = coordinateValues[nd]
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, coord)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, dx_ds1ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, dx_ds1ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, dx_ds2ds3)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, dx_ds1ds2ds3)
        
        smoothing = fm.createFieldsmoothing()
        coordinates.smooth(smoothing)                            
        fm.endChange()    

    def getRegion(self):
        return self.context.getDefaultRegion()
    
    def writeToFile(self,filename):
        region =self.context.getDefaultRegion()
        region.writeFile(filename)

    def writeToFileAtTime(self,filename,t):
        region =self.context.getDefaultRegion()
        sir = region.createStreaminformationRegion()
        tfile = sir.createStreamresourceFile(filename)
        sir.setResourceAttributeReal(tfile,sir.ATTRIBUTE_TIME,t)
        region.write(sir)
        
    def saveAsNumpy(self,filename,t):
        coordinates = np.zeros((len(self.nodes),24))
        region = self.context.getDefaultRegion()
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        cache.setTime(t)        
        for nd,node in self.nodes.items():
            cache.setNode(node)
            _,coord = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,3)
            _,dx_ds1 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            _,dx_ds2 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            _,dx_ds3 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            _,dx_ds1ds2 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
            _,dx_ds1ds3 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
            _,dx_ds2ds3 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)
            _,dx_ds1ds2ds3 = self.coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)
            coordinates[nd-1] = np.concatenate((coord,dx_ds1,dx_ds2,dx_ds3,dx_ds1ds2,dx_ds1ds3,dx_ds2ds3,dx_ds1ds2ds3))
        
        np.save(filename,coordinates)

if __name__ == '__main__':
    obj = Tube3d(16,16,1)
    obj.generateMesh(obj.getRegion())
    obj.loadStateFromFile(r'../volumecalculator/stage10.exregion')
    obj.writeToFile('tubeTest.ex2')