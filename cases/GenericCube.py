import sys

#Write Functions
def writeGiD_Header(f):
    f.write('MESH "model" dimension 3 ElemType Sphere Nnode 1'+'\n')
  
def writeGiD_CoordBegin(f):
    f.write('Coordinates'+'\n')
  
def writeGiD_CoordEnd(f):
    f.write('end coordinates'+'\n')
  
def writeGiD_ElemBegin(f):
    f.write('Elements'+'\n')

def writeGiD_Elems(f,s,r): 
    for i in range(1,s+1):
        f.write(str(i)+' '+str(i)+' '+str("%.4f"%r)+'\n')
  
def writeGiD_ElemEnd(f):
    f.write('end elements'+'\n')
  
def writeSize(f,s):
    f.write(str(s)+'\n')
  
def writeCoord(f,c):
    f.write(c+'\n')
    
#Utilities
def generateCluster(c_nodes,s,r):
  
    ss = int(round(pow(s,1.0/3)))
    
    center = [0,0,0]
    
    for i in range(0,ss):
        for j in range(0,ss):
            for k in range(0,ss):
              
                x = center[0] + (r-center[0])/ss * i
                y = center[1] + (r-center[1])/ss * j
                z = center[2] + (r-center[2])/ss * k
              
                c_nodes.append([x,y,z])

#Model Generators
def genericCube(xSize,ySize,zSize,cube):
  
    print 'Generating: ' + 'genericCube' + str(xSize) + 'x' + str(ySize) + 'x' + str(zSize) + '...'
    
    index = 1
    
    size = xSize * ySize * zSize
    radius = ((cube['top'][0]-cube['bot'][0])/(xSize-1.0))/2
    
    msh_filename = 'genericCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.msh'
    pts_filename = 'genericCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.pts'
    
    msh=open(msh_filename, 'w+')
    pts=open(pts_filename, 'w+')
    
    writeGiD_Header(msh)
    writeGiD_CoordBegin(msh)
    writeSize(pts,size)
    
    for i in range(0,xSize):
        for j in range(0,ySize):
            for k in range(0,zSize):
              
                x = cube['bot'][0] + ((cube['top'][0]-cube['bot'][0])/(xSize-1.0)) * i
                y = cube['bot'][1] + ((cube['top'][1]-cube['bot'][1])/(ySize-1.0)) * j
                z = cube['bot'][2] + ((cube['top'][2]-cube['bot'][2])/(zSize-1.0)) * k
                
                coords = str(index)+' '+str('%.4f'%x)+' '+str('%.4f'%y)+' '+str('%.4f'%z)
                
                writeCoord(msh,coords)
                writeCoord(pts,coords)
                
                index += 1
    
    writeGiD_CoordEnd(msh)
  
    writeGiD_ElemBegin(msh)
    writeGiD_Elems(msh,size,radius)
    writeGiD_ElemEnd(msh)

def offsetCube(xSize,ySize,zSize,cubes):
  
    print 'Generating: ' + 'offsetCube' + str(xSize) + 'x' + str(ySize) + 'x' + str(zSize) + '...'
        
    index = 1
    
    cube = cubes[0]
    
    size = xSize * ySize * zSize
    radius = ((cube['top'][0]-cube['bot'][0])/(xSize-1.0))/2
    
    msh_filename = 'offsetCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.msh'
    pts_filename = 'offsetCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.pts'
    
    msh=open(msh_filename, 'w+')
    pts=open(pts_filename, 'w+')
    
    writeGiD_Header(msh)
    writeGiD_CoordBegin(msh)
    writeSize(pts,size)
    
    for cube in cubes:
        for i in range(0,xSize):
            for j in range(0,ySize):
                for k in range(0,zSize):
                  
                    x = cube['bot'][0] + ((cube['top'][0]-cube['bot'][0])/(xSize-1.0)) * i
                    y = cube['bot'][1] + ((cube['top'][1]-cube['bot'][1])/(ySize-1.0)) * j
                    z = cube['bot'][2] + ((cube['top'][2]-cube['bot'][2])/(zSize-1.0)) * k
                    
                    coords = str(index)+' '+str('%.4f'%x)+' '+str('%.4f'%y)+' '+str('%.4f'%z)
                    
                    writeCoord(msh,coords)
                    writeCoord(pts,coords)
                    
                    index += 1
    
    writeGiD_CoordEnd(msh)

    writeGiD_ElemBegin(msh)
    writeGiD_Elems(msh,size*len(cubes),radius)
    writeGiD_ElemEnd(msh)
  
def clusterCube(xSize,ySize,zSize,cube,cluster,c_radiI):
  
    print 'Generating: ' + 'clusterCube' + str(xSize) + 'x' + str(ySize) + 'x' + str(zSize) + 'x' + str(len(cluster)) + '...'
    
    index = 1
    
    size = xSize * ySize * zSize * len(cluster)
    radius = (c_radiI / pow(len(cluster),1.0/3))/2
    
    msh_filename = 'clusterCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'x'+str(len(cluster))+'.'+str(int(round(radius*1000000)))+'.msh'
    pts_filename = 'clusterCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'X'+str(len(cluster))+'.'+str(int(round(radius*1000000)))+'.pts'
    
    msh=open(msh_filename, 'w+')
    pts=open(pts_filename, 'w+')
    
    writeGiD_Header(msh)
    writeGiD_CoordBegin(msh)
    writeSize(pts,size)
    
    for i in range(0,xSize):
        for j in range(0,ySize):
            for k in range(0,zSize):
                for n in cluster:
              
                    x = n[0] + cube['bot'][0] + ((cube['top'][0]-cube['bot'][0])/(xSize-1.0)) * i
                    y = n[1] + cube['bot'][1] + ((cube['top'][1]-cube['bot'][1])/(ySize-1.0)) * j
                    z = n[2] + cube['bot'][2] + ((cube['top'][2]-cube['bot'][2])/(zSize-1.0)) * k
                    
                    coords = str(index)+' '+str('%.4f'%x)+' '+str('%.4f'%y)+' '+str('%.4f'%z)
                    
                    writeCoord(msh,coords)
                    writeCoord(pts,coords)
                    
                    index += 1
    
    writeGiD_CoordEnd(msh)
  
    writeGiD_ElemBegin(msh)
    writeGiD_Elems(msh,size,radius)
    writeGiD_ElemEnd(msh)

def fanCube(xSize,ySize,zSize,cube):
  
    print 'Generating: ' + 'fanCube' + str(xSize) + 'x' + str(ySize) + 'x' + str(zSize) + '...'
        
    index = 1
    powfact = 1.1
    
    size = xSize * ySize * zSize
    radius = ((cube['top'][0]-cube['bot'][0]) / float(pow(powfact,xSize/4)) - (cube['top'][0]-cube['bot'][0]) / float(pow(powfact,xSize/4+1)))/2

    msh_filename = 'fanCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.msh'
    pts_filename = 'fanCube'+str(xSize)+'x'+str(ySize)+'x'+str(zSize)+'.'+str(int(round(radius*1000000)))+'.pts'
    
    msh=open(msh_filename, 'w+')
    pts=open(pts_filename, 'w+')
    
    writeGiD_Header(msh)
    writeGiD_CoordBegin(msh)
    writeSize(pts,size)
    
    for i in range(0,xSize):
        for j in range(0,ySize):
            for k in range(0,zSize):
                  
                x = (cube['top'][0]-cube['bot'][0]) / float(pow(powfact,xSize-i))
                y = (cube['top'][1]-cube['bot'][1]) / float(pow(powfact,ySize-j))
                z = (cube['top'][2]-cube['bot'][2]) / float(pow(powfact,zSize-k))
                
                coords = str(index)+' '+str('%.4f'%x)+' '+str('%.4f'%y)+' '+str('%.4f'%z)
                
                writeCoord(msh,coords)
                writeCoord(pts,coords)
                
                index += 1
    
    writeGiD_CoordEnd(msh)
    
    writeGiD_ElemBegin(msh)
    writeGiD_Elems(msh,size,radius)
    writeGiD_ElemEnd(msh)
    
def line(num,cube):
  
    print 'Generating: ' + 'line' + str(num) + '...'
        
    index = 1
    powfact = 1.1
    
    size = num
    radius = (2.0/num)/4.0 #Cube diagonal for 1x1x1 (TODO: Calculate this!)

    msh_filename = 'line'+str(num)+'.'+str(int(round(radius*1000000)))+'.msh'
    pts_filename = 'line'+str(num)+'.'+str(int(round(radius*1000000)))+'.pts'
    
    msh=open(msh_filename, 'w+')
    pts=open(pts_filename, 'w+')
    
    writeGiD_Header(msh)
    writeGiD_CoordBegin(msh)
    writeSize(pts,size)
    
    for i in range(0,num):
                  
        c = cube['bot'][0] + ((cube['top'][0]-cube['bot'][0])/(num-1.0)) * i
        
        coords = str(index)+' '+str('%.4f'%c)+' '+str('%.4f'%c)+' '+str('%.4f'%c)
        
        writeCoord(msh,coords)
        writeCoord(pts,coords)
                
        index += 1
    
    writeGiD_CoordEnd(msh)
    
    writeGiD_ElemBegin(msh)
    writeGiD_Elems(msh,size,radius)
    writeGiD_ElemEnd(msh)

#Main
args = sys.argv;

#Setting up default configuration
sizem = 1000
clust = 3

if(len(args) != 3):
  
    print "No parameters or wrong number of parameters, Using default configuration..."
    
else:
    
    sizem = int(args[1])
    clust = int(args[2])
    
#Normal Cube
#TODO: Merge this with offset cube
cube = {}

cube['bot'] = [0,0,0]
cube['top'] = [1,1,1]

sideI = int(round(pow(sizem,1.0/3)))

genericCube(sideI,sideI,sideI,cube)

#Cube with offsets
cubes = list()

cube1 = {}
cube2 = {}

cube1['bot'] = [0,0,0]
cube1['top'] = [0.25,0.25,0.25]

cube2['bot'] = [0.75,0.75,0.75]
cube2['top'] = [1,1,1]

cubes.append(cube1)
cubes.append(cube2)

sideI = int(round(pow(sizem/len(cubes),1.0/3)))

offsetCube(sideI,sideI,sideI,cubes)

#Cube with clusters
clusterI = list()

sideI   = clust
c_sizeI = sizem*0.999 / pow(sideI,3)
c_radiI = 0.10 / sideI

generateCluster(clusterI,c_sizeI,c_radiI)
clusterCube(sideI,sideI,sideI,cube,clusterI,c_radiI)

#Cube with fan
sideI = int(round(pow(sizem,1.0/3)))

fanCube(sideI,sideI,sideI,cube)

#Line
line(sizem,cube)

print "Done"
