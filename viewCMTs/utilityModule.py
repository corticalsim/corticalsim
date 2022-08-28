import numpy
import warnings 
warnings.filterwarnings("ignore")

def readrow(row, cols):
    a = numpy.fromstring(row, sep=' ') 
    return a

def sphereFit(spX,spY,spZ):
    spX = numpy.array(spX)
    spY = numpy.array(spY)
    spZ = numpy.array(spZ)
    A = numpy.zeros((len(spX),4))
    A[:,0] = spX*2
    A[:,1] = spY*2
    A[:,2] = spZ*2
    A[:,3] = 1
    f = numpy.zeros((len(spX),1))
    f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
    C, residules, rank, singval = numpy.linalg.lstsq(A,f)
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = numpy.sqrt(t)
    center = [C[0],C[1],C[2]]
    return radius,center
    
def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]    
   


