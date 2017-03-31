import numpy as np
#import matplotlib.pyplot as plt
import sys
from scipy.sparse import *
import math

# #First we define a function that takes in a matrix and computes the point where 
# #S(i) is minimized, splits the matrix at that position and makes two recursive 
# #calls to itself on both submatrices
# #This function will return an array of arrays, i.e. a 2D matrix
# #Each row of the matrix looks like this:
# #[position_in_bp_coordinates level S(i)_value(=Value_B) Value_A]

N_chr = int(sys.argv[3])
resolution = int(sys.argv[4])

def find_breaks(m,start_m,level,N_chr,global_mm):
	#We define the exiting condition here:
	#If the level is more than log_2(N_chr) then we should stop
	#Or if we're inside a chromosome
	break_point,m1,s1,m2,s2,S_min,Value_A,global_B = find_break_point(m,start_m,global_mm)

#3+math.log(N_chr,2)

	# print(m1.shape)
	# print(m2.shape)
	# print("aaa")

	# We run for 23 levels to ensure that we hit all chr boundaries
	# We also assume that chr sizes > 1Mb
	if( (level>=N_chr+1) or (m1.shape[0]<20) or (m2.shape[0]<20) ):
		return [break_point*resolution,level,S_min,Value_A,global_B,min(m1.shape[0],m2.shape[0])]
	#Check for being inside a chromosome(will come up with this later)	
	#elif()
	#	again 
	#	return [break_point*resolution,level,S_min,Value_A]
	else:
		return np.concatenate(([break_point*resolution,level,S_min,Value_A,global_B,min(m1.shape[0],m2.shape[0])],np.concatenate((find_breaks(m1,s1,level+1,N_chr,global_mm),find_breaks(m2,s2,level+1,N_chr,global_mm)))))




def find_break_point(m,start_m,global_mm):
	mm = np.cumsum(np.cumsum(m.toarray(),1),0)
	L = mm.shape[0]
	S = np.zeros(L)
	
	# print("The size of this matrix is",L)

	for i in range(L):
		S[i] = float(mm[i,L-1]-mm[i,i])/((i+1)*(L-i))
	
	#We don't want the first or last entry of S to be the breakpoint
	index_1 = min(20,L-1)

	S[:index_1] = [np.amax(S)]*(index_1)

	
	index_2 = min(20,L-1)
	S[-index_2:] = [np.amax(S)]*(index_2)



	# print(L)
	# print(mm)
	# print(S)

	# edge_point = int(7*L/9)
	# float(mm[int(8L/9),int(2L/9)]-mm[int()])/((break_point+1)*(L-break_point))

	#Now find the point where S is minimized
	break_point,min_val = np.argmin(S),np.amin(S) #+/- 1?
	
	#We also find the value of the region A described by IntegralImage(break_point)
	value_A = float(mm[break_point,break_point])/((break_point+1)*(break_point+1))
	min_val = float(mm[break_point,L-1]-mm[break_point,break_point])/((break_point+1)*(L-break_point))



	global_break_point = break_point+start_m
	gbp = global_break_point
	GL = global_mm.shape[0]
	global_B = float(global_mm[gbp,GL-1]-global_mm[gbp,gbp])/((gbp+1)*(GL-gbp))

	#Shifting the breakpoint forward by 1
	break_point = break_point+1
	
	# print("break_point is:",break_point)

	
	# print(break_point,min_val,value_A)


	#Now break the matrix about the break_point
	
	list1 = range(break_point)
	list2 = range(break_point,L)

	# m1 = mm[0:break_point,0:break_point]
	# m2 = mm[break_point:L,break_point:L]

	m1 = m.tocsr()[list1, :].tocsc()[:,list1]
	m2 = m.tocsr()[list2, :].tocsc()[:,list2]
	
	s1 = 0+start_m
	s2 = break_point+start_m
	break_point = break_point+start_m

	# print(mm)
	# print(break_point,np.shape(m1))
	# print(m1)
	# print(s1,np.shape(m2))
	# print(m2)
	# print(s2)

	# plt.plot(S)
	# plt.show()

	return break_point,m1,s1,m2,s2,min_val,value_A,global_B








def iter_loadtxt(filename, delimiter='\t', skiprows=0, dtype=float):
    def iter_func():
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.rstrip().split(delimiter)
                for item in line:
                    yield dtype(item)
        iter_loadtxt.rowlength = len(line)

    data = np.fromiter(iter_func(), dtype=dtype)
    data = data.reshape((-1, iter_loadtxt.rowlength))
    return data



#Load Hi-C Juicebox Dump

#i,j,value=np.loadtxt('jbdump.txt').T
data = iter_loadtxt(sys.argv[1])

i = np.asarray(data[:,0])
j = np.asarray(data[:,1])
value = np.asarray(data[:,2])

#Zero out the diagonal
for iter in range(i.shape[0]):
	if(i[iter]==j[iter]):
		value[iter]=0


#Convert into symmetric matrix with resolution level indices

i2 = [int(x/resolution) for x in i]
j2 = [int(y/resolution) for y in j]

i22 = [int(y/resolution) for y in j]
j22 = [int(x/resolution) for x in i]

value1 = [float(z) for z in value]
value2 = [float(z) for z in value]

i222 = i2 + i22
j222 = j2 + j22
value222 = value1 + value2


m2=coo_matrix((value222,(i222,j222)))

#print("The dimensions of the complete input matrix are:")
#print(m2.shape)
#print("\n")

# Take only a part of the matrix for testing
#list1 = range(40000)
#m = m2.tocsr()[list1, :].tocsc()[:,list1]

# Or Take the entire matrix for testing
m=m2

# print(m.todense())

#mm = np.cumsum(np.cumsum(m.toarray(),1),0)
# print("Integral Image computed! Now finding breaks...")

# plt.spy(m,precision=2,marker=".",markersize=1)
# plt.show()

#m = [[1,2,3,4,5,6,7],[8,9,10,11,12,13,14],[15,16,17,18,19,20,21],[22,23,24,25,26,27,28],[29,30,31,32,33,34,35],[36,37,38,39,40,41,42],[43,44,45,46,47,48,49]]

# print(find_break_point(mm,0))

# print(find_breaks(mm,0,0,23))

global_mm = np.cumsum(np.cumsum(m.toarray(),1),0)

Breaks = find_breaks(m,0,0,N_chr,global_mm);

n = int(np.shape(Breaks)[0]/6)

Breaks_Output = np.reshape(Breaks,(n,6))

Breaks_Sorted_by_Signal = Breaks_Output[Breaks_Output[:,2].argsort()]

np.set_printoptions(precision=8,threshold=sys.maxint,linewidth=250,suppress=True)

#Now we iterate through this list and choose the top 
Final_Breaks = np.zeros(N_chr-1)
curr_i = 0
m_size_threshold = 200

#Now: Potential change for the future
#m_size_threshold = 50

for i in range(n):
	if(curr_i>=N_chr-1):
		break
	if(Breaks_Sorted_by_Signal[i,5]>m_size_threshold):
		Final_Breaks[curr_i] = Breaks_Sorted_by_Signal[i,0]
		curr_i += 1

#print(np.sort(Final_Breaks))
np.savetxt(sys.argv[2],np.sort(Final_Breaks),fmt='%d 100',header="variableStep chrom=assembly span=10000",comments='')
#print(Breaks_Sorted_by_Signal)
