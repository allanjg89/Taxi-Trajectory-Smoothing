import sys  
reload(sys)  
sys.setdefaultencoding('utf8')
import matplotlib.pyplot as plt
import matplotlib as mt
import numpy as np
from scipy import *
from scipy.linalg import solve
from datetime import datetime
import time
import pdb
import traceback
from scipy.optimize import curve_fit




'''
Parameters:
	x0: number; intitial point
	A:list of numbers
	B:list of numbers
'''	
class poly3Fit:
	def __init__(self,x0,A,B):

		'''
		A try catch block is implemented here in the case
		that a singuar matrix error arises
		'''
		try:
			x = solve(np.array(A),np.array(B))
		except Exception as err:
			print err
			raise Exception("Smoothing unsuccesful. Values not determinant.")
		
		self.x0 = x0
		self.v0 = x[0]
		self.a0 = x[1]
		self.c = x[2]
		
	'''
	Parameters:
		x0:number
		x1:number
		x2:number
		t0:number
		t1:number
		t2:number
		
	The following function does the fitting as described by smoothTrajectory
	by solving for the coeficients of the polynomial x(t) through use of algorithm 1

	For algorithm 1:
		
		Consider points p0, p1, and p2

		   B					A				 X

		|x0-x0|		|t0	  t0^2/2   t0^3/6|     |v0| 
		|x1-x0| = 	|t1	  t1^2/2   t1^3/6|  *  |a0|
		|x2-x0|		|t2	  t2^2/2   t2^3/6|	    |c|
	  


	B = A*X
	'''	
	@classmethod	
	def algorithm1(cls, x0,x1,x2,t0,t1,t2):
		A = [
				[t0,t0**2/2,t0**3/6],\
				[t1,t1**2/2,t1**3/6],\
				[t2,t2**2/2,t2**3/6]
			 ]
	
					
		B = [0,x1-x0,x2-x0]
		
		return cls(x0,A,B)
		
	'''
	Parameters:
		x0:number
		v1:number
		a1:number
		deltaT:number
			
	The following function does the fitting as described by smoothTrajectory
	by solving for the coeficients of the polynomials x(t), v(t) and a(t) 
	through use of algorithm 2


	For algorithm 2:
		
		Consider points p0, and p1

		   B					   A			 X

		|x1-x0|		|t1	  t1^2/2   t1^3/6|     |v0| 
		|v1|     = 	|0 	  t1       t1^2/2|  *  |a0|
		|a1|		|0	  1        t1    |     |c|
	  

	B = A*X
	'''	
	@classmethod 
	def algorithm2(cls,x0,x1,v1,a1,t1):
		
		
		A = [
			 [t1, t1**2/2, t1**3/6],\
			 [1,  t1,      t1**2/2],\
			 [0,  1,       t1     ]
			]
	
					
		B = [x1-x0,v1,a1]
		
		return cls(x0,A,B)

		
	def newD(self,dt):
		return self.x0 + self.v0*dt + (self.a0/2.0)*dt**2 + (self.c/6.0)*dt**3
		
	def newV(self,dt):
		return self.v0 + self.a0*dt + (self.c/2.0)*dt**2
	
	def newA(self,dt):
		return self.a0 + self.c*dt
		

'''
Parameters:
	A:list of numbers
	B:list of numbers
	
	
	note that due to the limitations of numerical representaions very large values of t 
	cannot be computed
	
	
	
In matrix form:
	
				B							                   A                                    X
	
	|0                        |   |t0        t0^{2}    t0^{3}    t0^{4}    t0^{5}    t0^{6}   |   |b|
	|- a0*t0                  |   |t0^{2}/2  t0^{3}/3  t0^{4}/4  t0^{5}/5  t0^{6}/6  t0^{7}/7 |   |c|
	|- v0*t0 - a0/2*t0        | = |t0^{3}/6  t0^{4}/12 t0^{5}/20 t0^{6}/30 t0^{7}/42 t0^{8}/56| * |d|
	|a1 - a0                  |   |t1        t1^{2}    t1^{3}    t1^{4}    t1^{5}    t1^{6}   |   |e|
	|v1 - v0 - a0*t1          |   |t1^{2}/2  t1^{3}/3  t1^{4}/4  t1^{5}/5  t1^{6}/6  t1^{7}/7 |   |f|
	|x1 - x0 - v0*t1 - a0/2*t1|   |t1^{3}/6  t1^{4}/12 t1^{5}/20 t1^{6}/30 t1^{7}/42 t1^{8}/56|   |g|


B = A*X
'''	

class poly8Fit:
	def __init__(self, x0,v0,a0,A,B):
	# __init__(self,ref,prevT0,x0,v0,a0,A,B):

		'''
		A try catch block is implemented here in the case
		that a singuar matrix error arises
		'''
		
		try:
			x = solve(np.array(A),np.array(B))
		except Exception as err:
			print err
			raise Exception("Smoothing unsuccesful. Values not determinant.")
		
		#self.ref = ref
		#self.prevT0 = prevT0
		self.x0 = x0
		self.v0 = v0
		self.a0 = a0
		self.b= x[0]
		self.c = x[1]
		self.d = x[2]
		self.e = x[3]
		self.f = x[4]
		self.g = x[5]
		
	'''
	Parameters:
		x0:number
		x1:number
		v0:number
		v1:number
		a0:number
		a1:number
		t0:number
		t1:number
	'''	
	@classmethod	
	def algorithm3(cls, x0,x1,v0,v1,a0,a1,t0,t1):
		#tOff = 0
		#ref = 0
		#prevT0 = t0 - tOff
		#prevT0 = 0
		#t1 = t1-t0+tOff
		#t0 = tOff
		#initial = 0
		
		
		A = [
				[t0,       t0**2,    t0**3,    t0**4,    t0**5,     t0**6],\
				[t0**2/2,  t0**3/3,  t0**4/4,  t0**5/5,  t0**6/6,   t0**7/7],\
				[t0**3/6,  t0**4/12, t0**5/20, t0**6/30, t0**7/42,  t0**8/56],\
				[t1,       t1**2,    t1**3,    t1**4,    t1**5,     t1**6],\
				[t1**2/2,  t1**3/3,  t1**4/4,  t1**5/5,  t1**6/6,   t1**7/7],\
				[t1**3/6,  t1**4/12, t1**5/20, t1**6/30, t1**7/42,  t1**8/56]			 ]
	
		
		'''
		B =[initial + ref,\
			initial + ref - a0*t0,\
			initial + ref - v0*t0 - a0/2*t0**2,\
			a1 + ref - a0*t0,\
			v1 + ref - v0 - a0*t1,\
			x1 + ref - x0 - v0*t1 - a0/2*t1**2]
		'''
		
		B =[0,\
			0 - a0*t0,\
			0 - v0*t0 - a0/2*t0**2,\
			a1 - a0*t0,\
			v1 - v0 - a0*t1,\
			x1 - x0 - v0*t1 - a0/2*t1**2]

		
		return cls(x0,v0,a0,A,B)
		#cls(ref,prevT0,x0,v0,a0,A,B)
		
			
	def newD(self,t):
		#t = t-self.prevT0
		return self.x0 + self.v0*t +  self.a0/2*t**2 + self.b/6*t**3 + self.c/12*t**4 + self.d/20*t**5 +self.e/30*t**6 +self.f/42*t**7 + self.g/56*t**8
		#return self.x0-self.ref + self.v0*t +  self.a0/2*t**2 + self.b/6*t**3 + self.c/12*t**4 + self.d/20*t**5 +self.e/30*t**6 +self.f/42*t**7 + self.g/56*t**8
		
	def newV(self,t):
		#t = t-self.prevT0
		return self.v0 + self.a0*t + self.b/2*t**2 + self.c/3*t**3 +self.d/4*t**4 +self.e/5*t**5 +self.f/6*t**6 + self.g/7*t**7
		#return self.v0-self.ref + self.a0*t + self.b/2*t**2 + self.c/3*t**3 +self.d/4*t**4 +self.e/5*t**5 +self.f/6*t**6 + self.g/7*t**7
	
	def newA(self,t):
		#t = t-self.prevT0
		return  self.a0 + self.b*t + self.c*t**2 + self.d*t**3 + self.e*t**4 + self.f*t**5 + self.g*t**6
		#return  self.a0-self.ref + self.b*t + self.c*t**2 + self.d*t**3 + self.e*t**4 + self.f*t**5 + self.g*t**6
		

'''
The following is an attempt to use a least squares parameter optimization
to the position, velocity and acceleration function assuming they take the following 
form:
	a(t) = a0 + b*t +c*t^2
	v(t) = v0 + a0*t + (b/2)*t^2 + (c/3)*t^3
	x(t) = x0 + v0*t + (a0/2)*t^2 + (b/6)*t^3 + (c/12)*t^4
	
This would imply that the square difference equation would take the form:
	R2 = sum{ ( Xi+Vi+Ai - [x(ti)+v(ti)+a(ti)] )^2 }
	
	Xi, Vi, Ai are the taxi data
	
	Using scipy.optimize.curve_fit whose  inputs are a function to fit "f",
	time intervals to fit over "t" and data to fit to "y"
	
	f = x(ti)+y(ti)+a(ti)
	y = Xi+Vi+Ai

'''

'''
Parameters:
	x0:number
	v0:number
	a0:number
'''

class chi2Func:
	def __init__(self, x0,v0,a0):
		self.x0 = x0
		self.v0 = v0
		self.a0 = a0
	'''
	def func5(self,t,b,c,d):
		return self.x0 + self.v0 + self.a0 + t*(self.v0+self.a0) + self.a0/2*t**2 + b*(t**3/6+t**2/2+t) + c*(t**4/12+t**3/3+t**2) + d*(t**5/20+t**4/4+t**3)
	def func4(self,t,b,c):
		return self.x0 + self.v0 + self.a0 + t*(self.v0+self.a0) + self.a0/2*t**2 + b*(t**3/6+t**2/2+t) + c*(t**4/12+t**3/3+t**2)
	'''
	def func3(self,t,b):
		return self.x0 + self.v0 + self.a0 + t*(self.v0+self.a0) + self.a0/2*t**2 + b*(t**3/6+t**2/2+t)
    
 
'''
Parameters:
	t:list of numbers; time values corresponding to the data
	x:list of numbers; x positon values
	v:list of numbers; velocity values
	a:list of numbers; acceleration values
'''   
class myChi2Fit:
	def __init__(self,x,v,a,t):
		self.x0 = x[0]
		self.v0 = v[0]
		self.a0 = a[0]
		
		func = chi2Func(self.x0,self.v0,self.a0)
		
		popt, pcov = curve_fit(func.func3, t, [xi+vi+ai for xi,vi,ai in zip(x,v,a)])
		
		self.b = popt[0]
		#self.c = popt[1]
		#self.d = popt[2]

		
	def newD(self,t):
		return self.x0 + self.v0*t + (self.a0/2.0)*t**2 + (self.b/6.0)*t**3 #+ (self.c/12)*t**4 + (self.d/30)*t**5
			
	def newV(self,t):
		return self.v0 + self.a0*t + (self.b/2.0)*t**2 #+ (self.c/3)*t**3 + (self.d/4)*t**4
		
	def newA(self,t):
		return self.a0 + self.b*t #+ self.c*t**2 + self.d*t**3


		
		
		
			