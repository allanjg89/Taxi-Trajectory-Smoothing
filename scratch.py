class poly11Fit:
	def __init__(self,x0,v0,a0,A,B):

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
		self.v0 = v0
		self.a0 = a0
		self.b = x[0]
		self.c = x[1]
		self.d = x[2]
		self.e = x[3]
		self.f = x[4]
		self.g = x[5]
		self.h = x[6]
		self.i = x[7]
		self.j = x[8]
		
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
	def algorythm4(cls, x0,x1,x2,v0,v1,v2,a0,a1,a2,t0,t1,t2):
		A = [  
				[t0,       t0**2,    t0**3,    t0**4,    t0**5,     t0**6,    t0**7,    t0**8,     t0**9],\
				[t0**2/2,  t0**3/3,  t0**4/4,  t0**5/5,  t0**6/6,   t0**7/7,  t0**8/8,  t0**9/9,   t0**10/10],\
				[t0**3/6,  t0**4/12, t0**5/20, t0**6/30, t0**7/42,  t0**8/56, t0**9/72, t0**10/90, t0**11/110],\
				[t1,       t1**2,    t1**3,    t1**4,    t1**5,     t1**6,    t1**7,    t1**8,     t1**9],\
				[t1**2/2,  t1**3/3,  t1**4/4,  t1**5/5,  t1**6/6,   t1**7/7,  t1**8/8,  t1**9/9,   t1**10/10],\
				[t1**3/6,  t1**4/12, t1**5/20, t1**6/30, t1**7/42,  t1**8/56, t1**9/72, t1**10/90, t1**11/110],\
				[t2,       t2**2,    t2**3,    t2**4,    t2**5,     t2**6,    t2**7,    t2**8,     t2**9],\
				[t2**2/2,  t2**3/3,  t2**4/4,  t2**5/5,  t2**6/6,   t2**7/7,  t2**8/8,  t2**9/9,   t2**10/10],\
				[t2**3/6,  t2**4/12, t2**5/20, t2**6/30, t2**7/42,  t2**8/56, t2**9/72, t2**10/90, t2**11/110]
			 ]
	
		
		B =[0,\
			0- a0*t0,\
			0 - v0*t0 - a0/2*t0**2,\
			a1 - a0,\
			v1 - v0 - a0*t1,\
			x1 - x0 - v0*t1 - a0/2*t1**2,
			a2 - a0,\
			v2 - v0 - a0*t1,\
			x2 - x0 - v0*t1 - a0/2*t1**2]
		
		return cls(x0,v0,a0,A,B)
		
			
	def newD(self,t):
		return self.x0 + self.v0*t +  self.a0/2*t**2 + self.b/6*t**3 + self.c/12*t**4 + self.d/20*t**5 +self.e/30*t**6 +self.f/42*t**7 + self.g/56*t**8 + self.h/72*t**9 + self.i/90*t**10 + self.j/110*t**11
		
	def newV(self,t):
		return self.v0 + self.a0*t + self.b/2*t**2 + self.c/3*t**3 +self.d/4*t**4 +self.e/5*t**5 +self.f/6*t**6 + self.g/7*t**7 +  self.h/8*t**8 + self.i/9*t**9 + self.j/10*t**10 
	
	def newA(self,t):
		return  self.a0 + self.b*t + self.c*t**2 + self.d*t**3 + self.e*t**4 + self.f*t**5 + self.g*t**6 + self.h*t**7 + self.i*t**8 + self.j*t**9
		
		
		
		
elif algorythm == "alg4":
						i +=2 #we are starting with the third element since the numerically 
							  #computed acceleration is not well defined before that. 
						totT = t[i]
						funcsx = poly11Fit.algorythm4(x[i],x[i+1],x[i+2],vx[i],vx[i+1],vx[i+2],ax[i],ax[i+1],ax[i+2],t[i],t[i+1],t[i+2])
						funcsy = poly11Fit.algorythm4(y[i],y[i+1],y[i+2],vy[i],vy[i+1],vy[i+2],ay[i],ay[i+1],ay[i+2],t[i],t[i+1],t[i+2])
