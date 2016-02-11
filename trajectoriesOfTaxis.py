
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
from scipy.interpolate import UnivariateSpline
from fitFunctions import *
import cPickle


'''
Parameters:
	taxiID:int
	times: list of numbers: units in seconds
	x: list of numbers: units in meters
	y: list of numbers: units in meters
	deltaT: list of numbers: units seconds 
	deltaX: list of numbers: units meters 
	deltaY: list of numbers: units meters 
	velocityX: list of numbers: units meters/second
	velocityY: list of numbers: units meters/second
	deltaVx: list of numbers: units meters/second
	deltaVy: list of numbers: units meters/second
	accelerationX: list of numbers: units meters/second^2
	accelerationy: list of numbers: units meters/second^2
'''

#The following two quantities are for use in determining if any 
#given data point is valid.
MAXVELOCITY = 100 #m/s
MAXACCELERATION = 10 #m/s^2

class Taxi:
	#arrOfValues = [taxiID,times,x,y,deltaT,deltaX,deltaY,velocityX,deltaVx,velocityY,deltaVy,
	#				accelerationX, accelerationY]
	def __init__(self,arrayOfValues):
		if(len(arrayOfValues) != 13):
			raise Exception("Taxi instantiation error:\nLength of entry in constructor is not equal to 13\n")
		if(len(arrayOfValues[1]) != len(arrayOfValues[2]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[3]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[4]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[5]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[6]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[7]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[8]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[9]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[10]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[11]) or\
			len(arrayOfValues[1]) != len(arrayOfValues[12])):
				raise Exception("Taxi instantiation error:\nLength of the elements in entry for constructor do not match\n")
				
			
		self.taxiID = arrayOfValues[0]
		self.times = arrayOfValues[1]
		self.x = arrayOfValues[2]
		self.y = arrayOfValues[3]
		self.deltaT = arrayOfValues[4]
		self.deltaX = arrayOfValues[5]
		self.deltaY = arrayOfValues[6]
		
		self.velocityX = arrayOfValues[7]
		self.deltaVx = arrayOfValues[8]
		self.velocityY = arrayOfValues[9]
		self.deltaVy = arrayOfValues[10]
		 
		self.accelerationX = arrayOfValues[11]
		self.accelerationY = arrayOfValues[12]
		
		#the following saves to pickle file for future loading purposes
		file = open("./processedTaxis/pickledTaxi"+str(self.taxiID)+".plk","w")
		cPickle.dump(self,file)
		file.close()
		
	def __eq__(self, other) : 
		return self.__dict__ == other.__dict__
		
		
	
	@classmethod
	def fromFile(cls,file):
		taxiID,times,x,y = Taxi.readTaxiFile(file)
		#transforms x and y to meters relative to x[0] and y[0]
		Taxi.transformXY(x,y)
		times = [Taxi.convertDate(d) for d in times]
		
		return cls(Taxi.filterQuantities(taxiID,times,x,y,file))
	
	@classmethod
	def fromPickle(cls,file):
		try:
			return cPickle.load(open(file,'r'))
		except Exception as err:
			traceback.print_exc()
			print err

		
		
	@classmethod	
	def fromBasicQuantities(cls,taxiID,t,x,y):
		return cls(Taxi.filterQuantities(taxiID,times,x,y))

		
		
	@classmethod
	def default(cls):
		return cls([-1,[1E-10],[1E-10],[1E-10],[1E-10],\
		[1E-10],[1E-10],[1E-10],[1E-10],[1E-10],[1E-10],[1E-10],[1E-10]])
		
	@classmethod
	def smoothTaxi(cls,id,t,x,y,vx,vy,ax,ay):
		initialVal = [1E-10]
		deltaT = initialVal + list(np.diff(t))
		deltaX = initialVal + list(np.diff(x))
		deltaY = initialVal + list(np.diff(y))
		deltaVx = initialVal + list(np.diff(vx))
		deltaVy = initialVal + list(np.diff(vy))
				
		return cls([id,t,x,y,deltaT,deltaX,deltaY,vx,\
					deltaVx,vy,deltaVy,ax,ay])
					
	@classmethod
	def subTaxi(cls, taxi, timeRange):
		indices = findIndicesForValues(taxi.times, timeRange)
		s = indices[0]
		f = indices[1]
		
		taxiID = taxi.taxiID
		times = taxi.times[s:f]
		x = taxi.x[s:f]
		y = taxi.y[s:f]
		deltaT = taxi.deltaT[s:f]
		deltaX = taxi.deltaX[s:f]
		deltaY = taxi.deltaY[s:f]
		
		velocityX = taxi.velocityX[s:f]
		deltaVx = taxi.deltaVx[s:f]
		velocityY = taxi.velocityY[s:f]
		deltaVy = taxi.deltaVy[s:f]
		 
		accelerationX = taxi.accelerationX[s:f]
		accelerationY = taxi.accelerationY[s:f]
		
		return Taxi([taxiID,times,x,y,deltaT,deltaX,deltaY,velocityX,deltaVx,velocityY,deltaVy,accelerationX,accelerationY])
		
		
	#The following function returns an integer representation of the date and time
	@classmethod
	def convertDate(cls,date):
	    d = datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
	    return time.mktime(d.timetuple())
	
	
	#The following is to get position measurements in meters
	#---------------------------------------------------------------------------------
	
	@classmethod
	def asRadians(cls,degrees):
	    return degrees * np.pi / 180
	
	#http://stackoverflow.com/questions/3024404/transform-longitude-latitude-into-meters
	'''
	For latitude:
		Given that 360 degrees is  a full circle around the earth through the poles,
		and that distance is 40008000 meters.  
		
			y = deltaLatitude * 40008000 / 360
			
	For Longitude:
		There is dependence on latitude since the circumference of any cross section
		at a given latitude is different. The circumference at the equator 
		(latitude 0) is 40075160 meters. The 
		The circumference of a circle at a given latitude will be proportional to the 
		cosine of that latitude, so the formula will be 
			
			x = (deltaLongitude * 40075160 / 360) * cos(latitude) 
	'''
	@classmethod
	def getPos(cls,nullLatitude, nullLongitude, latitude, longitude):
	    deltaLatitude = latitude - nullLatitude
	    deltaLongitude = longitude - nullLongitude
	    latitudeCircumference = 40075160 * np.cos(Taxi.asRadians(nullLatitude))
	    resultX = deltaLongitude * latitudeCircumference / 360
	    resultY = deltaLatitude * 40008000 / 360
	    return resultX, resultY
	
	@classmethod
	def transformXY(cls,x,y):
		nullLatitude = x[0]
		nullLongitude = y[0]
		for i in range(len(x)):
			newX,newY = Taxi.getPos(nullLatitude,nullLongitude,x[i],y[i])
			x[i] = newX
			y[i] = newY
	#---------------------------------------------------------------------------------
	'''
	The folowing function removes bad data points and generates derived 
	quantities for a a particualr set of taxi data. 
	
	Given points p1, p2, and p3:
		
		deltaT_12 = t2-t1
		deltaX_12 = x2-x1
		velocityX_12 =  deltaX_12/deltaT_12; note that this would be 
											  the velocity at p2
		deltaVx_23 = veolcityX_23-velocityX_12
		accelerationX_23 = deltaVx_23/deltaT_23^2; note that this would 
												be the acceleration at p3 
		
		Bad Data points are defined in the following way:
			
			1) t2 = t1 (identical time stamps)
			2) t2 < t1 (decreasing time)
			3) velocity is greater than some predefined velocity;
			this condition would correspond to points that moved 
			and erroneous distance in s short amount of time.
			4) acceleration is greater than some predefined value.
			This is a bad point for analogous reasons to 3)

			
		The "bad" data point along with the data point right before
		are written to a file named badPointsForTaxi{taxiID}.txt.
		The point right before is also written because
		the bad point was determined to be bad only context to the 
		point that came right before. Refer to the criteria for 
		bad points right above.
			
		
	The function returns a list of lists containing the derived 
	and filtered quantities
	'''
	@classmethod		
	def filterQuantities(cls,taxiID,times,x,y,file=None):
			
		t = [1E-10]; xF = [1E-10]; yF = [1E-10] ;deltaT = [1E-10]
		deltaX = [1E-10]; deltaY = [1E-10]; velocityX = [1E-10]
		velocityY = [1E-10]; deltaVx = [1E-10]; deltaVy = [1E-10]
		accelerationX = [1E-10]; accelerationY = [1E-10]
		
		prevI = 0
		badPointIndices = []
		badPointCount = 0
		i = 1
		initialBadCount = 0 
		upperBadCount = 3
		'''
							#initialBadCount is a measure that is meant for the 
							case where the first data point in the file is actually 
							a bad point. Such a case could cause a miss 
							classification of every other point in the 
							file as bad. Thus by keeping count of how many 
							bad points where registered at the beginning with 
							respect to this first point one can set an upper 
							bound which would cause the loop to restart with 
							the next data point (i=1) in the file, essentially 
							classifying the first data point as bad.
		'''
		
		while i < len(x):
			currT = times[i]-times[0]
			
			if currT != t[prevI] and currT > t[prevI]:
					
				dt = currT - t[prevI]
				dx = x[i] - xF[prevI]
				dy = y[i] - yF[prevI]
				vx = dx/dt
				vy = dy/dt
				
				
				if np.sqrt(vx**2+vy**2) <  MAXVELOCITY: 
					dVx = vx-velocityX[prevI]
					dVy = vy-velocityY[prevI]
					ax = dVx/dt**2
					ay = dVy/dt**2
					
					if np.sqrt(ax**2+ay**2) < MAXACCELERATION:
						initialBadCount = 0
						t.append(currT); xF.append(x[i]); yF.append(y[i])
						deltaT.append(dt); deltaX.append(dx); deltaY.append(dy)
						velocityX.append(vx); velocityY.append(vy); deltaVx.append(dVx)
						deltaVy.append(dVy); accelerationX.append(ax); accelerationY.append(ay)
						prevI += 1
			else:
				initialBadCount +=1
				badPointCount += 1
				badPointIndices.append(i)
			
			if initialBadCount == upperBadCount:
				i -= upperBadCount-1
				initialBadCount = 0
				badPointIndices = badPointIndices[0:len(badPointIndices)-upperBadCount+1]
				badPointIndices.append(i-upperBadCount+1)
			i+=1
				
				
				
		badPointFile = "./badData/badPointsForTaxi%d.txt"%taxiID
		msg = "%d of %d data points where rejected.\n"+\
				"These bad points are written to file:\n %s\n"
		print msg % (badPointCount,len(x),badPointFile)

		
		#========================================================		
		#Writing to the files Here
		#========================================================		
		fout = open(badPointFile, "w")
		
		if file:
			count = 0 #keeps track of the current line being read
			next = 0 #keeps track of teh current bad indice being considered
			fin = open(file)
			prevL = fin.readline()
			for l in fin.readlines():
				if count == badPointIndices[next]:
					next += 1; fout.write(prevL); fout.write(l)
				if next >= len(badPointIndices):
						break
				prevL = l
				count += 1
		else:
			for i in badPointIndices:
				line0 = "%d\t%d\t%d\t%d\n"%taxiID,times[i-1],x[i-1],y[i-1]
				line1 = "%d\t%d\t%d\t%d\n"%taxiID,times[i],x[i],y[i]
				fout.write(line)
						
		fout.close()
		#=============================================================

				
		return [taxiID,t,xF,yF,deltaT,deltaX,deltaY,velocityX,\
			deltaVx,velocityY,deltaVy,accelerationX,accelerationY]
	
	
		
	'''
	Parameters:
		file:string; path to file
		
	Reads in a taxi file and returns the folowing values:
		taxiID:int
		times:list of strings
		longitute: list of numbers
		latitude: list of numbers
	'''
	@classmethod
	def readTaxiFile(cls,file):
		fin = open(file)
		times = []
		longitude = []
		latitude = []
		for l in fin.readlines():
			arrValues = l.split(',')
			if len(arrValues) != 4:
				raise Exception("Data for file "+file+\
						" is corrupted. One or more lines contains"+\
						 "an incorrect number of entries.\n")
						
			times.append(arrValues[1])
			longitude.append(float(arrValues[2]))
			latitude.append(float(arrValues[3]))
			
		taxiID = int(arrValues[0]) 
		
		return taxiID,times,longitude,latitude
	
		
		
		
'''
Parameters:      
	taxi: taxi object 
	
	keyword arguments:
		desiredDeltaT:number; desired time difference between two points
		algorithm: string; "alg1" or "alg2" for algorithm 1 an 2 respectively
		match: string: "match" or "no"; if set to "match" the first entries in 
			   the x and y position in every fit will be set to x0 and y0 respectively
	
	calling example:
		smoothTaxi2Fit = smoothTrajectory(taxi,desiredDeltaT=1,algorithm ="alg2", match = "yes" )

Returns a new taxi object


The fundemental assumptions are:
	1) The data inputed is free of error such that the postions and time 
	intervals are in fact the true values
	2)The acceleration between any two points is linear in t such that
	a(t) = a0 + c*t
	v(t) = v0 + a0*t + (c/2)*t^2
	x(t) = x0 + v0*t + (a0/2)*t^2 + (c/6)*t^3
	
algorithm 1:
	
	To smooth the trajectory we will fit every 3 points to a cubic funtion in form of x(t).
	As an exampe consider points p1, p2, p3, and p4:
		iteration 1: fit to p1, p2 and p3; generate more points between p1 and p2 with fit
		iteration 2: fit to p2, p3 and p4; generate more points between p2 and p3 with fit
		
		continue in this fashion...
		
	The motivation for applying the fit to only two of the points and not all three 
	is influenced by the mechanics of driving. The velocity of the vehicle
	is determined by the driver not by a natural phenomenom that may 
	adhere to a nice function. Thus at face value it seems logical to refit
	at every new point
	
algorithm 2:
	To smooth the trajectory we will fit every 2 points. We will solve for the coefficients
	to the cubic function x(t), by solving the 3 sets of equations above using the points as 
	boundary conditions.
	
	Consider points p1 and p2
	
	a2 = a1 + c*t1
	v2 = v1 + a1*t1+ (c/2)*t1^2
	x2 = x1 + v1*t1 + (a1/2)*dt1^2 + (c/6)*dt1^3
	
algorithm 3:
	
	algorithm 3  is an attempt to fit the highest degree polynomial given
	all the data (ie boundary conditions) we have at our disposal for 2 points.
		
	NOTE: For very large values of time the numerical limits of computing do not allow 
	for satisfactory fits. The time values could be normalize with respect to the first 
	point (i.e. t1 = t1-10; t0 ~= 0) but from experimentation this results in a very 
	poor fit to the velocity and acceleration when compared to un-normalized time values. 
	This poor fit when normalizing could be due to some error I failed to identify or
	to the limited precision fo arithmetic operations done numerically.
		
	The equations of motion are:
		a(t) = a0 + b*t  + c*t^{2} + d*t^{3} +e*t^{4} +f*t^{5} g*t^{6}
		v(t) = v0 + a0*t + b/2*t^{2} + c/3*t^{3} + d/4*t^{4} +e/5*t^{5} +f/6*t^{6} g/7*t^{7}
		x(t) = x0 + v0*t + a0/2*t^{2} + b/6*t^{3} + c/12*t^{4} + d/20*t^{5} +e/30*t^{6} +f/42*t^{7} + g/56*t^{8}

	Consider points p0 and p1:
		a0 - a0                       = b*t0 + c*t0^{2} + d*t0^{3} +e*t0^{4} +f*t0^{5} g*t0^{6}
		v0 - v0 - a0*t0               = b/2*t0^{2} + c/3*t0^{3} + d/4*t0^{4} +e/5*t0^{5} +f/6*t0^{6} + g/7*t0^{7}
		x0 - x0 - v0*t0 - a0/2*t0^{2} = b/6*t0^{3} + c/12*t0^{4} + d/20*t0^{5} +e/30*t0^{6} +f/42*t0^{7} + g/56*t0^{8}
		
		a1 - a0                       = b*t1 + c*t1^{2} + d*t1^{3} +e*t1^{4} +f*t1^{5} g*t1^{6}
		v1 - v0 - a0*t1               = b/2*t1^{2} + c/3*t1^{3} + d/4*t1^{4} +e/5*t1^{5} +f/6*t1^{6} + g/7*t1^{7}
		x1 - x0 - v0*t1 - a0/2*t1^{2} = b/6*t1^{3} + c/12*t1^{4} + d/20*t1^{5} +e/30*t1^{6} +f/42*t1^{7} + g/56*t1^{8}


algorithm 4:
	
	Algorithm 4 is an attempt to use a least squares parameter optimization
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
		
		f_i = x(ti)+y(ti)+a(ti)
		y_i = Xi+Vi+Ai


'''	
def smoothTrajectory(taxi,**kwargs): 

	deltaT = taxi.deltaT
	x = taxi.x
	y = taxi.y
	vx = taxi.velocityX
	vy = taxi.velocityY
	ax = taxi.accelerationX
	ay = taxi.accelerationY
	t = taxi.times
	
	if "desiredDeltaT" in kwargs:
		desiredDeltaT = kwargs["desiredDeltaT"]
	else:
		desiredDeltaT= np.mean(deltaT)/5.0
					
	if "match" in kwargs:
		match = kwargs["match"]
	else:
		match = "no"
		
	if "algorithm" in kwargs:
		algorithm = kwargs["algorithm"]
	else:
		algorithm = "alg1"

	#In the case that the desired time step is greater that or 
	#equal to the maximum time step the function simply returns 
	#the inputed taxi object
	if desiredDeltaT>=max(taxi.deltaT):
		return taxi
		
	n = 2 #the number of points to inlcude in the chi squared fit in the case of algorithm 4

	smoothT = []; smoothX = []; smoothY = []; smoothVx = []; smoothVy = []
	smoothAx = []; smoothAy = []; 			
	for i in range(len(x)-1):
			
			try:
				#if i<len(x)-3:
				if algorithm == "alg1" and i<len(x)-3:
					totT = t[i]
					funcsx = poly3Fit.algorithm1(x[i],x[i+1],x[i+2],t[i],t[i+1],t[i+2])
					funcsy = poly3Fit.algorithm1(y[i],y[i+1],y[i+2],t[i],t[i+1],t[i+2])
				elif algorithm == "alg2" and i<len(x)-2:
					i +=2 #we are starting with the third element since the numerically 
						  #computed acceleration is not well defined before that. 
					totT = t[i]
					funcsx = poly3Fit.algorithm2(x[i],x[i+1],vx[i+1],ax[i+1],t[i+1])
					funcsy = poly3Fit.algorithm2(y[i],y[i+1],vy[i+1],ay[i+1],t[i+1])
				elif algorithm == "alg3"and i<len(x)-2:
					i +=2 #we are starting with the third element since the numerically 
						  #computed acceleration is not well defined before that. 
					totT = t[i]
					funcsx = poly8Fit.algorithm3(x[i],x[i+1],vx[i],vx[i+1],ax[i],ax[i+1],t[i],t[i+1])
					funcsy = poly8Fit.algorithm3(y[i],y[i+1],vy[i],vy[i+1],ay[i],ay[i+1],t[i],t[i+1])
				elif algorithm == "alg4" and i<len(x)-(n+1):
					i +=2 #we are starting with the third element since the numerically 
						  #computed acceleration is not well defined before that. 
					totT = t[i]
					funcsx = myChi2Fit(x[i:i+n],vx[i:i+n],ax[i:i+n],t[i:i+n])
					funcsy = myChi2Fit(y[i:i+n],vy[i:i+n],ay[i:i+n],t[i:i+n])

				else:
					raise Exception("algorithm not properly specified. PLease specify:\n alg1, alg2, alg3, or alg4")
				
				count = 0
				while (totT<t[i+1]):
					if match == "yes" and count == 0: #a measure intended for the second 
												   #algorithm to match the intial position point for 
												   #each fit 
						smoothX.append(x[i])
						smoothY.append(y[i])
					else:
						smoothX.append(funcsx.newD(totT))
						smoothY.append(funcsy.newD(totT))
						
					smoothVx.append(funcsx.newV(totT))
					smoothVy.append(funcsy.newV(totT))
					smoothAx.append(funcsx.newA(totT))
					smoothAy.append(funcsy.newA(totT))
					smoothT.append(totT)
					totT += desiredDeltaT
					count += 1
					
			except Exception as err:
				traceback.print_exc()
				print err
				
				errmsg  = "Note values corresponding to elements "+ \
						"%d, %d, and %d could not be fitted.\n" + \
						"Smooth trajectory between elements %d and %d "+\
						"could not be created.\n"
				print errmsg % (i, i+1, i+2, i, i+2)
				continue
			
						
		
	return Taxi.smoothTaxi(taxi.taxiID,smoothT,smoothX,smoothY,smoothVx,smoothVy,smoothAx,smoothAy)
	

'''
Parameters:
	taxi:Taxi Object
Optional Parameters:
	desiredDeltaT: number: specifies the time interval between 
				   consecutive data points
				
Returns:
	Taxi Object
	scipy.interpolate.UnivariateSpline Object
	scipy.interpolate.UnivariateSpline Object
	
Here an order three spline is fit to the x and y coordinates and let the 
software determine the number of knots. The velocity and acceleration 
are determined via numerical differentiation. (NOTE: the software has 
built in methods that preform differentiation of a given order.)
	

The number of knots are selected until the following condition is met
squared error < N where N is the number of points.

\sqrt{(y_{i}-\hat{f}(x_{i}))^{2}} <= N

'''


def smoothingSplineTaxi(taxi,**kwargs):
	degree = 3
	sX = UnivariateSpline(taxi.times, taxi.x, s=degree)
	sY = UnivariateSpline(taxi.times, taxi.y, s=degree)
	init = 1E-14
	
	if "desiredDeltaT" in kwargs:
		desiredDeltaT = kwargs["desiredDeltaT"]
	else:
		desiredDeltaT= np.mean(taxi.deltaT)/5.0
	
	smoothT = []; smoothX = []; smoothY = []; smoothVx = []; smoothVy = []
	smoothAx = []; smoothAy = []; 			
	for i in range(len(taxi.x)-1):
		
		totT = taxi.times[i]
		
		while (totT<taxi.times[i+1]):
			x = sX(totT)
			x2 = sX(totT+desiredDeltaT)
			x3 = sX(totT+2*desiredDeltaT)
			y = sY(totT)
			y2 = sY(totT+desiredDeltaT)
			y3 = sY(totT+2*desiredDeltaT)
			vX = (x2-x)/desiredDeltaT 
			vX2 =(x3-x2)/desiredDeltaT 
			vY = (y2-y)/desiredDeltaT
			vY2 = (y3-y2)/desiredDeltaT
			aX = (vX2-vX)/desiredDeltaT
			aY = (vY2-vY)/desiredDeltaT
			 
			smoothX.append(x)
			smoothY.append(y)
			smoothVx.append(vX)
			smoothVy.append(vY)
			smoothAx.append(aX)
			smoothAy.append(aY)
			smoothT.append(totT)
			totT += desiredDeltaT
			
	return Taxi.smoothTaxi(taxi.taxiID,smoothT,smoothX,smoothY,smoothVx,smoothVy,smoothAx,smoothAy), sX, sY

	
	
	
	
	
'''
The following computed the standard deviation of the distance intervals
'''
def stdOfDeltaD(taxi):
	return np.std([np.sqrt(x**2+y**2) for x,y in zip(taxi.deltaX,taxi.deltaY)])
	
	
#==============================================================================================
#CREATION OF FIGURES
#==============================================================================================
def createHistograms(taxis):
	allTimeIntervals = []
	allDistanceIntervals = []
	
	for taxi in taxis:
		allTimeIntervals.extend([t/60 for t in taxi.deltaT])
		allDistanceIntervals.extend([np.sqrt(x*x+y*y) for x,y in zip(taxi.deltaX,taxi.deltaY)])

	n, bins, patches = plt.hist(allTimeIntervals, bins=10, range=[0,.2], normed=1)
	plt.xlabel('minutes', fontsize=18)
	plt.ylabel('proportion', fontsize=16)
	plt.show()
	
    
	n, bins, patches = plt.hist(allDistanceIntervals, bins=10, range = [0,7], normed=1)
	plt.xlabel('meters', fontsize=18)
	plt.ylabel('proportion', fontsize=16)
	plt.show()	
    
    
def createTrajectoryPlots(taxi,smoothTaxi,rangeT=None):
	if rangeT == None:
		rangeT = [taxi.times[0],taxi.times[len(taxi.times)-1]]
	indices = findIndicesForValues(taxi.times, rangeT)
	sIndices = findIndicesForValues(smoothTaxi.times, rangeT)
	start = indices[0]
	end = indices[1]
	sStart = sIndices[0]
	sEnd = sIndices[1]
	divFactor = 5
	tstart = taxi.times[start]
	tend = taxi.times[end]
	xstart = taxi.x[start]
	xend = taxi.x[end]
	ystart = taxi.y[start]
	yend = taxi.y[end]
	t = taxi.times[start:end]
	tepsilon  = (max(t)-min(t))/divFactor
	x = taxi.x[start:end]
	xepsilon  = (max(x)-min(x))/divFactor
	y = taxi.y[start:end]
	yepsilon  = (max(y)-min(y))/divFactor
	st = smoothTaxi.times[sStart:sEnd]
	sx = smoothTaxi.x[sStart:sEnd]
	sy = smoothTaxi.y[sStart:sEnd]

	fig, ax = plt.subplots()
	ax.plot(t, x, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sx,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(x)-xepsilon,max(x)+xepsilon])
	plt.xlabel('seconds', fontsize=18)
	plt.ylabel('meters', fontsize=16)
	plt.title('X Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()

	#plot of y coordiate
	fig, ax = plt.subplots()
	ax.plot(t, y, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sy,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(y)-xepsilon,max(y)+xepsilon])
	plt.xlabel('seconds', fontsize=18)
	plt.ylabel('meters', fontsize=16)
	plt.title('Y Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()

	#plot of x(t) vs y(t)
	fig, ax = plt.subplots()
	ax.plot(x, y, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(sx, sy,'.-', label='Smoothed')
	plt.xlim([min(x)-xepsilon,max(x)+xepsilon])
	plt.ylim([min(y)-yepsilon,max(y)+yepsilon])
	plt.xlabel('meters', fontsize=18)
	plt.ylabel('meters', fontsize=16)
	plt.title('X(t) vs Y(t) Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()	
	
def createVelocityPlots(taxi,smoothTaxi,rangeT=None):
	if rangeT == None:
		rangeT = [taxi.times[0],taxi.times[len(taxi.times)-1]]
	indices = findIndicesForValues(taxi.times, rangeT)
	sIndices = findIndicesForValues(smoothTaxi.times, rangeT)
	start = indices[0]
	end = indices[1]
	sStart = sIndices[0]
	sEnd = sIndices[1]
	divFactor = 5
	tstart = taxi.times[start]
	tend = taxi.times[end]
	xstart = taxi.velocityX[start]
	xend = taxi.velocityX[end]
	ystart = taxi.velocityY[start]
	yend = taxi.velocityY[end]
	t = taxi.times[start:end]
	tepsilon  = (max(t)-min(t))/divFactor
	x = taxi.velocityX[start:end]
	xepsilon  = (max(x)-min(x))/divFactor
	y = taxi.velocityY[start:end]
	yepsilon  = (max(y)-min(y))/divFactor
	st = smoothTaxi.times[sStart:sEnd]
	sx = smoothTaxi.velocityX[sStart:sEnd]
	sy = smoothTaxi.velocityY[sStart:sEnd]

	fig, ax = plt.subplots()
	ax.plot(t, x, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sx,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(x)-xepsilon,max(x)+xepsilon])
	plt.ylabel('meters/seconds', fontsize=16)
	plt.title('Velocity of X Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()

	#plot of y coordiate
	fig, ax = plt.subplots()
	ax.plot(t, y, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sy,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(y)-xepsilon,max(y)+xepsilon])
	plt.xlabel('seconds', fontsize=18)
	plt.ylabel('meters/seconds', fontsize=16)
	plt.title('Velocity of Y Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()

	
	
	
def createAccelerationPlots(taxi,smoothTaxi,rangeT=None):
	if rangeT == None:
		rangeT = [taxi.times[0],taxi.times[len(taxi.times)-1]]
	indices = findIndicesForValues(taxi.times, rangeT)
	sIndices = findIndicesForValues(smoothTaxi.times, rangeT)
	start = indices[0]
	end = indices[1]
	sStart = sIndices[0]
	sEnd = sIndices[1]	
	divFactor = 5
	tstart = taxi.times[start]
	tend = taxi.times[end]
	xstart = taxi.accelerationX[start]
	xend = taxi.accelerationX[end]
	ystart = taxi.accelerationY[start]
	yend = taxi.accelerationY[end]
	t = taxi.times[start:end]
	tepsilon  = (max(t)-min(t))/divFactor
	x = taxi.accelerationX[start:end]
	xepsilon  = (max(x)-min(x))/divFactor
	y = taxi.accelerationY[start:end]
	yepsilon  = (max(y)-min(y))/divFactor
	st = smoothTaxi.times[sStart:sEnd]
	sx = smoothTaxi.accelerationX[sStart:sEnd]
	sy = smoothTaxi.accelerationY[sStart:sEnd]

	fig, ax = plt.subplots()
	ax.plot(t, x, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sx,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(x)-xepsilon,max(x)+xepsilon])
	plt.ylabel('meters/seconds^{2}', fontsize=16)
	plt.title('Acceleration of X Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()

	#plot of y coordiate
	fig, ax = plt.subplots()
	ax.plot(t, y, 'ro-', label='Data', fillstyle = 'none')
	ax.plot(st, sy,'.-', label='Smoothed')
	plt.xlim([min(t)-tepsilon,max(t)+tepsilon])
	plt.ylim([min(y)-xepsilon,max(y)+xepsilon])
	plt.xlabel('seconds', fontsize=18)
	plt.ylabel('meters/seconds^{2}', fontsize=16)
	plt.title('Acceleration of Y Coordinate', fontsize=20)
	legend = ax.legend(loc='best', shadow=True)
	plt.show()



	

'''
Parameters:
	x:list of numbers
	values: list of numbers
	
returns indices: list of numbers
	
The following function finds and returns the indices
within a given list x corresponding to the values in "values".
Note that that values and x must be sorted in increasing order
'''
def findIndicesForValues(x,values):
	
	indices = []
	prevV = values[0]
	indice = 0
	
	for v in values:
		if v < prevV:
			raise Exception("Invalid Values: List values in findIndiceForValues should be increasing.")
			
		theMin = abs(v-x[indice])
		start = indice+1
		
		for i in range(start,len(x)):
			currDiff = abs(v-x[i])
			if currDiff<theMin:
				theMin = currDiff
				indice = i
			prevV= v
			
		indices.append(indice)
		
	
	return indices


taxi = Taxi.fromFile("./taxiData/1131.txt")
taxi2 = Taxi.fromPickle("./processedTaxis/pickledTaxi1131.plk")

print taxi.x[0:10]
print taxi2.x[0:10]


