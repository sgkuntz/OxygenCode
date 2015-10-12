# Curve-fitting for surface: 
# try multiple variable fitting
# options: 'graph' and 'compare'

# python CurveFitting_v5.py Me R303 17.5 22.5 27.5 graph

import math, os, re, string, sys
import matplotlib as mpl
import pylab as pyl
import numpy as np
from numpy.linalg import inv
import scipy.stats as st
	
def surfaceRegression(mat,tmat): #mat is [1 T O . . . ], and tmat is [t]
	# invMat = inv(mat)
	transMat = mat.transpose()
	productMat = np.dot(transMat,mat)
	invMat = inv(productMat)
	productMat2 = np.dot(invMat,transMat)
	betamat = np.dot(productMat2,tmat)
	return betamat	
			
def calcRsqValues(yMatrix,xMatrix,bMatrix,meany): # y is output, x is input, b is surface
	SStot, SSreg, SSres = 0,0,0			
	for value in range(len(yMatrix)):
		timing = yMatrix[value]	
		SStot += (timing-meany)**2
		estimatedTiming = 0
		for variable in range(len(bMatrix)):
			estimatedTiming += bMatrix[variable]*xMatrix[value][variable]
		SSres += (timing-estimatedTiming)**2
		SSreg += (estimatedTiming - meanTiming)**2
	n = len(yMatrix)
	p = len(bMatrix)			
	# print 'SSres=',SSres, 'SSreg=',SSreg, 'SStot=',SStot
	Rsquared = 1 - ((SSres)/(SStot))
	Rsquared2 = SSreg/SStot
	Rbarsquared = 1 - ((SSres*(n-1))/(SStot*(n--p-1)))
	return Rsquared, Rsquared2, Rbarsquared

def changeName(i):
	for (shortName, fullName) in speciesNames:
		if i == shortName:
			return fullName
def changeTemp(j):
	if len(j) == 4: return j
	elif len(j) == 2: return j+'.0'

speciesNames, speciesOrdering = [],[]			
speciesNamesFile = open('speciesNamesList.txt','r')
for j in speciesNamesFile:
	abbreviatedName, traditionalName, filler = re.split("[\s]",j)
	speciesNames.append((abbreviatedName, traditionalName))
	speciesOrdering.append(traditionalName)	
	
temperatureListDefault = ['17.5','22.5','27.5']

SpList = []; Tlist = []; Olist = []
for species in speciesOrdering:
	for i in sys.argv:
		if changeName(i) == species:
			SpList.append(changeName(i))
for i in sys.argv:
	if i in temperatureListDefault:
		Tlist.append(changeTemp(i))
	elif i[2:3] == '%':
		Olist.append(i[:2])

comboDict = {}
animalList = []

timeLapseInfoFile = open('TimeLapseInstructions_run_trial15.txt','r')
for i in timeLapseInfoFile:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temp, scaling, zStack, filler] = re.split("[\s]",i)
	if mryTime <> '' and tfTime <> '':
		if changeName(strain) in SpList:
			speciesFull = changeName(strain) # change the species name
			if changeTemp(temp) in Tlist:
				temperature = changeTemp(temp)
				if date+position not in animalList:
					animalList.append(date+position)
					if temperature not in comboDict:
						comboDict[temperature] = {}
					if strain[:2] not in comboDict[temperature]:
						comboDict[temperature][strain[:2]] = []
					comboDict[temperature][strain[:2]].append((float(scaling)*(int(tfTime)-int(mryTime)))/60) # record time in hours

## Curve-fitting
print comboDict.keys()
print Olist

outputFile = open('curveFitting/Curvefitting_outputData_'+sys.argv[1]+'.txt','w')
fig = pyl.figure(figsize=(25.5,33))
timingSum = 0
inputArray, outputArray = [], []
inputArrayLarge, inputArrayLarge2, outputArrayLarge = [], [], []
for temperature in Tlist: #species will include oxygen concentration
	for oxygen in Olist:
		for timing in comboDict[temperature][oxygen]:
			inputArray.append((1.,1/float(temperature),1/float(oxygen)))
			# inputArrayLarge.append((1.,1/float(temperature),1/float(oxygen),1/(float(temperature)*float(oxygen)),1/(float(temperature)**2),1/(float(oxygen)**2),1/(float(temperature)**3),1/(float(oxygen)**3)))
			# inputArrayLarge.append((math.log(1/(float(oxygen))),1/float(temperature)))
			##inputArrayLarge.append((1.,1./float(temperature)))
			inputArrayLarge2.append((1.,1./float(temperature),1./float(oxygen)))
			inputArrayLarge.append((1.,math.exp(1./float(temperature)), math.exp(1./float(oxygen)))) #, math.exp(1./float(temperature))/float(oxygen)))
			outputArray.append(math.log(timing))
			##outputArrayLarge.append(math.log(timing)+math.log(float(oxygen)))
			outputArrayLarge.append(timing)
			timingSum += timing
			
meanTiming = timingSum/len(outputArray)			
inputMatrix = np.array(inputArray)
outputMatrix = np.array(outputArray)
inputMatrixLarge = np.array(inputArrayLarge)
inputMatrixLarge2 = np.array(inputArrayLarge2)
outputMatrixLarge = np.array(outputArrayLarge)

surfaceMatrix = surfaceRegression(inputMatrix,outputMatrix)
surfaceMatrixExpanded = surfaceRegression(inputMatrixLarge,outputMatrixLarge)
surfaceMatrixExpanded2 = surfaceRegression(inputMatrixLarge2,outputMatrixLarge)

(Rsq, Rsq2, RbarSq) = calcRsqValues(outputMatrix,inputMatrix,surfaceMatrix,meanTiming)
(RsqL, Rsq2L, RbarSqL) = calcRsqValues(outputMatrix,inputMatrixLarge,surfaceMatrixExpanded,meanTiming)

print 'simple=', Rsq, Rsq2, RbarSq
print 'long=',RsqL, Rsq2L, RbarSqL
print surfaceMatrix, surfaceMatrixExpanded, surfaceMatrixExpanded2


sanityCheckFile = open('curveFitting/curveFitting3D_'+sys.argv[1]+'sanityCheck.txt','w')
for value in range(len(inputMatrix)):
	actualTiming = outputMatrix[value] #log(t)
	exponentTiming = math.exp(outputMatrix[value]) # t
	calculatedTiming = math.exp(surfaceMatrix[0])*math.exp(inputMatrix[value][1]*surfaceMatrix[1] + inputMatrix[value][2]*surfaceMatrix[2]) # est. t
	sanityCheckFile.write(str(actualTiming)+'\t'+str(exponentTiming)+'\t'+str(calculatedTiming)+'\n')

if 'stats' in sys.argv:
	for temperature in Tlist:
		# z = np.arange(17, 28.5,.1)
		z = np.arange(10, 30,.1)
		w, w1, w2 = [], [], []
		confidenceGap50pos1, confidenceGap95pos1, confidenceGap50neg1, confidenceGap95neg1 = [], [], [], []
		talpha = st.t.isf(0.05/2,(n-2)) # isf is one-tailed, so divide by 2 for 2-tailed critical value
		tbeta = st.t.isf(0.5/2,(n-2))
		standardErrorOfTheEstimate1 = math.sqrt((findVariance(y1Values)-findVariance(x1Values)*b1**2)*(n-1)/(n-2))
		for u in z:
			w1.append(math.exp(a1+(b1/u)))
			confidenceGap95pos1.append(math.exp(a1+(b1/u) + talpha*standardErrorOfTheEstimate1*math.sqrt(1.+(1./n)+(1/u-sum(x1Values)/n)**2/(findVariance(x1Values)*(n-1)))))
			confidenceGap95neg1.append(math.exp(a1+(b1/u) - talpha*standardErrorOfTheEstimate1*math.sqrt(1.+(1./n)+(1/u-sum(x1Values)/n)**2/(findVariance(x1Values)*(n-1)))))

if 'graph' in sys.argv:
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
	
	fig = plt.figure(figsize = (14,8))
	ax = fig.gca(projection='3d')
	
	xValues = []
	yValues = []
	zValues = []
	
	# plot the points
	for value in range(len(outputMatrix)):
		xValues.append(1/inputMatrix[value][1]) # temperature, (1/T)**-1
		yValues.append(1/inputMatrix[value][2]) # oxygen, 1/1/O
		zValues.append(math.exp(outputMatrix[value])) # timing
	ax.scatter(xValues, yValues, zValues, cmap = mpl.cm.spectral, marker = '.')

	# plot the surface	
	X = np.arange(17,28,0.25)
	Y = np.arange(9,30,0.25)
	X,Y = np.meshgrid(X,Y)
	# print surfaceMatrix[2]
	##Z1 = math.exp(surfaceMatrixExpanded[0])*(1/Y)*np.exp(surfaceMatrixExpanded[1]/X)
	Z1 = surfaceMatrixExpanded[0] + surfaceMatrixExpanded[1]*np.exp(1/X) + surfaceMatrixExpanded[2]*np.exp(1/Y)# + (1/Y)*surfaceMatrixExpanded[3]*np.exp(1/X)
	Z2 = ((65.407/Y)+0.9972)*np.exp((0.4729*Y+31.)/X) # + (1/Y)*surfaceMatrixExpanded[3]*np.exp(1/X)
#	Z = math.exp(surfaceMatrix[0])*np.exp((surfaceMatrix[1]/X)+(surfaceMatrix[2]/Y))
	#surf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,linewidth=0,antialiased=False, cmap = mpl.cm.spectral)
	#fig.colorbar(surf, shrink=0.5, aspect=5)
	#ax.plot_wireframe(X,Y,Z,rstride=10,cstride=10,linewidth=0.5,antialiased=False)
#	ax.plot_wireframe(X,Y,Z1,rstride=10,cstride=1,linewidth=0.75,antialiased=False)
	ax.plot_wireframe(X,Y,Z2,rstride=10,cstride=1,linewidth=0.75,antialiased=False)
	ax.contourf(X,Y,Z2,zdir='y',offset=9, cmap = mpl.cm.cool)
	ax.contourf(X,Y,Z2,zdir='x',offset=17, cmap = mpl.cm.autumn)
	#ax.invert_yaxis()
	ax.set_xlabel('Temperature (C)', color = 'k')
	ax.set_ylabel('Oxygen (%)', color = 'k')
	ax.set_zlabel('Developmental Time', color = 'k')
	plt.savefig('curveFitting/curveFitting3D_'+sys.argv[1]+'_combinedOxygenAndTemp_rev.pdf', transparent = True)
	plt.close()
	
	
# print '& $t_{',temperature,'} =','%.2f' % math.exp(a1),'e^{','%.2f' % b1,' / T}$ & ','%.3f' %Rsquared_1,'&$t_{',temperature,'} \plusminus', '%.3f' % (talpha*standardErrorOfTheEstimate/math.sqrt((findVariance(x1Values)*(n-1)))), '\sqrt{','%.2f' % (1.+(1./n)*(findVariance(x1Values)*(n-1))),'+ (\frac{1}{T}-','%.2f' % (sum(x1Values)/n),')^{2}}$'

