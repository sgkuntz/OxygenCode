# Curve-fitting: 
# 1. use least-squares approach to find curves for each animal across temperatures
# 2. extend this approach to different events?
# 3. determine if there are trends in the normalized data (see if they are stat sig different from vertical) -> can try F statistic for these
# in v4, try generating curve fitting of the form t = a+ c*exp(a-bT) and t = c*exp(a/kT)
# options: 'graph' and 'compare'

# python CurveFitting_v5.py Me R303 17.5 22.5 27.5 graph

import math, os, re, string, sys
import matplotlib as mpl
import pylab as pyl
import numpy as np
import scipy.stats as st

def SxySqrd(coordinates):
	n = len(coordinates)
	b,a = findSlopeAndIntercept(coordinates)
	runningSum = 0
	for (x,y) in coordinates:
		runningSum += (y-(b*x+a))**2
	return runningSum/(n-2)

def findSlopeAndIntercept(coordinates):
	n = len(coordinates)
	vectorX = []
	vectorY = []
	runningProduct = 0
	runningSquare = 0
	for (x,y) in coordinates:
		vectorX.append(x)
		vectorY.append(y)
		runningProduct += x*y
		runningSquare += x**2
	b = (n*runningProduct - sum(vectorX) * sum(vectorY))/(n*runningSquare-(sum(vectorX))**2) # slope or m
	a = (sum(vectorY) * runningSquare - sum(vectorX) * runningProduct)/(n*runningSquare-(sum(vectorX))**2) # intercept or b
	return b, a
	
def findVariance(vectorA):
	n = len(vectorA)
	runningSum = 0
	for i in vectorA:
		runningSum += i**2
	return (runningSum-(sum(vectorA)**2/n))/(n-1)

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
	
tempNames = [['12.5','12.5'],['15','15.0'],['17.5','17.5'],['20','20.0'],['22.5','22.5'],['25','25.0'],['27.5','27.5'],['30','30.0'],['32.5','32.5'],['17','17.0']]
temperatureListDefault = ['12.5','15','17.5','17','20','22.5','25','27.5','30','32.5','35']
	
	
SpList = []; Tlist = []
for species in speciesOrdering:
	for i in sys.argv:
		if changeName(i) == species:
			SpList.append(changeName(i))
for i in sys.argv:
	if i in temperatureListDefault:
		Tlist.append(changeTemp(i))
		
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
					if speciesFull not in comboDict:
						comboDict[speciesFull] = {}
					if temperature not in comboDict[speciesFull]:
						comboDict[speciesFull][temperature] = []
					comboDict[speciesFull][temperature].append((float(scaling)*(int(tfTime)-int(mryTime)))/60) # record time in hours
	else:
		print strain, temp


popList = []
if 'merge' in sys.argv:
	for species in SpList:
		if len(species)>4:
			popList.append(species)
	for species in popList:
		SpList.pop(SpList.index(species))



## Curve-fitting
# set up tValues (log transform the temperature)

print comboDict.keys()

outputFile = open('curveFitting/Curvefitting_outputData_'+sys.argv[1]+'.txt','w')

fig = pyl.figure(figsize=(16,20))
counter = 0
for species in SpList:
	counter +=1
	tValues, yValues, xValues = [], [], []
	t1Values, y1Values, x1Values = [], [], []
	t2Values, y2Values, x2Values = [], [], []
	xOriginals = []
	for T in Tlist:
		for t in comboDict[species][T]:
			tValues.append((1/(float(T)),t))
			xValues.append(1/(float(T)))
			yValues.append(t)
			xOriginals.append(float(T))
	b,a = findSlopeAndIntercept(tValues)	
	n = len(tValues)
	Rsquared = 1 - ((n-2)*SxySqrd(tValues)/((n-1)*findVariance(yValues))) # Pearson Product-Moment Correlation Coefficient

	for T in Tlist:
		for t in comboDict[species][T]:
			t1Values.append((1/float(T),math.log(float(t))))
			x1Values.append(1/float(T))
			y1Values.append(math.log(float(t)))
	b1,a1 = findSlopeAndIntercept(t1Values)	
	n = len(t1Values)
	Rsquared_1 = 1 - ((n-2)*SxySqrd(t1Values)/((n-1)*findVariance(y1Values))) # Pearson Product-Moment Correlation Coefficient
	
	print Rsquared, Rsquared_1

	# z = np.arange(17, 28.5,.1)
	z = np.arange(15, 33.5,.1)
	w, w1, w2,w3, w4 = [], [], [], [], []
	confidenceGap50pos, confidenceGap95pos, confidenceGap50neg, confidenceGap95neg = [], [], [], []
	confidenceGap50pos1, confidenceGap95pos1, confidenceGap50neg1, confidenceGap95neg1 = [], [], [], []
	talpha = st.t.isf(0.05/2,(n-2)) # isf is one-tailed, so divide by 2 for 2-tailed critical value
	tbeta = st.t.isf(0.5/2,(n-2))
	standardErrorOfTheEstimate = math.sqrt((findVariance(yValues)-findVariance(xValues)*b**2)*(n-1)/(n-2))
	standardErrorOfTheEstimate1 = math.sqrt((findVariance(y1Values)-findVariance(x1Values)*b1**2)*(n-1)/(n-2))
	for u in z:
		#print species[-2:]
		w.append(b*(1/u)+a)
		w1.append(math.exp(a1+(b1/u))) # red
		w2.append(-1219.68843958+1182.25510561*math.exp(1/u)+ 228.55053286/float(species[-2:])) # green
		w3.append(((65.41/float(species[-2:]))+1.)*math.exp((0.47*float(species[-2:])+31.)/u)) # blue
	#	w4.append(281.5*math.exp(-0.126*u)+(3*(u*u)-142.84*u+1861)/21.) #[  -38.72701231  1238.8803309    228.6350275 ]
		
		#confidenceGap95pos.append(b*math.log(u-15)+a + talpha*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(math.log(u-15)-sum(xValues)/n)**2/(findVariance(xValues)*(n-1))))
		#confidenceGap95neg.append(b*math.log(u-15)+a - talpha*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(u-sum(xOriginals)/n)**2/(findVariance(xOriginals)*(n-1))))

		confidenceGap95pos1.append(math.exp(a1+(b1/u) + talpha*standardErrorOfTheEstimate1*math.sqrt(1.+(1./n)+(1/u-sum(x1Values)/n)**2/(findVariance(x1Values)*(n-1)))))
		confidenceGap95neg1.append(math.exp(a1+(b1/u) - talpha*standardErrorOfTheEstimate1*math.sqrt(1.+(1./n)+(1/u-sum(x1Values)/n)**2/(findVariance(x1Values)*(n-1)))))


	if 'graph' in sys.argv:
		pyl.subplot(6,1,counter)
		pyl.plot(yValues, xOriginals, marker = 's', markersize = 5, c = 'cyan', linestyle = '')
		pyl.title(species, fontsize = 30)
		pyl.ylabel('Temperature ('+unichr(176)+'C)', fontsize = 20)
		pyl.xlabel('time (hours post-cellularization)', fontsize = 20)
#		pyl.xlim(10,55)
#		pyl.ylim(17,28)
		pyl.xlim(0,80)
		pyl.ylim(16,30)
		# find confidence interval for future estimates (lower line and upper line)
	#	pyl.plot(w, z, marker = '',linestyle = '-', lw = 1.5, c = 'red')
		pyl.plot(w1, z, marker = '',linestyle = '-', lw = 1.5, c = 'red')
		pyl.plot(w2, z, marker = '',linestyle = '-', lw = 1.5, c = 'green')
		pyl.plot(w3, z, marker = '',linestyle = '--', lw = 1.5, c = 'blue')
	#	pyl.plot(w4, z, marker = '',linestyle = '--', lw = 1.5, c = 'black')
		
		pyl.plot(confidenceGap95pos1, z, marker = '', lw = 0.5, linestyle = '--', c = 'orange')
		pyl.plot(confidenceGap95neg1, z, marker = '', lw = 0.5, linestyle = '--', c = 'orange')
		meansDict = {}
		outputFile.write(species)
		for i in range(len(yValues)):
			if xOriginals[i] not in meansDict:
				meansDict[xOriginals[i]] = []
			meansDict[xOriginals[i]].append(yValues[i])
		for T in Tlist:
			# print len(meansDict[float(T)])
			outputFile.write('\t' + str(np.mean(meansDict[float(T)])))
		outputFile.write('\n')

	print '$t_{',species,'} =', '%.2f' % a,'%.2f' % b,' \ln({T-15})$ & ','%.3f' %Rsquared,'&$t_{',species,'} \plusminus', '%.3f' % (talpha*standardErrorOfTheEstimate/math.sqrt((findVariance(xValues)*(n-1)))), '\sqrt{','%.2f' % (1.+(1./n)*(findVariance(xValues)*(n-1))),'+ (\ln(T-15)-','%.2f' % (sum(xValues)/n),')^{2}}$'
	print '& $t_{',species,'} =','%.2f' % math.exp(a1),'e^{','%.2f' % b1,' / T}$ & ','%.3f' %Rsquared_1,'&$t_{',species,'} \plusminus', '%.3f' % (talpha*standardErrorOfTheEstimate/math.sqrt((findVariance(x1Values)*(n-1)))), '\sqrt{','%.2f' % (1.+(1./n)*(findVariance(x1Values)*(n-1))),'+ (\frac{1}{T}-','%.2f' % (sum(x1Values)/n),')^{2}}$'

if 'graph' in sys.argv:
	pyl.savefig('curveFitting/curveFitting_'+sys.argv[1]+'_logTemp_TwoFits.pdf', transparent = True)
	pyl.close()

## compare all vs. all (statistical backing to independence between climate groups)

# go through each species and generate the relevant statistics
if 'compare' in sys.argv:
	for species1 in SpList:
		for species2 in SpList:
			tValues1, tValues2, tValuesAll = [], [], []
			for T in Tlist:
				for t in comboDict[species1][T]:
					tValues1.append((math.log((float(T)-15)),t))
					tValuesAll.append((math.log((float(T)-15)),t))
				for t in comboDict[species2][T]:
					tValues2.append((math.log((float(T)-15)),t))
					tValuesAll.append((math.log((float(T)-15)),t))
			n1 = len(tValues1)
			n2 = len(tValues2)
			var_pooled = ((n1-2)*SxySqrd(tValues1) + (n2-2)*SxySqrd(tValues2))/(n1+n2-4)
			var_common = SxySqrd(tValuesAll)
			var_improvement = ((n1 + n2 - 2)*var_common - (n1 + n2 -4)* var_pooled)/2
			dof_common = 2
			dof_improvement = n1 + n2 -4
			F = var_improvement/var_pooled
			F_critical = st.f.isf(0.05,dof_common, dof_improvement)
			if F > F_critical:
				# print species1, species2,'F=',F
				pass
			elif species1 == species2:
				pass
			else:
				print species1, species2, 'similar'
