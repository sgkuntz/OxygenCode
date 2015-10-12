# example command: python devEventAdjDataProcessing_vLinePlot5.py graphX Pe Me 27.5 30 white
# must be run from Dropbox folder
# analyzes only manually verified data, including mry, tf, and GUI-verified data

import os, sys, re, string, math, array, random
import scipy.stats as st
import pylab as pyl
import matplotlib as mpl
import numpy as np

### Part 0 - Can read in info from the TimeLapseInstructions file and the TimeLapseCorrections files
global backgroundSetting
if 'white' in sys.argv:
	backgroundSetting = 'white'
	markerShape = 's'
	lineColor = 'k'
	faceColor = 'w'
elif 'black' in sys.argv:
	backgroundSetting = 'black'
	markerShape = 'o'
	lineColor = 'w'
	faceColor = 'k'


speciesNames = [['Mj', 'Dmoj'],['An', 'Dana'],['Er','Dere'],['Me', 'Dmel'],['Me18','Dm18'],['Ps','Dpse'],['Pe','Dper'],['Si','Dsim'],['Se','Dsec'],['Wi','Dwil'],['Vi','Dvir'],['Ya','Dyak'],['Psch','DpseCH'],['Mjch','DmojCH'],['Mech','DmelCH'],['Mehc','DmelHC'],['10%O2','10%Oxygen'],['14%O2','14%Oxygen'],['17%O2','17%Oxygen'],['21%O2','21%Oxygen'],['25%O2','25%Oxygen'],['29%O2','29%Oxygen'],['10%O2Ps','10%OxygenDpse'],['14%O2Ps','14%OxygenDpse'],['17%O2Ps','17%OxygenDpse'],['21%O2Ps','21%OxygenDpse'],['25%O2Ps','25%OxygenDpse'],['29%O2Ps','29%OxygenDpse'],['18maleX25female','Dmel18m25f@25'],['Thor','Thor'],['chrb','chrb'],['CG16758','CG16758']]#,['18femaleX25male','21%Oxygen']]
tempNames = [['12.5','12.5'],['15','15.0'],['17.5','17.5'],['20','20.0'],['22.5','22.5'],['25','25.0'],['27.5','27.5'],['30','30.0'],['32.5','32.5'],['17','17.0'],['12.7','12.7'],['21.7','21.7']]

RaleighStrains = ['R303','R329','R379','R380','R437','R705']
HeterochronicStrains = ['12414','12592','12622','16532']
melStrains = RaleighStrains + HeterochronicStrains
CYyakStrains = ['CY01A','CY04B','CY08A','CY15B4','CY17C','CY22B','CY23A','CY21B','CY27A9']
NYyakStrains = ['NY20','NY48','NY55','NY56','NY60','NY62','NY65','NY73','NY81','NY85']
yakStrains = CYyakStrains + NYyakStrains 
NSsimStrains = ['NS37','NS67','NS79','NS113','NS137']
MDsimstrains = ['MD15','MD63','MD73','MD105','MD106','MD199','MD221','MD233','MD251']
simStrains = NSsimStrains + MDsimstrains
otherSpecies = ['An','Er','Me','MeAcc','Me18','Mj','Pe','Ps','Se','Si','Vi','Wi','Ya','10%O2','14%O2','17%O2','21%O2','25%O2','29%O2','10%O2Ps','14%O2Ps','17%O2Ps','21%O2Ps','25%O2Ps','29%O2Ps','18maleX25female','18femaleX25male','Thor','chrb','CG16758']
extraExperiments = ['Psch','Mjch','Mech','Mehc']
speciesListDefault = melStrains + yakStrains + simStrains + otherSpecies + extraExperiments
temperatureListDefault = ['12.5','15','17.5','17','20','22.5','25','27.5','30','32.5','35','12.7','21.7']
temperatureListReference = ['12.5','15.0','17.5','17.0','20.0','22.5','25.0','27.5','30.0','32.5','35','12.7','21.7']
possibleEventList = ['posterior_gap','pole_bud_appears','nuclei_at_periphery','pole_cells_form','yolk_contraction','cellularization_begins','membrane_reaches_yolk','pole_cells_migrate','cephalic_furrow_formation','pole_cell_invagination','cephalic_furrow_reclines','amnioproctodeal_invagination','amnioproctodeal_invaginationA','amnioproctodeal_invagination_begins','amnioproctodeal_invagination_ends','transversal_fold_formation','anterior-midgut_primordium','stomodeal_plate_forms','stomodeum_invagination','clypeolabral_lobe_forms','clypeolabrum_rotates','posterior_gap_before_germband_shortens','germband_retraction_starts','amnioserosa','germband_retracted','dorsal_divot','clypeolabrum_retracts','anal_plate_forms','midgut_unified','clypeolabrum_ventral_lobe_even','gnathal_lobes_pinch','head_involution_done','heart-shaped_midgut','convoluted_gut','muscle_contractions','trachea_fill','hatch', 'temperature_shift']

speciesOrdering = ['Dsim', 'Dsec','Dmel','Dere','Dyak','Dana','Dper','Dpse','Dwil','Dmoj','Dvir', 'DpseCH', 'DmojCH', 'DmelCH', 'DmelHC','10%Oxygen','14%Oxygen','17%Oxygen','Dmel18m25f@25','21%Oxygen','25%Oxygen','29%Oxygen','10%OxygenDpse','14%OxygenDpse','17%OxygenDpse','21%OxygenDpse','25%OxygenDpse','29%OxygenDpse','Thor','chrb','CG16758']

SpList = []; Tlist = []

def changeName(i):
	for (shortName, fullName) in speciesNames:
		if i == shortName:
			return fullName
	newName = ''
	if '-merge' not in sys.argv:
		newName = i
	if i in melStrains:
		return 'Dmel'+newName
	elif i in yakStrains:
		return 'Dyak'+newName
	elif i in simStrains:
		return 'Dsim'+newName
def changeTemp(j):
	for (temp, tempInt) in tempNames:
		if j == temp:
			return tempInt

for i in sys.argv:
	if i in speciesListDefault:
		SpList.append(changeName(i))
	elif i in temperatureListDefault:
		Tlist.append(changeTemp(i))
		
if SpList == []:
	for species in speciesListDefault:
		SpList.append(changeName(species))
if Tlist == []:
	for temp in temperatureListDefault:
		Tlist.append(changeTemp(temp))
print SpList

comboDict = {} ; animalDict = {}; eventListPrelim = ['membrane_reaches_yolk','trachea_fill']

for (dirpath, dirnames, filenames) in os.walk("."):
	for filename in filenames:
		if 'TimeLapseCorrections_' in filename:
			if 'deleted' not in filename:
				try:
					[title, trial, species, temperature] = re.split("_",filename)
				except:
					species = 'fail'
				if changeName(species) in SpList:
					if temperature[:4] in Tlist:
						correctionsFile = open(dirpath+'/'+filename,'r')
						for i in correctionsFile:
							line = re.split("[\s]",i)
							[condition, event, datePosition, timepoint, dilation, filler] = re.split("[\s]",i)
							if condition <> 'conditions':
								if timepoint <> 'NaN':
									if event not in eventListPrelim:
										eventListPrelim.append(event)
									if '-merge' in sys.argv:
										if condition[4:] in yakStrains:
											condition = condition[:4] + 'Ya'
										elif condition[4:] in melStrains:
											condition = condition[:4] + 'Me'
										elif condition[4:] in simStrains:
											condition = condition[:4] + 'Si'
									if condition[:4] in Tlist:
										condition = condition[:4] + changeName(species)
									elif condition[:2] in temperatureListDefault:
										condition = changeTemp(condition[:2]) + changeName(species)
									if condition not in comboDict:
										comboDict[condition] = {}
									if datePosition not in comboDict[condition]:
										comboDict[condition][datePosition] = {}
									if event not in comboDict[condition][datePosition]:
										comboDict[condition][datePosition][event] = [changeName(species), int(dilation), timepoint]
									if event not in eventListPrelim:
										eventListPrelim.append(event)
		
timeLapseInfoFile = open('TimeLapseInstructions_run_trial15.txt','r')
for i in timeLapseInfoFile:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temp, scaling, zStack, filler] = re.split("[\s]",i)
	if changeName(strain) in SpList:
		speciesFull = changeName(strain) # change the species name
		if changeTemp(temp) in Tlist:
			temperature = changeTemp(temp)
			if temperature in temperatureListReference:
				if temperature+speciesFull not in comboDict:
					comboDict[temperature+speciesFull] = {}
				if date+'_'+position not in comboDict[temperature+speciesFull]:
					comboDict[temperature+speciesFull][date+'_'+position] = {}
				if 'membrane_reaches_yolk' not in comboDict[temperature+speciesFull][date+'_'+position]:
					comboDict[temperature+speciesFull][date+'_'+position]['membrane_reaches_yolk'] = [speciesFull, int(scaling), mryTime]
				if 'trachea_fill' not in comboDict[temperature+speciesFull][date+'_'+position]:
					comboDict[temperature+speciesFull][date+'_'+position]['trachea_fill'] = [speciesFull, int(scaling), tfTime]
	
timeLapseInfoFile2 = open('TimeLapseInstructions_supplemental_trial15.txt','r')
for i in timeLapseInfoFile2:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temp, scaling, zStack, filler] = re.split("[\s]",i)
	if changeName(strain) in SpList:
		speciesFull = changeName(strain) # change the species name
		if changeTemp(temp) in Tlist:
			temperature = changeTemp(temp)
			if temperature in temperatureListReference:
				if temperature+speciesFull not in comboDict:
					comboDict[temperature+speciesFull] = {}
				if date+'_'+position not in comboDict[temperature+speciesFull]:
					comboDict[temperature+speciesFull][date+'_'+position] = {}
				if 'membrane_reaches_yolk' not in comboDict[temperature+speciesFull][date+'_'+position]:
					comboDict[temperature+speciesFull][date+'_'+position]['membrane_reaches_yolk'] = [speciesFull, int(scaling), mryTime]
				if 'trachea_fill' not in comboDict[temperature+speciesFull][date+'_'+position]:
					comboDict[temperature+speciesFull][date+'_'+position]['trachea_fill'] = [speciesFull, int(scaling), tfTime]


comboListPrime = comboDict.keys()
comboListPrime.sort()
comboList = []
for species in speciesOrdering:
	for combo in comboListPrime:
		if species == combo[4:]:
			comboList.append(combo)

eventList = []			
for event in possibleEventList:
	if event in eventListPrelim:
		eventList.append(event)

print comboList
			
meanDict = {}

numberOfTemperatures = len(Tlist) # number of Temperatures
numberOfSpecies = len(SpList) # number of Species

### Part 0b - Dot Line Graph (event determines color, condition determines y-value, time determines x-value
fig = pyl.figure(figsize=(18,12))#8, 10)) # [8, 14)) for poster] # 15,10))  8, 14)) [ 12, 8)) for figures] [15, 15)) for other figures]
ax = fig.add_subplot(111)
yCounter = 0

revisedEventList = ['pole_bud_appears','membrane_reaches_yolk','pole_cell_invagination','amnioproctodeal_invaginationA','amnioproctodeal_invagination_begins','amnioproctodeal_invagination_ends','amnioserosa','clypeolabrum_retracts','clypeolabrum_ventral_lobe_even','heart-shaped_midgut','trachea_fill', 'temperature_shift']
#revisedEventList = ['pole_bud_appears','membrane_reaches_yolk','pole_cell_invagination','amnioproctodeal_invagination_begins','amnioserosa','clypeolabrum_retracts','clypeolabrum_ventral_lobe_even','heart-shaped_midgut','trachea_fill', 'temperature_shift']
#revisedEventList = ['pole_bud_appears','membrane_reaches_yolk','pole_cell_invagination','amnioproctodeal_invagination','amnioserosa','clypeolabrum_retracts','clypeolabrum_ventral_lobe_even','heart-shaped_midgut','trachea_fill', 'temperature_shift']

colorDict = {}
colorCounter = 0
for event in revisedEventList:
	colorCounter += 1.5
	if colorCounter > 9.5:
		colorCounter = 0.5
	elif colorCounter == 7.5:
		colorCounter = 8
	if colorCounter <> 6:
		colorDict[event] = colorCounter/10.
	else:
		colorDict[event] = 7/10.

conditionInterval = 150
legendLabelList = []
averageEventValue = {}

for conditions in comboList: # cycle through each T+Sp combo
	normalizedTiming_allEvents, normalizedTimingRev_allEvents = [], []
	numberOfAnimals_allEvents = []
	animalListTemp = comboDict[conditions].keys()
	animalListTemp.sort()
	yCounter += conditionInterval
	posCounter = 0
	normalizedTiming = []
	for animal in animalListTemp:
		normalizedTiming = []
		posCounter += 3
		eventListTemp = comboDict[conditions][animal].keys()
		eventListTemp.sort()
		for event in eventList:
			if event <> 'germband_maxima':
				try:
					normalizedTiming = (comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal][event][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2])))
					if conditions+event not in averageEventValue:
						averageEventValue[conditions+event] = []
					averageEventValue[conditions+event].append(normalizedTiming)
					yCoor = (int(posCounter/2 + yCounter-3*len(animalListTemp)/4))
					if 'white' in sys.argv:
						ax.plot(normalizedTiming, yCoor, marker = markerShape, label = event, linestyle = '', markersize = 2, c=mpl.cm.spectral(colorDict[event]), markeredgecolor=mpl.cm.spectral(colorDict[event]))						
					elif 'black' in sys.argv:
						ax.plot(normalizedTiming, yCoor, marker = markerShape, label = event, linestyle = '', markersize = 4, c=mpl.cm.hsv(colorDict[event]), markeredgecolor=mpl.cm.hsv(colorDict[event]))
				except:
					pass
indexVector = []
for x in range(len(comboList)+2):
	indexVector.append(x*conditionInterval)
if len(Tlist) > 1:
	labelList = ['']
	for combo in comboList:
		labelList.append((combo[:4]+ unichr(176)+'C'))
	labelList.append('')
	figTitle = SpList[0]
if len(SpList) > 1:
	labelList = ['']
	for combo in comboList:
		labelList.append(combo[4:])
	labelList.append('')
	figTitle = str(Tlist[0]) + unichr(176)+'C'
index = np.array(indexVector)
handles, labels = ax.get_legend_handles_labels()
legendHLList = zip(handles, labels)
revisedLegendHLList = []
for event in eventList:
	for (handle, label) in legendHLList:
		if label == event:
			if event not in legendLabelList:
				legendLabelList.append(event)
				revisedLegendHLList.append((handle,label))
				yCoor = 0 
				for conditions in comboList:
					tempValueList = []
					yCoor += conditionInterval
					for item in averageEventValue[conditions+event]:
						if math.isnan(item) == False:
							tempValueList.append(item)
					ax.plot(float(sum(tempValueList))/len(tempValueList), yCoor, marker = 'd', label = event, linestyle = '', markersize = 14, c=mpl.cm.spectral(colorDict[event]), markeredgecolor=mpl.cm.spectral(colorDict[event]))
newList = zip(*revisedLegendHLList)
ax.legend(*newList)
pyl.xticks(color = lineColor, fontsize = 20)
pyl.yticks(index, labelList, color = lineColor, fontsize = 25)
pyl.xlabel('time (min post-cellularization)', color = lineColor, fontsize = 25)
pyl.title(figTitle, color = lineColor, fontsize = 30)
pyl.xlim(-500,2500)
ax.spines['bottom'].set_color(lineColor)
ax.spines['left'].set_color(lineColor)
ax.yaxis.label.set_color(lineColor)
ax.tick_params(axis='x', colors=lineColor)
ax.tick_params(axis='y', colors=lineColor)
pyl.savefig('dotLinePlots/dotLinePlot_' + str(sys.argv[1]) + '_' + figTitle + '_' + backgroundSetting +  '.pdf', transparent = True, facecolor = faceColor, linecolor = 'w') # to change from black to white background, make non-transparent
pyl.close()

# Normalized graph
textOutput = open('dotLinePlots/dotLineNumbers_' + str(sys.argv[1]) + '_' + figTitle + '.txt', 'w')
#diagnosticOutput = open('dotLinePlots/dotLineNumbers_'+ str(sys.argv[1])+'_'+figTitle + '_diagnostic.txt','w')

fig = pyl.figure(figsize=(11,14)) # [14, 6)) for poster]
ax = fig.add_subplot(111)
yCounter = 0
statsDict = {}
statsDictSpecial = {}
legendLabelList = []
for conditions in comboList: # cycle through each T+Sp combo
	conditionDict = {}
	normalizedTiming_allEvents, normalizedTimingRev_allEvents = [], []
	numberOfAnimals_allEvents = []
	animalListTemp = comboDict[conditions].keys()
	animalListTemp.sort()
	yCounter += conditionInterval
	posCounter = 0
	normalizedTiming = []
	for animal in animalListTemp:
		normalizedTiming = []
		posCounter += 3
		eventListTemp = comboDict[conditions][animal].keys()
		eventListTemp.sort()
		for event in eventList:
			try:
				normalizedTiming = (comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal][event][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2])))
				yCoor = (int(posCounter + yCounter-3*len(animalListTemp)/2))
				scaledNormalizedTiming = (normalizedTiming/(comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal]['trachea_fill'][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2]))))
				if 'white' in sys.argv:
					ax.plot(scaledNormalizedTiming, yCoor, marker = markerShape, label = event, linestyle = '', markersize = 5, c=mpl.cm.spectral(colorDict[event]), markeredgecolor=mpl.cm.spectral(colorDict[event]))
				elif 'black' in sys.argv:
					ax.plot(scaledNormalizedTiming, yCoor, marker = markerShape, label = event, linestyle = '', markersize = 5, c=mpl.cm.hsv(colorDict[event]), markeredgecolor=mpl.cm.hsv(colorDict[event]))
				#ax.plot(scaledNormalizedTiming, yCoor, marker = markerShape, label = event, linestyle = '', markersize = 5, c=mpl.cm.gist_rainbow(colorDict[event]), markeredgecolor=mpl.cm.gist_rainbow(colorDict[event]))
				if event not in conditionDict:
					conditionDict[event] = [0,0, []]
				conditionDict[event][0] += scaledNormalizedTiming
				conditionDict[event][1] += 1
				conditionDict[event][2].append(scaledNormalizedTiming)
				if event not in statsDict: # write to a dictionary of vectors for statistical analysis (one key for each event)
					statsDict[event] = []
					statsDictSpecial[event] = []
				statsDict[event].append((float(conditions[:4]), scaledNormalizedTiming)) # vector needs two dimensions: time and temperature/species -> record conditions
				statsDictSpecial[event].append((float(conditions[4:6]), scaledNormalizedTiming)) # vector needs two dimensions: time and temperature/species -> record conditions
				diagnosticOutput.write(str(event)+'\t'+str(conditions[4:6]) + '\t'+ str(scaledNormalizedTiming) + '\n')
			except:
				pass
	for event in eventList:
		try:
			meanValue = conditionDict[event][0]/conditionDict[event][1]
			varianceValue = []
			for scaledValue in conditionDict[event][2]:
				varianceValue.append((scaledValue-meanValue)**2)
			stDevValue = math.sqrt(sum(varianceValue)/(len(varianceValue)-1))
			textOutput.write(str(conditions[:4]) + '\t' + str(conditions[4:]) + '\t' + str(event) + '\t' + str(meanValue) + '\t' + str(stDevValue) + '\n') 
		except:
			pass
index = np.array(indexVector)
handles, labels = ax.get_legend_handles_labels()
legendHLList = zip(handles, labels)
revisedLegendHLList = []
for event in eventList:
	for (handle, label) in legendHLList:
		if label == event:
			if event not in legendLabelList:
				legendLabelList.append(event)
				revisedLegendHLList.append((handle,label))
newList = zip(*revisedLegendHLList)
pyl.yticks(index, labelList, color = lineColor, fontsize = 20)
pyl.xticks(color = lineColor, fontsize = 20)
pyl.title(figTitle, color = lineColor, fontsize = 25)
pyl.xlim(-0.2,1.1)
ax.spines['bottom'].set_color(lineColor)
ax.spines['left'].set_color(lineColor)
ax.yaxis.label.set_color(lineColor)
ax.tick_params(axis='x', colors=lineColor)
ax.tick_params(axis='y', colors=lineColor)
pyl.xlabel('proportion of development', color = lineColor, fontsize = 20)
pyl.savefig('dotLinePlots/dotLinePlot_' + str(sys.argv[1]) + '_' + figTitle + '_' + backgroundSetting +  '_Normalized.pdf', transparent = True, facecolor = faceColor, linecolor = 'w') # to change from black to white background, make non-transparent
pyl.close()



statsList = []
statsDictList = statsDict.keys()
for stage in revisedEventList:
	if stage in statsDictList:
		statsList.append(stage)
exciseList = ['membrane_reaches_yolk', 'trachea_fill']#, 'amnioproctodeal_invaginationA','amnioproctodeal_invagination_ends'] # 'germband_maxima', 'clypeolabrum_ventral_lobe_even']
for event in exciseList:
	try:
		statsList.pop(statsList.index(event))
	except:
		pass


if 'distribution' in sys.argv:
	# make a distribution of the occurrences, x axis is time, y axis is number of occurrences
	fig = pyl.figure(figsize=(25,6))
	ax = fig.add_subplot(121)

	for event in statsList:
		xValues = []
		for (temperature, time) in statsDict[event]:
			xValues.append(time)
		n = len(xValues)/4
		pyl.hist(xValues, n, label = event, edgecolor = mpl.cm.spectral(colorDict[event]), color = mpl.cm.spectral(colorDict[event])) #, histtype = 'stepfilled')
	pyl.xlim(-.2, 1.)
	handles, labels = ax.get_legend_handles_labels()
	legendHLList = zip(handles, labels)
	revisedLegendHLList = []
	for event in eventList:
		for (handle, label) in legendHLList:
			if label == event:
				if event not in legendLabelList:
					legendLabelList.append(event)
					revisedLegendHLList.append((handle,label))
	newList = zip(*revisedLegendHLList)
	# print zip(*revisedLegendHLList)
	pyl.title('Distribution across all conditions'+sys.argv[1], color = lineColor, fontsize = 25)
	ax.legend(*newList, bbox_to_anchor=(1.05,1), loc = 2, borderaxespad = 0.)
	pyl.xlabel('proportion of development', color = lineColor, fontsize = 20)
	pyl.xticks(fontsize = 20)
	pyl.yticks(fontsize = 20)
	pyl.ylabel('number of occurrences', color = lineColor, fontsize = 20)
	pyl.savefig('curveFitting/eventDistribution_'+sys.argv[1]+'.pdf',transparent = True)
	pyl.close()



# Statistical analysis
if 'stats' in sys.argv:

	def dictionarySelection(DummyEvent):
		if len(Tlist) == 1:
			return statsDictSpecial[DummyEvent]
		else:
			return statsDict[DummyEvent]

	def SyxSqrd(coordinates):
		vectorX = []
		b,a = findSlopeAndIntercept(coordinates)
		runningSum = 0
		for (x,y) in coordinates:
			if math.isnan(y) == False:
				runningSum += (y-(b*x+a))**2
				vectorX.append(x)
		return runningSum/(len(vectorX)-2)

	def findSlopeAndIntercept(coordinates):
		n = len(coordinates)
		vectorX = []
		vectorY = []
		runningProduct = 0
		runningSquare = 0
		for (x,y) in coordinates:
			if math.isnan(y) == False:
				vectorX.append(x)
				vectorY.append(y)
				runningProduct += x*y
				runningSquare += x**2
		b = ((len(vectorX)*runningProduct) - (sum(vectorX) * sum(vectorY)))/(len(vectorX)*runningSquare-(sum(vectorX))**2)
		a = (sum(vectorY) * runningSquare - sum(vectorX) * runningProduct)/(len(vectorX)*runningSquare-(sum(vectorX))**2)
		return b, a
	
	def findVariance(vectorA):
		runningSum = 0
		vectorB = []
		for (i,j) in vectorA:
			if math.isnan(j) == False:
				runningSum += i**2
				vectorB.append(i)
		return (runningSum-(sum(vectorB)**2/len(vectorB)))/(len(vectorB)-1)

	numberOfComparisons = 0

	statsOutputFile = open('curveFitting/curveFittingStats_'+sys.argv[1]+'_'+'eventTrends.txt','w')
	statsOutputFile.write('Event1' + '\t'+ 'Slope1' + '\t'+ 'Intercept1' + '\t'+ 'Event2' + '\t'+ 'Slope2' + '\t'+ 'Intercept2' + '\t'+ 'tValue' + '\t'+ 'DOF' + '\t'+ 'criticalValue' + '\t'+ 'pValue' + '\n')

	maxPValue = 0.0
	heatMapList = []
	redundancyList = []
	for event1 in statsList:
		for event2 in statsList:
			if event1 <> event2:
				numberOfComparisons +=1
	for event1 in statsList:
		redundancyList.append(event1)
		for (timing, temperature) in dictionarySelection(event1):
			pyl.figure(0)
			pyl.plot(timing, temperature, marker = markerShape, linestyle = '', markersize = 3, c=mpl.cm.spectral(colorDict[event1]))
			pyl.figure(1,figsize = (8,5)) # good for 3 temperature graphs across 5 oxygen concentrations each
			# pyl.figure(1,figsize = (10,2)) # good for 5 oxygen graphs across 3 temperatures each
			pyl.plot(temperature, timing-0.25+(random.random()/2.), marker = markerShape, linestyle = '', markersize = 3, c=mpl.cm.spectral(colorDict[event1]),markeredgecolor=mpl.cm.spectral(colorDict[event1]))
		for event2 in statsList:
			if event1 <> event2 and event2 not in redundancyList:
				slope1, intercept1 = findSlopeAndIntercept(dictionarySelection(event1))
				slope2, intercept2 = findSlopeAndIntercept(dictionarySelection(event2))
				n1 = len(dictionarySelection(event1))
				n2 = len(dictionarySelection(event2))
				# find if their slopes are significantly different (assume the intercepts will be)
				Syx1 = SyxSqrd(dictionarySelection(event1))
				Syx2 = SyxSqrd(dictionarySelection(event2))
				var1 = findVariance(dictionarySelection(event1))
				var2 = findVariance(dictionarySelection(event2))		
				var_pooled = (Syx1*(n1-2)+Syx2*(n2-2))/(n1+n2-4)
				tTest = (slope1-slope2)/(math.sqrt((var_pooled/((n1-1)*var1))+(var_pooled/((n2-1)*var2))))
				dof = n1+n2-4
				criticalValue = st.t.isf(0.001/2,dof)
				pValue = (numberOfComparisons)*(1.0 - (st.t.cdf(math.fabs(tTest),dof) - st.t.cdf(-math.fabs(tTest), dof)))
			
				statsOutputFile.write(str(event1) + '\t'+ str(slope1) + '\t'+ str(intercept1) + '\t'+ str(event2) + '\t'+ str(slope2) + '\t'+ str(intercept2) + '\t'+ str(math.fabs(tTest)) + '\t'+ str(dof) + '\t'+ str(criticalValue) + '\t'+ str(pValue) + '\n')
				# record values for heat map: event1, event2, pValue
				if pValue > maxPValue:
					maxPValue = pValue
				# heatMapList.append((event1, event2, pValue))
				if pValue > 0.05:
					pValue = 0.05
				heatMapList.append((event1, event2, pValue))
				

	heatMapArray = []
	for event1 in statsList:
		heatMapArrayTemp = []
		for event2 in statsList:
			if event1 == event2:
				# heatMapArrayTemp.append(maxPValue)
				heatMapArrayTemp.append(0.05)
			else:
				for (eventA, eventB, pValue) in heatMapList:
					if event1 == eventA and event2 == eventB:
						heatMapArrayTemp.append(pValue)
					elif event1 == eventB and event2 == eventA:
						heatMapArrayTemp.append(pValue)
						# print pValue
		heatMapArray.append(heatMapArrayTemp)
	
				
#	cdict = {'red': ((0.0, 0.0, 0.0),(0.5, 0.0, 0.0),(1.0, 1.0, 1.0)),'green': ((0.0, 0.0, 0.0),(0.5, 0.0, 0.0),(1.0, 1.0, 1.0)),'blue': ((0.0, 1.0, 1.0),(0.5, 0.0, 0.0),(1.0, 0.0, 0.0))}
#	my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	Z = np.array(heatMapArray)
	pyl.figure(figsize=(16,12))
	pyl.pcolor(Z,cmap=pyl.get_cmap("spectral"))
	pyl.colorbar()
	pyl.xlabel('event', color = 'k', fontsize = 15)
	indexx = np.arange(len(statsList))
	pyl.xticks(indexx+0.5, statsList, color = 'k', fontsize = 10, rotation = 90)
	pyl.ylabel('event', color = 'k', fontsize = 15)
	pyl.title(sys.argv[1], color = 'k', fontsize = 20)
	indexy = np.arange(len(statsList))
	pyl.yticks(indexy+0.5, statsList, color = 'k', fontsize = 10)
	pyl.savefig('curveFitting/curveFitting_'+sys.argv[1]+'_heatMap.pdf', transparent = True)
	pyl.close()

	for event in statsList:
		if len(Tlist) == 1:
			z = np.arange(9,30,2)
		else:
			z = np.arange(17,35,2)
		slope, intercept = findSlopeAndIntercept(dictionarySelection(event)) 
		# print event, slope, intercept
		pyl.figure(0)
		pyl.plot(z, z*slope+intercept, marker = '', linestyle = '--', lw = 1, c=mpl.cm.spectral(colorDict[event]))
		pyl.figure(1)
		pyl.plot(z*slope+intercept, z, marker = '', linestyle = '--', lw = 1, c=mpl.cm.spectral(colorDict[event]))
	
		n = len(statsDict[event])
		talpha = st.t.isf(0.05/2,(n-2))
		var = findVariance(dictionarySelection(event))
		runningSum = 0
		stError = []
		for (x,y) in dictionarySelection(event):
			runningSum += x
		mean = runningSum/n
		upGap = []
		downGap = []
		for u in z:
			stError.append(talpha*math.sqrt(SyxSqrd(dictionarySelection(event)))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var)))
			upGap.append((u*slope+intercept - talpha*math.sqrt(SyxSqrd(dictionarySelection(event)))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var))))
			downGap.append((u*slope+intercept + talpha*math.sqrt(SyxSqrd(dictionarySelection(event)))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var))))
		try:
			pyl.figure(0)
			pyl.fill_between(z, downGap, upGap, edgecolor = mpl.cm.spectral(colorDict[event]), facecolor = mpl.cm.spectral(colorDict[event]))
		except:
			print "error in figure 1"
	pyl.figure(0)
	pyl.plot(z,z*0., marker = '', linestyle = '--', lw = 1, c = 'k')
	pyl.figure(1)
	pyl.plot(z*0.,z, marker = '', linestyle = '--', lw = 1, c = 'k')

	# establish confidence intervals for the regression lines

	spNames = ''
	Tnames = ''
	for species in SpList:
		spNames += species+'_'
	for temperature in Tlist:
		Tnames +=temperature+'_'

	pyl.figure(0)
	if len(Tlist) == 1:
		pyl.xlabel('oxygen (%)')
		pyl.xlim(9,30)		
	else:
		pyl.xlabel('temperature (C)')
		pyl.ylim(17,28)
	pyl.ylabel('proportion of development')
	pyl.ylim(-0.2,1.)
	pyl.title(sys.argv[1])
	pyl.savefig('curveFitting/curveFitting_'+sys.argv[1]+'_eventTrends_inv.pdf', transparent = True)
	pyl.close()
	
	pyl.figure(1)
	pyl.xlim(-0.2,1.)
	pyl.xlabel('proportion of development')
	if len(Tlist) == 1:
		pyl.ylabel('oxygen (%)')
		pyl.ylim(9,30)
	else:
		pyl.ylabel('temperature (C)')
		pyl.ylim(17,28)	
	pyl.title(sys.argv[1])
	pyl.savefig('curveFitting/curveFitting_'+sys.argv[1]+'_'+spNames+Tnames+'eventTrends.pdf', transparent = True)
	pyl.close()
