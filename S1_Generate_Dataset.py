'''
Created on 2017. 11. 7.

@author: mjkwon
'''
import random
import os
import time
from datetime import datetime

#==============================================================================================
# ProcessingTime: Calculate the processing time
#==============================================================================================
def ProcessingTime(startTime):
    endTime=time.time()
    print ('Processing Time:',str(datetime.timedelta(seconds=(endTime-startTime))),'\n')
    startTime=time.time()
    return startTime

#==============================================================================================
# MakeDicHormoneReceptor: Get dictionary whose keys are hormones and values are receptors
#==============================================================================================
def MakeDicHormoneReceptor():
    dicResult = {}
    f = open('./raw_data/Hormone(name)_receptor(symbol).txt')
    for line in f:
        hormone, receptors = line.strip().split('\t')
        dicResult[hormone.lower()] = receptors.split('|')
    return dicResult

#==============================================================================================
# MakeDicDrugTarget: Get dictionary whose keys are drugs and values are target
#=======================================================================================
def MakeDicDrugTarget():
    dicResult = {}
    f = open('./raw_data/Drug(name)_target(symbol).txt')
    for line in f:
        drug, targets = line.strip().split('\t')
        dicResult[drug.lower()] = targets.split('|')
    return dicResult

#==============================================================================================
# MakeDicDiseaseDrug: Get dictionary whose keys are diseases and values are drugs
#==============================================================================================
def MakeDicDiseaseDrug():
    ### Disease (id) : Drugs (name) ###
    dicDssDrg = {}
    f = open('./raw_data/Drug(name)_disease(name).txt')
    for line in f:
        dr, dsss = line.strip().split('\t')
        for dss in dsss.split('|'):
            if dss.lower() not in dicDssDrg.keys(): dicDssDrg[dss.lower()] = [dr.lower()]
            else: dicDssDrg[dss.lower()].append(dr.lower())
    f.close()
    return dicDssDrg

#==============================================================================================
# CallGoldStandardSet: Get dictionary whose keys are diseases and values are hormone-drug pairs
#==============================================================================================
def CallGoldStandardSet(disease):
    dicGoldSet[disease] = set([])
    f = open('./raw_data/Gold standard set.txt')
    for line in f:
        if '#' in line: continue
        dss, hr, rel, dr = line.strip().split('\t')
        if dss.lower() == disease.lower():
            dicGoldSet[disease].add((hr.lower(), dr.lower()))

#==============================================================================================
# GenerateUnlabeledSet: Generate total 25 unlabeled sets for each of 20 target diseases
#==============================================================================================
def GenerateUnlabeledSet(disease):
    # Five sizes of unlabeled samples (i.e. 1,3,5,7 and 10 times of the size of the corresponding gold standard sets)
    for size in [1, 3, 5, 7, 10]:
        # Five unlabeled sets per size
        for i in range(1, 6):

            if os.path.exists("processing") == False: os.mkdir("processing")

            fo = open('./processing/UnlabeledSet/Unlabeled set_%s_%s_%s.txt' % (disease, size, i), 'w+')
            history = {}
            num = size * len(dicGoldSet[disease])
            cnt = 0
            while (True):
                hormone = random.choice(list(dicHormonReceptor.keys()))
                drug = random.choice(list(dicDiseaseDrug[disease.lower()]))
                if drug.lower() not in dicDrugTargets.keys(): continue

                if (hormone.lower(), drug.lower()) in dicGoldSet[disease]: continue
                if (hormone.lower(), drug.lower()) in history.keys(): continue

                fo.write('\t'.join([hormone, 'rel', drug]) + '\n')
                history[(hormone.lower(), drug.lower())] = ''
                cnt += 1
                if cnt == num: break
            fo.close()

def GetListDisease():
    result = []
    f = open('./raw_data/Target diseases(name).txt')
    f.readline()  # Header
    for line in f:
        disease = line.strip().replace('"', '')
        result.append(disease)
    f.close()
    return result

#==============================================================================================
# GenerateDataset: Generate datasets for diseases
#==============================================================================================
def GenerateDataset():
    global dicHormonReceptor, dicDrugTargets, dicDiseaseDrug, dicGoldSet
    dicHormonReceptor = MakeDicHormoneReceptor()
    dicDrugTargets = MakeDicDrugTarget()
    dicDiseaseDrug = MakeDicDiseaseDrug()
    listDisease = GetListDisease()
    dicGoldSet = {}

    for disease in listDisease:

        ### Gold standard set ###
        CallGoldStandardSet(disease)

        ### Unlabeled set ###
        GenerateUnlabeledSet(disease)


    
if __name__ == '__main__':

    startTime = datetime.now()
    print("Start: " + str(startTime))

    GenerateDataset()

    endTime = datetime.now()
    print("End: " + str(endTime))
    print("Process time: " + str(endTime - startTime))