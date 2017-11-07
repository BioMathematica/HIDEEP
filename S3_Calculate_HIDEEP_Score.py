'''
Created on 2017. 11. 7.

@author: mjkwon
'''
from S1_Generate_Dataset import MakeDicHormoneReceptor, GetListDisease
import numpy as np
from sklearn.metrics import roc_auc_score
from datetime import datetime
import os

#==============================================================================================
# GetScore: Get a score for a hormone-drug pair
#==============================================================================================
def GetScore(min_hep, heps):
    starts = set([])
    ends = set([])
    for hep in list(set(heps.split('//'))):
        hep = hep.split('|')
        starts.add(hep[0])
        ends.add(hep[-1])
    score = pow(7, -float(min_hep)) * len(ends) * len(starts)
    return score

#==============================================================================================
# GetAnalysisResult: Call results of network analysis
#==============================================================================================
def GetAnalysisResult(Type,disease, size, i,tg):
    pair, data, target = [], [], []

    if Type == 'Gold standard set':
        f = open('./processing/NetworkAnalysis/Analysis_%s_%s.txt' % (Type,disease))
    else:
        f = open('./processing/NetworkAnalysis/Analysis_%s_%s_%s_%s.txt' % (Type,disease, size, i))

    f.readline()
    for line in f:
        listLine = line.strip().split('\t')

        hormone, drug = listLine[:2]
        # if hormone.lower() not in dicHormonReceptor.keys(): continue

        ## The length of HEP
        if listLine[2] == 'None': lengnHEP = 100.0 # '100.0' represents infinite length
        else: lengnHEP = float(listLine[2]); heps = listLine[3]
        
        ### Score ###
        if lengnHEP == 100.0: score = 0
        else: score = GetScore(lengnHEP, heps)

        pair.append((hormone,drug)); data.append(score), target.append(tg)
    f.close()
    return pair, data, target

#==============================================================================================
# PerformanceEvaluation: Evaluate performance of HIDEEP
#==============================================================================================
def PerformanceEvaluation():
    global dicHormonReceptor
    dicHormonReceptor = MakeDicHormoneReceptor()
    listDisease = GetListDisease()

    if os.path.exists("processing/Scoring") == False: os.mkdir("processing/Scoring")

    fo2 = open('./processing/Scoring/AUC_For_All_Disease_Size.txt', 'w+')
    fo2.write('\t'.join(['Disease', 'Size of unlabeled set', 'iteration', 'AUC']) + '\n')
    for disease in listDisease:
        for size in [1,3,5,7,10]:
            aucScores = []

            for iter in range(1, 6):
                pair_g, data_g, labels_g= GetAnalysisResult('Gold standard set',disease, size,iter,1)
                pair_u, data_u, labels_u= GetAnalysisResult('Unlabeled set',disease, size,iter,0)
    
                auroc_score = roc_auc_score(np.array(labels_g + labels_u), np.array(data_g + data_u))
                aucScores.append(auroc_score)
                fo2.write('\t'.join([disease, str(size), str(iter), str(auroc_score)]) + '\n')
                
                ###  Write scoring result ###
                fo = open(('./processing/Scoring/Scoring_%s_%s_%s.txt') % (disease, str(size), str(iter)), 'w+')
                for i in range(len(pair_g)):
                    h, d = pair_g[i]
                    fo.write('\t'.join([h, d, str(data_g[i]), str(1)]) + '\n')
                for i in range(len(pair_u)):
                    h, d= pair_u[i]
                    fo.write('\t'.join([h, d, str(data_u[i]), str(0)]) + '\n')
                fo.close()

            fo2.write('\t'.join([disease,str(size),"Average",str(float(sum(aucScores) / len(aucScores)))])+'\n')
    fo2.close()

if __name__ == '__main__':
    startTime = datetime.now()
    print ('Start: ' + str(startTime))

    PerformanceEvaluation()

    endTime = datetime.now()
    print ('End: ' + str(endTime))
    print ('Process time: ' + str(endTime - startTime))
