'''
Created on 2017. 11. 7.

@author: mjkwon
'''
from S1_Generate_Dataset import MakeDicHormoneReceptor, MakeDicDrugTarget, MakeDicDiseaseDrug, GetListDisease
from datetime import datetime
import os

#==============================================================================================
# MakeDicNetworkGraph: Get dictionary for the molecular network
#==============================================================================================
def MakeDicNetworkGraph():
    graph = {}
    f = open('./raw_data/Molecular_network.txt')
    for line in f:
        leftType, left, rel, rightType, right, source = line.strip().split('\t')

        # Uni-directional interaction
        if left in graph.keys():
            graph[left].add(right)
        else:
            graph[left] = set([right])

        # Bi-directional interaction
        if rel == 'interact':
            if right in graph.keys():
                graph[left].add(left)
            else:
                graph[right] = set([left])
    f.close()
    return graph

#==============================================================================================
# MakeDicDiseaseGene: Get dictionary whose keys are diseases and values are hormone-drug pairs
#==============================================================================================
def MakeDicDiseaseGene():
    dicResult = {}
    f = open('./raw_data/Disease(name)_gene(symbol).txt')
    for line in f:
        disease, genes = line.strip().split('\t')
        dicResult[disease.lower()] = genes.split('|')
    f.close()
    return dicResult

#==============================================================================================
# GetTargetAndReceptor: Get list for hormone's receptors or drug's targets
#==============================================================================================
def GetTargetAndReceptor(drug_targets,hormone_receptors,entity):
    result = []
    entity = entity.strip().lower().replace('"','')
    if entity in hormone_receptors.keys(): result += hormone_receptors[entity]
    if entity in drug_targets.keys(): result += drug_targets[entity]
    result = list(set(result))
    return result

#==============================================================================================
# bfs & backtrace: Get list for shortest paths from a start node to an end node
#==============================================================================================
def bfs(graph, start, end):
    visited = set([start])
    parent = {}
    parent[start] = set()
    queue = set([start])
    while queue:
        queue2 = set()
        for node in queue:
            if node == end:
                return backtrace(parent, start, end)
            for adjacent in graph.get(node, []):
                if (adjacent not in visited) or (adjacent in visited and parent[adjacent].isdisjoint(queue) == False):
                    if adjacent in parent.keys():
                        parent[adjacent].add(node)
                    else:
                        parent[adjacent] = set([node])
                    queue2.add(adjacent)
                    visited.add(adjacent)
        queue = queue2

def backtrace(parent, start, end):
    paths = [[end]]
    while paths[0][-1] != start:
        paths2 = []
        for each_path in paths:
            for each_parent in parent[each_path[-1]]:
                new_path = each_path + [each_parent]
                paths2.append(new_path)
        paths = paths2
    for each_path in paths:
        each_path.reverse()
    return paths

#==============================================================================================
# AnalyzeDEP: Get list for shortest paths from a start node to an end node
#==============================================================================================
def AnalyzeDEP(library, graph, drugTargets, diseaseGenes):
    drug_effective_paths = []
    path_lengths = []
    for target in drugTargets:
        closest_distance = 100
        paths = []
        for diseaseGene in diseaseGenes:
            # if target == diseaseGene: paths = [[target]]; break
            if (target, diseaseGene) in library.keys():
                tg_shortest_paths = library[(target, diseaseGene)]
            else:
                tg_shortest_paths = bfs(graph, target, diseaseGene)
                library[(target, diseaseGene)] = tg_shortest_paths

            if tg_shortest_paths == None: continue
            else: spl = len(tg_shortest_paths[0]) - 1

            if closest_distance == spl:
                paths += tg_shortest_paths

            elif closest_distance > spl:
                closest_distance = spl
                paths = tg_shortest_paths

        if paths == []: continue
        drug_effective_paths += paths
        path_lengths.append(closest_distance)

    if drug_effective_paths == []:
        return [None, drug_effective_paths, str(100)]
    else:
        return [len(drug_effective_paths), drug_effective_paths, str(float(sum(path_lengths)) / float(len(path_lengths)))]

#==============================================================================================
# AnalyzeHEP: Get list for the shortest paths from receptors to the nearest molecule of the DEPs
#==============================================================================================
def AnalyzeHEP(library, graph, deps, receptors):
    final_SPL, numOfHEP, hep, dep = 100, 0, [], []

    for one_DEP in deps:
        for end in one_DEP:
            if type(end) == list: continue
            for receptor in receptors:
                if (receptor, end) in library.keys():
                    all_paths = library[(receptor, end)]
                else:
                    all_paths = bfs(graph, receptor, end)
                    library[(receptor, end)] = all_paths

                if all_paths == None:
                    continue
                else:
                    spl = len(all_paths[0]) - 1

                if spl < final_SPL:
                    final_SPL = spl
                    hep = ['|'.join(path) for path in all_paths]
                    dep = ['|'.join(one_DEP)]
                    numOfHEP = len(all_paths)
                elif spl == final_SPL:
                    hep += ['|'.join(path) for path in all_paths]
                    dep += ['|'.join(one_DEP)]
                    numOfHEP += len(all_paths)

    if final_SPL == 100:
        return [None, None, None, None]
    else:
        lenOfHEP = final_SPL
        return [numOfHEP, lenOfHEP, hep, dep]

#==============================================================================================
# CallGoldStandardSamples: Get dictionary whose keys are gold standard hormone-drug pairs
#==============================================================================================
def CallGoldStandardSamples(disease):
    dicResult = {}
    f = open('./raw_data/Gold standard set.txt')
    for line in f:
        if '#' in line: continue
        dss, hr, rel, dr = line.strip().split('\t')
        if dss.lower() == disease.lower():
            dicResult[(hr.lower(), dr.lower())] = ''
    f.close()
    return dicResult

#==============================================================================================
# CallUnlabeledSamples: Get dictionary whose keys are unlabeled hormone-drug pairs
#==============================================================================================
def CallUnlabeledSamples(disease,fileName):
    dicResult = {}
    f = open(fileName)
    for line in f:
        hr, rel, dr = line.strip().split('\t')
        dicResult[(hr.lower(),dr.lower())] = ''
    f.close()
    return dicResult

#==============================================================================================
# NetworkAnalysis: Analyze DEPs and HEPs for each of gold/unlabeled hormone-drug pairs
#==============================================================================================
def NetworkAnalysis():
    global dicHormonReceptor, dicDrugTargets, dicDiseaseDrug, dicGoldSet, dicDiseaseGene, graph, library
    dicHormonReceptor = MakeDicHormoneReceptor()
    dicDrugTargets = MakeDicDrugTarget()
    dicDiseaseDrug = MakeDicDiseaseDrug()
    dicDiseaseGene = MakeDicDiseaseGene()
    listDisease = GetListDisease()
    graph = MakeDicNetworkGraph()
    library = {}

    for disease in listDisease:
        for Type in ['Gold standard set','Unlabeled set']:
            listDssGenes = dicDiseaseGene[disease.lower()]

            if Type == 'Gold standard set':
                sizeRange = range(1,2)
                iterRange = range(1,2)
            else:
                sizeRange = [1,3,5,7,10]
                iterRange = range(1,6)
                fName = './processing/UnlabeledSet/%s_%s_%s_%s.txt'

            if os.path.exists('./processing/NetworkAnalysis/')==False: os.mkdir('./processing/NetworkAnalysis/')

            for size in sizeRange:
                for i in iterRange:
                    if Type == 'Gold standard set':
                        fo = open('./processing/NetworkAnalysis/Analysis_%s_%s.txt'%(Type,disease), 'w+')
                        dicSamples = CallGoldStandardSamples(disease)
                    else:
                        fo = open('./processing/NetworkAnalysis/Analysis_%s_%s_%s_%s.txt' % (Type,disease,str(size),str(i)), 'w+')
                        dicSamples = CallUnlabeledSamples(disease,fName % (Type,disease,size,i))

                    fo.write('hormone\tdrug\tlength of hep\n')
                    for hormone,drug in dicSamples.keys():

                        receptors = GetTargetAndReceptor(dicDrugTargets,dicHormonReceptor,hormone)
                        targets = GetTargetAndReceptor(dicDrugTargets,dicHormonReceptor,drug)

                        ### DEP ###
                        numOfDEP, deps, mean_pathLength = AnalyzeDEP(library, graph, targets, listDssGenes)

                        if deps == [] or deps == None:
                            fo.write('\t'.join([hormone,drug,'None', 'None'])+'\n')
                            continue
                        if receptors == []: fo.write('\t'.join([hormone, drug, 'None', 'None']) + '\n'); continue

                        ### HEP ###
                        numOfHEP, lenOfHEP, hep, dep = AnalyzeHEP(library,graph,deps,receptors)

                        if [numOfHEP, lenOfHEP, hep, dep] == [None,None,None,None]:
                            fo.write('\t'.join([hormone,drug,'None', 'None'])+'\n')
                            continue

                        ### Writing###
                        fo.write('\t'.join([hormone,drug,str(lenOfHEP),'//'.join(hep)])+'\n')
                    fo.close()

if __name__ == '__main__':
    startTime = datetime.now()
    print ('Start: '+str(startTime))

    NetworkAnalysis()

    endTime = datetime.now()
    print ('End: '+str(endTime))
    print ('Process time: '+str(endTime-startTime))
    