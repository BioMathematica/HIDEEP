[Synopsis]
HIDEEP: a systems approach to predict hormone impacts on drug efficacy based on effect paths

[Environment]
 python-3.6.1

[How to use]
 Input Files
 - Get molecular interactions (e.g. "./raw_data/Molecular_network.txt")
  "Left type [delimiter:tab] Left [delimiter:tab] Relation [delimiter:tab] Right type [delimiter:tab] Right [delimiter:tab] Source DB"
  (e.g. "G [delimiter:tab] FLNC [delimiter:tab] interact [delimiter:tab] G [delimiter:tab] MAP2K4 [delimiter:tab] BioGRID")
 - Get hormone-receptor relations (e.g. "./raw_data/Hormone(name)_receptor(symbol).txt")
  "Hormone Name [delimiter:tab] Gene Symbol"
  (e.g. "Corticosterone NR3C1|NR3C2")
 - Get drug-target relations (e.g. "./raw_data/Drug(name)_target(symbol).txt")
  "Drug Name [delimiter:tab] Gene Symbol"
  (e.g. "Lapatinib [delimiter:tab] ERBB2|EGFR")
 - Get drug-disease relations (e.g. "./raw_data/Drug(name)_disease(name).txt")
  "Drug Name [delimiter:tab] Disease Name"
  (e.g. "Desflurane Reperfusion [delimiter:tab] Injury|Hypertension|Myocardial Ischemia")
 - Get disease-gene relations (e.g. "./raw_data/Disease(name)_gene(symbol).txt")
  "Disease Name [delimiter:tab] Gene Symbol"
  (e.g. "Postoperative Complications [delimiter:tab] NPPB|BCHE|POMC")
 - Get gold standard samples (e.g. "./raw_data/Gold standard set.txt")
  "Disease Name [delimiter:tab] Hormone Name [delimiter:tab] Interaction [delimiter:tab] Drug Name"
  (e.g. "Hypertension [delimiter:tab] Melatonin [delimiter:tab] inhibits [delimiter:tab] Amlodipine")
 - Get target disease list (e.g. "./raw_data/Target diseases(name).txt")
  "Disease Name"
  (e.g. "Hypertension")

 Run Code
 - "python HIDEEP.py"

 Output Files
 - Generate five unlabeled datasets per size (e.g. "./processing/UnlabeledSet/Unlabeled set_Hypertension_1_1.txt")
  "Hormone [delimiter:tab] relation [delimiter:tab] Drug"
  (e.g. "Serotonin [delimiter:tab] rel [delimiter:tab] Carteolol")
 - Write analysis results of hormone-drug pairs for each disease (e.g. "./processing/NetworkAnalysis/Analysis_Unlabeled set_Hypertension_1_1.txt")
  "Hormone [delimiter:tab] Drug [delimiter:tab] length of HEP [delimiter:tab] HEPs"
  (e.g. "Serotonin [delimiter:tab] Carteolol [delimiter:tab] 1 [delimiter:tab] HTR5A|GNAS//HTR6|GNAS//HTR7|GNAS//HTR4|GNAS")
 - Write scores of hormone-drug pairs for each disease (e.g. "./processing/Scoring/Scoring_Hypertension_1_1.txt")
  "Hormone [delimiter:tab] Drug [delimiter:tab] Score [delimiter:tab] Label"
  (e.g. "Serotonin [delimiter:tab] Carteolol [delimiter:tab] 0.57 [delimiter:tab] 0")
 - Write AUROCs for five datasets per size (i.e. "./processing/Scoring/AUC_For_All_Disease_Size.txt")
  "Disease [delimiter:tab] Size of unlabeled sets [delimiter:tab] Iteration [delimiter:tab] AUC"
  (e.g. "Hypertension [delimiter:tab] 1 [delimiter:tab] 1 [delimiter:tab] 0.89")