'''
Created on 2017. 11. 7.

@author: mjkwon
'''
import os

if __name__ == '__main__':

    print ("Start of S1_Generate_Dataset.py")
    os.system("python S1_Generate_Dataset.py")
    print ("End of S1_Generate_Dataset.py")

    print ("Start of S2_Analyze_DEP_And_HEP.py")
    os.system("python S2_Analyze_DEP_And_HEP.py")
    print("End of S2_Analyze_DEP_And_HEP.py")

    print ("Start of S3_Calculate_HIDEEP_Score.py")
    os.system("python S3_Calculate_HIDEEP_Score.py")
    print("End of S3_Calculate_HIDEEP_Score.py")