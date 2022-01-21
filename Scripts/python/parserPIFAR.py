import pandas as pd
import numpy as np
import argparse
import re
import os

ap = argparse.ArgumentParser()
ap.add_argument("-i", metavar="<Input File >", dest="InputFile",
    help="Enter the path to Input file ", required=True)
ap.add_argument("-o", metavar="<OutputFileName>", dest="OutFile",
    help="Enter the File Name for Output", required=True)
ap.add_argument("-p", metavar="<FileForPieChart>", dest="FileForPie",
    help="Enter name of file for pie chart", required=True)
args = vars(ap.parse_args())

def ParsePfarResult(ReceiveFile,RecieveOutFilepath,FileForPie):
    df = pd.DataFrame()
    
    fileHandle=open(ReceiveFile,"r")
    for line in fileHandle:
        line=line.strip("\n")
        if (('#' in line) or ('   =' in line) or ("   --" in line)):
            continue
        elif 'POSSIBLE FACTORS PRESENT IN THE GENOME:' in line:
            break
        else:
            if 'NAME INTERACTION AND VIRULENCE FACTOR' in line:
                level2=line.split(": ")[1]
                # print(level2)
            elif 'TYPE INTERACTION AND VIRULENCE FACTOR' in line:
                level1=line.split(": ")[1]
                # print(level1)
            elif 'GENE_FACTOR' in line:
                continue
            elif 'PVF' in line:
                level3=re.split(r'\s{2,}', line)
                # print(level1,level2,level3[1],level3[2],level3[3],level3[4])
                df = df.append({'level1': level1, 'level2': level2, 'level3': level3[1],'level4': level3[2],'level5': level3[3],'level6':level3[4]}, ignore_index=True)
            else:
                level4=re.split(r'\s{2,}', line)
                # print(level4)
                if len(level4)>1:
                    df = df.append({'level1': level1, 'level2': level2, 'level3': level3[1],'level4': level4[1],'level5': level4[2],'level6':level4[3]}, ignore_index=True)
    df=df.replace('-',np.NaN)
    df=df[df.level4.notnull()]
    VirulenceNameAndType=df.groupby(['level1','level2','level3']).size().reset_index(name='freq')
    VirulenceNameDataFrame=df[['level1']]
    VirulenceNameDataFrame['freq'] = VirulenceNameDataFrame.groupby('level1')['level1'].transform('count')
    VirulenceNameDataFrame=VirulenceNameDataFrame.drop_duplicates(['level1'])
    VirulenceNameDataFrame['freq']=(100. * VirulenceNameDataFrame['freq'] / VirulenceNameDataFrame['freq'].sum()).round(0)
    VirulenceNameDataFrame.to_csv(FileForPie,sep="\t",index=False,header=False)
    VirulenceNameAndType.to_csv(RecieveOutFilepath,sep="\t",index=False)
    strict = df[df.level5.notnull() & df.level6.notnull()]
    lessStrictPF = df[df.level5.notnull() & df.level6.isnull()]
    # lessStrictBlast = df[df.level5.isnull() & df.level6.notnull()]
    lessStrictBlast = df[df.level4.notnull() & df.level6.notnull()] #####Update 11/oct/2021
    forBlastkrona=lessStrictBlast.groupby(['level1','level2','level3']).size().reset_index(name='freq')
    forBlastkrona.to_csv(RecieveOutFilepath.replace("FileForKrona","FileForKronaBlast"),sep="\t",index=False)
    if lessStrictBlast.empty:
        #print('DataFrame is empty!')
        f = open(FileForPie.replace("blhm","bl"), 'w')
        f.write("null"+"\t"+"100"+"\n")
        f.close() 
    else:
        lessStrictBlastForPie=lessStrictBlast[['level1']]
        lessStrictBlastForPie
        lessStrictBlastForPie['freq'] = lessStrictBlastForPie.groupby('level1')['level1'].transform('count')
        lessStrictBlastForPie=lessStrictBlastForPie.drop_duplicates(['level1'])
        lessStrictBlastForPie['freq']=(100. * lessStrictBlastForPie['freq'] / lessStrictBlastForPie['freq'].sum()).round(0)
        lessStrictBlastForPie.to_csv(FileForPie.replace("blhm","bl"),sep="\t",index=False,header=False)
 
    strict.to_csv(os.path.join(os.path.dirname(RecieveOutFilepath),"strict.txt"),sep="\t",index=False)
    lessStrictPF.to_csv(os.path.join(os.path.dirname(RecieveOutFilepath),"lessStrictPF.txt"),sep="\t",index=False)
    lessStrictBlast.to_csv(os.path.join(os.path.dirname(RecieveOutFilepath),"lessStrictBlast.txt"),sep="\t",index=False)


       
if __name__=="__main__":
    ParsePfarResult(args['InputFile'],args['OutFile'],args['FileForPie'])
