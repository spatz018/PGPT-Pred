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
ap.add_argument("-m", metavar="<mode>", dest="mode",
    help="for img kegg mapper enter IMK while for Blast+Hmmer use BH", required=True)
ap.add_argument("-p", metavar="<FileForPieChart>", dest="FileForPie",
    help="Enter name of file for pie chart", required=True)
args = vars(ap.parse_args())

class fileOutForPlot:
      def __init__(self, df,outFileName,outFileNameforPie):
        self.df = df
        self.outFileName = outFileName
        self.outFileNameforPie = outFileNameforPie

      def pie(self):
        Level2df=self.df[['level2']]
        Level2df['freq'] = Level2df.groupby('level2')['level2'].transform('count')
        Level2df=Level2df.drop_duplicates(['level2'])
        Level2df['freq']=(100. * Level2df['freq'] / Level2df['freq'].sum()).round(0)
        Level2df.to_csv(self.outFileNameforPie,sep="\t",index=False,header=False)

      def Krona(self):
        NameAndType=self.df.groupby(['level1','level2','level3','level4','level5','level6']).size().reset_index(name='freq')
        NameAndType.to_csv(self.outFileName,sep="\t",index=False)

def IMKParser(ReceiveInputFile,ReceiveOutFile,ReceiveFileForPie):
    df = pd.DataFrame()
    fileHandle=open(ReceiveInputFile,"r")
    for line in fileHandle:
        line=line.strip("\n")
       # print(line)
        if ((line.startswith('#')) or ('   =' in line) or ("   --" in line) or (len(line.strip()) == 0 ) or (line.startswith('   #'))):
            continue
        elif 'POSSIBLE PGPTs PRESENT IN THE GENOME:' in line:
            break
        else:
            if line.startswith('   NAME PGPT:'):
                # level2=line.split(": ")[1].split("|")
                level2=line.split(": ")[1]
                # print(level2)
            elif '   TYPE PGPT:' in line:
                listToHold=[]
                level1=line.split(": ")[1].split("#")
                # print(level1)
                # for TypePGPT in range(0,len(level1)):
                    # for NamePGPT in range(0,len(level2)):
                        # listToHold.append(level1[TypePGPT]+";"+level2[NamePGPT])
                for TypePGPT in range(0,len(level1)):
                    listToHold.append(level1[TypePGPT]+";"+level2)
                # print(listToHold)
            elif 'GENE_FACTOR' in line:
                continue
            elif '   PGPT' in line:
                listToHoldUpdate=[]
                level3=re.split(r'\s{2,}', line)
                level3.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        df = df.append({'level1': item[2], 'level2': item[3], 'level3': item[4],'level4': item[5],'level5': item[6],'level6':str(item[7]+'->'+level3[0]),'level7':level3[1]}, ignore_index=True)
                    except:
                        continue
            else:
                level4=re.split(r'\s{2,}', line)
                level4.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        df = df.append({'level1': item[2], 'level2': item[3], 'level3': item[4],'level4': item[5],'level5': item[6],'level6':str(item[7]+'->'+level3[0]),'level7':level4[0]}, ignore_index=True)
                    except:
                        continue
    #df.to_csv(ReceiveOutFile+"main.txt",sep="\t",index=False)
    return(df)


def BHParser(ReceiveInputFile,ReceiveOutFile,ReceiveFileForPie):
    df = pd.DataFrame()
    fileHandle=open(ReceiveInputFile,"r")
    for line in fileHandle:
        line=line.strip("\n")
       # print(line)
        if ((line.startswith('#')) or ('   =' in line) or ("   --" in line) or (len(line.strip()) == 0 ) or (line.startswith('   #'))):
            continue
        elif 'POSSIBLE PGPTs PRESENT IN THE GENOME:' in line:
            break
        else:
            if 'NAME PGPT:' in line:
                # level2=line.split(": ")[1].split("|")
                level2=line.split(": ")[1]
            elif 'TYPE PGPT:' in line:
                listToHold=[]
                level1=line.split(": ")[1].split("#")
                # for TypePGPT in range(0,len(level1)):
                    # for NamePGPT in range(0,len(level2)):
                        # listToHold.append(level1[TypePGPT]+";"+level2[NamePGPT])
                for TypePGPT in range(0,len(level1)):
                    listToHold.append(level1[TypePGPT]+";"+level2)
            elif 'GENE_FACTOR' in line:
                continue
            elif '   PGPT' in line:
                listToHoldUpdate=[]
                level3=re.split(r'\s{2,}', line)
                level3.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        if len(level3)==4:
                            df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level3[1],'level8':level3[2],'level9':level3[3]}, ignore_index=True)
                        elif len(level3)==3:
                            if 'PF' in level3[2]:
                                df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level3[1],'level8':level3[2],'level9':'-'}, ignore_index=True)
                            else:
                                df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level3[1],'level8':'-','level9':level3[2]}, ignore_index=True)
                        else:
                            print("issue")
                    except:
                        continue
            else:
                level4=re.split(r'\s{2,}', line)
                level4.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        
                        if len(level4)==3:
                            df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level4[0],'level8':level4[1],'level9':level4[2]}, ignore_index=True)
                        elif len(level4)==2:
                            if 'PF' in level4[1]:
                                df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level4[0],'level8':level4[1],'level9':'-'}, ignore_index=True)
                            else:
                                df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[6]+'->'+level3[0]),'level7':level4[0],'level8':'-','level9':level4[1]}, ignore_index=True)
                        else:
                            print("issue")
                    except:
                        continue
    df=df.replace('-',np.NaN)
    df=df[df.level7.notnull()]
    return(df)
    
def IMKOldParser(ReceiveInputFile,ReceiveOutFile,ReceiveFileForPie):
    df = pd.DataFrame()
    fileHandle=open(ReceiveInputFile,"r")
    for line in fileHandle:
        line=line.strip("\n")
       # print(line)
        if ((line.startswith('#')) or ('   =' in line) or ("   --" in line) or (len(line.strip()) == 0 ) or (line.startswith('   #'))):
            continue
        elif 'POSSIBLE PGPTs PRESENT IN THE GENOME:' in line:
            break
        else:
            if line.startswith('   NAME PGPT:'):
                # level2=line.split(": ")[1].split("|")
                level2=line.split(": ")[1]
                # print(level2)
            elif '   TYPE PGPT:' in line:
                listToHold=[]
                level1=line.split(": ")[1].split("#")
                # print(level1)
                # for TypePGPT in range(0,len(level1)):
                    # for NamePGPT in range(0,len(level2)):
                        # listToHold.append(level1[TypePGPT]+";"+level2[NamePGPT])
                for TypePGPT in range(0,len(level1)):
                    listToHold.append(level1[TypePGPT]+";"+level2)
                # print(listToHold)
            elif 'GENE_FACTOR' in line:
                continue
            elif '   PGPT' in line:
                listToHoldUpdate=[]
                level3=re.split(r'\s{2,}', line)
                level3.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[7]+'->'+level3[0]),'level7':level3[1]}, ignore_index=True)
                    except:
                        continue
            else:
                level4=re.split(r'\s{2,}', line)
                level4.pop(0)
                for item in listToHold:
                    try:
                        item=item.split(";")
                        df = df.append({'level1': item[1], 'level2': item[2], 'level3': item[3],'level4': item[4],'level5': item[5],'level6':str(item[7]+'->'+level3[0]),'level7':level4[0]}, ignore_index=True)
                    except:
                        continue
    #df.to_csv(ReceiveOutFile+"main.txt",sep="\t",index=False)
    return(df)


def ParsePGPTResult(TakeInputFile,TakeOutFile,TakeMode,TakeFileForPie):
    if TakeMode=='IMK':
        ReturnDataframe=IMKParser(TakeInputFile,TakeOutFile,TakeFileForPie)
        fileOutForPlotHandle=fileOutForPlot(ReturnDataframe,TakeOutFile,TakeFileForPie)
        fileOutForPlotHandle.pie()
        fileOutForPlotHandle.Krona()
        
        
        
        
    elif TakeMode=='BH':
        ReturnDataframe=BHParser(TakeInputFile,TakeOutFile,TakeFileForPie)
        strict = ReturnDataframe[ReturnDataframe.level8.notnull() & ReturnDataframe.level9.notnull()]
        lessStrictBlast = ReturnDataframe[ReturnDataframe.level6.notnull()]
        fileOutForPlotstrict=fileOutForPlot(strict,TakeOutFile+"-blhm.txt",TakeFileForPie+"-blhm.txt")
        fileOutForPlotstrict.pie()
        fileOutForPlotstrict.Krona()
        fileOutForPlotlessStrictBlast =fileOutForPlot(lessStrictBlast,TakeOutFile+"-bl.txt",TakeFileForPie+"-bl.txt")
        fileOutForPlotlessStrictBlast.pie()
        fileOutForPlotlessStrictBlast.Krona()
        
    elif TakeMode=='IMKOld':
        ReturnDataframe=IMKOldParser(TakeInputFile,TakeOutFile,TakeFileForPie)
        fileOutForPlotHandle=fileOutForPlot(ReturnDataframe,TakeOutFile,TakeFileForPie)
        fileOutForPlotHandle.pie()
        fileOutForPlotHandle.Krona()
    



if __name__=="__main__":
    ParsePGPTResult(args['InputFile'],args['OutFile'],args['mode'],args['FileForPie'])
