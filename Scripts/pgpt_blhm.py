# -*- coding: utf-8 -*-
################### pgpt_blhm.py ############################
# AUTHOR: Sascha Patz at university of Tuebingen, 2022-01-11
# DESCRIPT: Takes protein fasta input file and runs blastp agains pgpt ontology 
#           ane executes hmmsearch. Generate a pfar-like output format with all
#           hits
# CMD: see python3 pgpt_blhm.py -h
##################################################################
vers = "0.1"
# Include standard modules
import argparse,sys,os,io
from pyfasta import Fasta
from concurrent.futures.thread import ThreadPoolExecutor


def readfactors(factor):
	fi = io.open(factor,"r",encoding='utf8').readlines()
	fac_d = {}
	for line in fi:
		lin = line.rstrip().split(",")
		fid = lin[0] # factor id
		#print(fid)
		fna = lin[1] # factor name
		fty = lin[2] # factor type (classification path, paths sepatrated by \t, single path classes separated by ;)
		fpf = lin[3] # factors pfam domains (not all factors associated to pfams)
		fac_d[fid]=[fna,fty,fpf]
	#print(d)
	return(fac_d)

def hmmscan(seq_file,op,query_id,pf_hmm,eval,threads):
	os.system("hmmsearch --tblout "+str(op+"/").replace("//","/")+query_id+".hmo -E "+eval+" --cpu "+threads+" "+pf_hmm+" "+seq_file)
	fh_l = io.open(str(op+"/").replace("//","/")+query_id+".hmo","r",encoding='utf8').readlines()[3:]
	return(gethmmhits(fh_l))

def blastseq(seq_file,op,fac_db,max_tseq,al_n,cov,eval,threads):
	os.system("blastp -query "+seq_file+" -db "+fac_db+" -out "+seq_file+".blast7 -num_alignments "+str(al_n)+"-qcov_hsp_perc "+str(cov)+" -evalue "+str(eval)+" -outfmt 7 -num_threads "+str(threads))

def reafblastfile(query_id,op):
	fb_l = io.open(str(op+"/").replace("//","/")+query_id+".blast7","r",encoding='utf8').readlines()
	return(getblasthits(fb_l)) # returns dictionary of PGPT ids containing each a list of two lists, one for protein hits and another for e-Value

def getblasthits(bl_l): #takes the blast output file and generates a dictionary of PGPT ids containing each a list of two lists, one for protein hits and another for e-Value
	bl_d = {}
	for line in bl_l:
		if "PGPT0" in line:
			lin_l = line.rstrip().split("\t")
			pgpt = lin_l[1].split("-")[0]
			if pgpt not in bl_d.keys():
				bl_d[pgpt] = [[lin_l[0]],[lin_l[10]]]
			else:
				bl_d[pgpt][0].append(lin_l[0])
				bl_d[pgpt][1].append(lin_l[10])
	return(bl_d)

def gethmmhits(hmm_l):
	hmm_d = {}
	for line in hmm_l:
		if line[0] != "#":
			lin_l = line.rstrip().split("  PF")
			#print(lin_l)
			prot = lin_l[0].split(" ")[0]
			pfam = "PF"+lin_l[1].split(".")[0]
			if prot not in hmm_d.keys():
				hmm_d[prot] = [pfam]
			else:
				hmm_d[prot].append(pfam)
	#print(hmm_d)
	return(hmm_d)

def checkPFAMS(pgpt,hmm_d,acc,factors_db):
	try:
		pfam_l = hmm_d[acc]
		rp_l = []
		for pfam in pfam_l:
			if pfam in factors_db[pgpt][2]:
				 rp_l.append(pfam)
		return(rp_l)
	except KeyError:
		return([])

def minLengpos(el,wid):
	n = wid-len(el)
	if n <= 0:
		return(2)
	else:
		return(n)

def result(blast_db,hmm_db,factors_db,query_id,op):
	foc = io.open(str(op+"/").replace("//","/")+query_id+"_PGPT-blhm.txt", "w",encoding='utf8')
	foc.write("# ["+query_id+"]\n# This analysis shows the genomic genes that are related to plant growth-promotion (PGPTs)\n") # assuming only one input file for pfar output, currently
	pgpt_l = []
	a=0
	for pgpt in sorted(blast_db.keys()):
		PGPT_name = factors_db[pgpt][0]
		PGPT_type = factors_db[pgpt][1]
		foc.write("\n   "+"#"*115+"\n   NAME PGPT: "+PGPT_name+"\n   TYPE PGPT: "+PGPT_type+"\n   "+"#"*115+"\n   GENE_FACTOR    GENE_HITS"+" "*16+"PFAMS"+" "*65+"BLAST\n   "+"-"*115+"\n")
		foc.write("   "+pgpt+" "*4+blast_db[pgpt][0][0]+" "*minLengpos(blast_db[pgpt][0][0],25)+" ".join(checkPFAMS(pgpt,hmm_db,blast_db[pgpt][0][0],factors_db))+" "*minLengpos(" ".join(checkPFAMS(pgpt,hmm_db,blast_db[pgpt][0][0],factors_db)),70)+blast_db[pgpt][1][0]+"\n")
		c = 1
		for acc in blast_db[pgpt][0][1:]:
			foc.write(" "*18+acc+" "*minLengpos(acc,25)+" ".join(checkPFAMS(pgpt,hmm_db,acc,factors_db))+" "*minLengpos(" ".join(checkPFAMS(pgpt,hmm_db,acc,factors_db)),70)+blast_db[pgpt][1][blast_db[pgpt][0].index(acc)]+"\n")
			c+= 1
		pgpt_l.append(PGPT_name+"\t"+str(c))
		a += c
	foc.write("\n POSSIBLE "+str(a)+"GENES ASSOCIATED TO  PGPTs (BY ALL BLASTp HITS) PRESENT IN THE GENOME:\n ===================================================================================\n")
	print("... summarizing PGPT hits ...")
	for PGPT_name in pgpt_l:
		foc.write(PGPT_name+"\n")
	foc.close()

if __name__ == "__main__":
	# Initiate the parser
	parser = argparse.ArgumentParser(prog='pgpt_blhm.py',formatter_class=argparse.MetavarTypeHelpFormatter)

	# Add long and short argument
	parser.add_argument("--input", "-in", type=str, help="input protein fasta file", required=True)
	parser.add_argument("--outpath", "-op", default = ".", type=str, help="path to output folder, def.: current directory")
	parser.add_argument("--version", "-v", help="PGPTpy version", action="store_true")

	# Read arguments from the command line
	args = parser.parse_args()
	# Check for arguments
	if args.version:
		print("Version  "+vers)
	seq_file = args.input
	op = args.outpath
	query_id = seq_file.split("/")[-1].split("_")[0]
	fac_db = "./factors/PGPT-blastpdb/nitrogenasePGPT"
	fac_inf = "./factors/PlantGrowthPromotingTraits.csv"
	pfam_hmm = "./factors/Pfam-A.hmm"
	o_hfname = seq_file+".hmo"
	threads = str(8)
	pthr = str(2)
	eval = str(1e-5)
	cov = str(80)
	max_tseq = str(1)
	al_n =str(1)
	fac_dic = readfactors(fac_inf)
	os.system("faSplit sequence "+seq_file+" 12 "+str(op+"/").replace("//","/")+query_id+"_")
	fi_l = []
	for file in os.listdir(str(op+"/").replace("//","/")):
		if file.endswith(".fa"):
			fi_l.append(os.path.join(str(op+"/").replace("//","/"), file))
	hmm_dic = hmmscan(seq_file,op,query_id,pfam_hmm,cov,eval,threads)
	with ThreadPoolExecutor(max_workers=20) as executor:
		#hmm_dic = hmmscan(seq_file,op,query_id,pfam_hmm,eval,threads)
		for s_file in fi_l:
			executor.submit(blastseq,s_file,op,fac_db,max_tseq,al_n,eval,pthr)
	os.system("find "+str(op+"/").replace("//","/")+" -type f -name '*.fa.blast7' -exec cat {} + >> "+str(op+"/").replace("//","/")+query_id+".blast7")
	os.system("rm "+str(op+"/").replace("//","/")+"*.fa.blast7")
	bla_dic = reafblastfile(query_id,op)
	result(bla_dic,hmm_dic,fac_dic,query_id,op)

