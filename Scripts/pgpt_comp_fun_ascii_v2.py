# -*- coding: utf-8 -*-
################### pgpt_comp_kegg.py ############################
# AUTHOR: Sascha Patz at university of Tuebingen, 2022-01-11
# DESCRIPT: Takes blastKOALA, KOfamScan, IMG KEGG, MEGAN KEGG or 
#           PICRUST KEGG files and merge them into one file and returns gene
#           count table reporting PGPTs only.
#CMD: see python3 pgpt_comp_fun_ascii_v2.py -h
##################################################################
vers = "0.1"
# Include standard modules
import argparse, sys, os, io
from ete3 import Tree
from  src.FUNtree import *

# functions

def check_file_exist(file,info):
	if info == "inpfarfact":
		print("Checking file path of factors/traits ...")
	else:
		print("Checking path input data ...")
	if os.path.exists(file):
		if os.path.isfile(file):
			print("... single file found")
			return("F")
		elif os.path.isdir(file):
			print("... path to files found")
			return("P")
	else:
		print("ERROR: File, directory or path does not exist for:"+ file)
		sys.exit()

def read_classif(classif):
	if classif == "PGPT":
		fi = io.open("/nfs/wsi/ab/projects/anupam/server/pb/tools/PGPTpy_v1.4/src/pgpt_class_kegg.txt", "r",encoding='utf8').readlines()
		cl = convert_tsv(fi)
		cl_t = Tree(cl, format = 1)
	return(cl_t)

def read_pfarFactors(pfar):
	pfar_l = io.open(pfar, "r",encoding='utf8').readlines()[1:]
	pf_d = {}
	for line in pfar_l:
		lin_l = line.split(",")
		lin_l_l = lin_l[2].split("#")
		pf_d[lin_l[0]] = str("PGPT;"+"#PGPT;".join(lin_l_l))
	return(pf_d)

def getleafpath(leaf, pfar_d):
	leaf = leaf.split("-")[0]
	return(pfar_d[leaf])

def img2dic(file,form): #form = file/folder
	cnt_d = {}
	tag_d ={}
	fail_l=[]
	if form == "F":
		print("Analyzing IMG KEGG annotations for single file ...")
		fi_l = open(file,"r").readlines()
		if fi_l[0].startswith("gene_oid") or fi_l[0].startswith("Gene ID"):
			fi_l = fi_l[1::]
		fi_s = str(os.path.splitext(os.path.basename(file))[0])
		org_l = [fi_s]
		for line in  fi_l:
			#print(line.split("\t"))
			if "KO:" in line.split("\t")[2]:
				ko = line.split("\t")[2].split(":")[1]
				tag = line.split("\t")[0]
			elif "KO:" in line.split("\t")[6]:
				ko = line.split("\t")[6].split("=")[0].split(":")[1]
				tag = line.split("\t")[1]
			elif "KO:" in line.split("\t")[9]:
				ko = line.split("\t")[9].split(":")[1]
				tag = line.split("\t")[0]
			if ko not in cnt_d:
				cnt_d[ko] = [1]
				tag_d[ko] = [tag]
			elif tag not in tag_d[ko]:
				cnt_d[ko][0] += 1
				tag_d[ko].append(tag)
	else:
		print("Analyzing IMG KEGG annotations for multiple files ...")
		org_l = []
		i_i = 0
		for fi in os.listdir(file):
			fi_l = io.open(os.path.join(file, fi), "r",encoding='utf8').readlines()[1::]
			fi_s = str(os.path.splitext(os.path.basename(fi))[0])
			k_b  = False
			try:
				for line in  fi_l:
					if "KO:" in line.split("\t")[2]:
						ko = line.split("\t")[2].split(":")[1]
						tag = line.split("\t")[0]
					elif "KO:" in line.split("\t")[6]:
						ko = line.split("\t")[6].split("=")[0].split(":")[1]
						tag = line.split("\t")[1]
					elif "KO:" in line.split("\t")[9]:
						ko = line.split("\t")[9].split(":")[1]
						tag = line.split("\t")[0]
					if ko not in cnt_d:
						cnt_d[ko] = create_default_list(i_i,"i")
						cnt_d[ko].append(1)
						tag_d[ko]= create_default_list(i_i,"s")
						tag_d[ko].append([tag])
					elif ko in cnt_d:
						if len(tag_d[ko]) < i_i+1:
							cnt_d[ko].append(1)
							tag_d[ko].append([tag])
						elif tag not in tag_d[ko][i_i]:
							cnt_d[ko][i_i] += 1
							tag_d[ko][i_i].append(tag)
				for k in cnt_d.keys():
					if len(cnt_d[k]) < i_i+1:
						cnt_d[k].append(0)
						tag_d[k].append("na")
				print(" --> File "+fi+" loaded")
				i_i += 1
				org_l.append(fi_s)
			except IndexError:
				fail_l.append(fi)
				print("--> FILE"+fi+" not loaded, as an Error has occurred. Filename written to unprocessed_files.txt")
				pass
	print("... KEGG annotation merged for all files.")
	return(cnt_d,tag_d,org_l,fail_l)

def create_default_list(index,t):
	if index == 0:
		return([])
	elif index > 0:
		if t == "i":
			return([0]*index)
		elif t == "s":
			return(["na"]*index)

def map_annot2class(anc_d,ant_d,class_t,l):
	print("Mapping KEGG annotations to PGPT classification ...")
	for leaf in class_t.traverse():
		#print(leaf.name)
		l_s = leaf.name[-6:] #KEGG id
		if l_s in anc_d: # ID in KEGG annotation count dic
			leaf.add_features(anno_cnt = anc_d[l_s], anno_tag = ant_d[l_s])
		else:
			leaf.add_features(anno_cnt = [0]*l, anno_tag = ["na"]*l)
	for n in class_t.traverse("levelorder"):
		if n.is_leaf() == False:
			anc_l = []
			leaf_l = []
			for leaf in n:
				if leaf.name not in leaf_l:
					leaf_l.append(leaf.name)
					anc_l.append(leaf.anno_cnt)
			n.add_features(anno_cnt = sum_lists(anc_l))
	print("... mapping DONE")
	return(class_t)

def fun_annot(file,classif,inf, form): #form = file/folder
	if inf == "img":
		anc_d, ant_d, org_l, fail_l = img2dic(file, form)
	else:
		print("ERROR: Annotation file format unknonw")
		sys.exit()
	return(map_annot2class(anc_d,ant_d,classif,len(org_l)),org_l,fail_l)

def report_annot(form,classif,orgs,fail,op,pfarfact,query_id):
	fof = io.open(op+"/unprocessed_files.txt", "w",encoding='utf8')
	fof.write("\n".join(fail))
	fof.close()
	pid = []
	for lev in ["Level0","Level1","Level2","Level3","Level4","Level5","Level6"]:
		if form == "tab":
			if lev == "Level6":
				foc = io.open(op+"/level6_kegg2pgpt_cnt.txt", "w",encoding='utf8')
				fot = io.open(op+"/level6_kegg2pgpt_tag.txt", "w",encoding='utf8')
				foc.write("#CLASS\t " + "\t".join(orgs)+"\n")
				fot.write("#CLASS\t " + "\t".join(orgs)+"\n")
				for leaf in classif:
					if leaf.name not in pid:
						foc.write(leaf.name + "\t " + "\t".join([str(x) for x in leaf.anno_cnt])+"\n")
						fot.write(leaf.name + "\t " + "\t".join(str(x) for x in leaf.anno_tag)+"\n")
						pid.append(leaf.name)
				foc.close()
				fot.close()
			else:
				foc = io.open(op+"/"+lev.lower()+"_kegg2pgpt_cnt.txt", "w",encoding='utf8')
				foc.write("#HEAD\t " + "\t".join(orgs)+"\n")
				lev_l = []
				for node in classif.traverse("levelorder"):
					if lev in node.name:
						if node.name not in lev_l:
							foc.write(node.name + "\t " + "\t".join([str(x) for x in node.anno_cnt])+"\n")
							lev_l.append(node.name)
				foc.close()
		elif form == "roary":
			if lev == "Level6":
				foc = io.open(op+"/level6_kegg2pgpt_roary.txt", "w", encoding='utf8')
				fot = io.open(op+"/level6_kegg2pgpt_tag.txt", "w", encoding='utf8')
				foc.write('Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome fragment,Order within fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,' + ','.join(orgs)+'\n')
				fot.write("#CLASS\t " + "\t".join(orgs)+"\n")
				n=1
				for leaf in classif:
					if leaf.name not in pid:
						foc.write(leaf.name + ',,A PGPT gene,100,100,1,1,'+str(n)+',,,,1000,1000,1000,' + ",".join([str(x) for x in leaf.anno_cnt])+"\n")
						fot.write(leaf.name + "\t " + "\t".join(str(x) for x in leaf.anno_tag)+"\n")
						pid.append(leaf.name)
						n+=1
				foc.close()
				fot.close()
			else:
				foc = io.open(op+"/"+lev.lower()+"_kegg2pgpt_roary.txt", "w", encoding='utf8')
				foc.write('Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome fragment,Order within fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,' + ','.join(orgs)+'\n')
				lev_l = []
				n=1
				for node in classif.traverse("levelorder"):
					if lev in node.name:
						if node.name not in lev_l:
							foc.write(node.name + ',,A PGPT gene,100,100,1,1,'+str(n)+',,,,1000,1000,1000,' + ",".join([str(x) for x in node.anno_cnt])+"\n")
							lev_l.append(node.name)
							n+=1
				foc.close()
		elif form =="pfar_kegg":
			query_id = query_id.split("/")[-1].split("_")[0]
			if lev == "Level6":
				print("... converting to pfar_kegg format ...")
				foc = io.open(op+"/"+query_id+"_PGPT-img_kegg.txt", "w",encoding='utf8')
				foc.write("# ["+orgs[0].split("_")[0]+"]\n# This analysis shows the genomic genes that are related to plant growth-promotion (PGPTs)\n") # assuming only one input file for pfar output, currently
				pgpt_d = {}
				pgpt_l = []
				for leaf in classif:
					if leaf.name not in pid:
						#print(leaf.name,leaf.anno_tag)
						pgpt_d[leaf.name] = (leaf.anno_tag,getleafpath(leaf.name, read_pfarFactors(pfarfact))) #dic entry as tuple of 2 lists: 1_GeneTagID, 2_PGPT-paths associated
						#foc.write(leaf.name + "\t " + "\t".join([str(x) for x in leaf.anno_cnt])+"\n")
						#fot.write(leaf.name + "\t " + "\t".join(str(x) for x in leaf.anno_tag)+"\n")
						pid.append(leaf.name)
				for PGPT in sorted(pgpt_d.keys()):
					PGPT_id = PGPT.split("-")[0]
					PGPT_name = PGPT.split("-")[1]
					PGPT_type = pgpt_d[PGPT][1]
					PGPT_hits = pgpt_d[PGPT][0]
					if "na" not in PGPT_hits:
						pgpt_l.append(PGPT_name)
						foc.write("\n   "+"#"*115+"\n   NAME PGPT: "+PGPT_name+"\n   TYPE PGPT: "+PGPT_type+"\n   "+"#"*115+"\n   GENE_FACTOR    GENE_HITS\n   "+"-"*115+"\n")
						foc.write("   "+PGPT_id+" "*4+PGPT_hits[0]+"\n")
						for p in PGPT_hits[1:]:
							 foc.write(" "*18+p+"\n")
				foc.write("\n POSSIBLE PGPTs PRESENT IN THE GENOME:\n =======================================\n")
				print("... summarizing PGPT hits ...")
				for PGPT_name in pgpt_l:
					foc.write(PGPT_name+"\n")
		elif form =="pfar_blast":
			query_id = query_id.split("/")[-1].split("_")[0]
			if lev == "Level6":
				print("... converting to pfar_kegg format ...")
				foc = io.open(op+"/"+query_id+"_PGPT-blhm.txt", "w",encoding='utf8')
				foc.write("# ["+orgs[0].split("_")[0]+"]\n# This analysis shows the genomic genes that are related to plant growth-promotion (PGPTs)\n") # assuming only one input file for pfar output, currently
				pgpt_d = {}
				pgpt_l = []
				for leaf in classif:
					if leaf.name not in pid:
						#print(leaf.name,leaf.anno_tag)
						pgpt_d[leaf.name] = (leaf.anno_tag,getleafpath(leaf.name, read_pfarFactors(pfarfact))) #dic entry as tuple of 2 lists: 1_GeneTagID, 2_PGPT-paths associated
						#foc.write(leaf.name + "\t " + "\t".join([str(x) for x in leaf.anno_cnt])+"\n")
						#fot.write(leaf.name + "\t " + "\t".join(str(x) for x in leaf.anno_tag)+"\n")
						pid.append(leaf.name)
				for PGPT in sorted(pgpt_d.keys()):
					PGPT_id = PGPT.split("-")[0]
					PGPT_name = PGPT.split("-")[1]
					PGPT_type = pgpt_d[PGPT][1]
					PGPT_hits = pgpt_d[PGPT][0]
					if "na" not in PGPT_hits:
						pgpt_l.append(PGPT_name)
						foc.write("\n   "+"#"*115+"\n   NAME PGPT: "+PGPT_name+"\n   TYPE PGPT: "+PGPT_type+"\n   "+"#"*115+"\n   GENE_FACTOR    GENE_HITS\n   "+"-"*115+"\n")
						foc.write("   "+PGPT_id+" "*4+PGPT_hits[0]+"\n")
						for p in PGPT_hits[1:]:
							 foc.write(" "*18+p+"\n")
				foc.write("\n POSSIBLE PGPTs PRESENT IN THE GENOME:\n =======================================\n")
				print("... summarizing PGPT hits ...")
				for PGPT_name in pgpt_l:
					foc.write(PGPT_name+"\n")
def remove_empt(liste):
    while("" in liste):
        liste.remove("")
    while("*" in liste):
        liste.remove("*")
    return(liste)

def sum_lists(data):
	l_l = [sum(j) for j in zip(*data)]
	return(l_l)

if __name__ == "__main__":
	# Initiate the parser
	parser = argparse.ArgumentParser(prog='pgpt_com_fun_ascii-test.py',formatter_class=argparse.MetavarTypeHelpFormatter)

	# Add long and short argument
	parser.add_argument("--input", "-in", type=str, help="input protein/gene kegg file or folder", required=True)
	parser.add_argument("--inpfarfact", "-ipf", type=str, help="PGPTs csv file",default=None)
	parser.add_argument("--classif", "-cl", default = "PGPT", type=str, help="functional classification, def.: PGPT")
	parser.add_argument("--outpath", "-op", default = ".", type=str, help="path to output folder, def.: current directory")
	parser.add_argument("--inform", "-if", type=str, default="img", help="choose input format: [img], def.: img")
	parser.add_argument("--outform", "-of", type=str, default="tab", help="choose output format: [tab/roary/pfar_kegg/pfar_blast], def.: tab")
	parser.add_argument("--version", "-v", help="PGPTpy version", action="store_true")

	# Read arguments from the command line
	args = parser.parse_args()
	# Check for arguments
	if args.version:
		print("Version  "+vers)
	if args.inpfarfact != None:
		form_s = check_file_exist(args.input,"inpfarfact")
	if args.input:
		form_s = check_file_exist(args.input,"input")
	#report_annot(fun_annot(args.input,read_classif(args.classif),args.inform,form_s))
	cla_t = read_classif(args.classif)
	acla_t, o_l, f_l = fun_annot(args.input,cla_t,args.inform,form_s)
	report_annot(args.outform,acla_t, o_l, f_l, args.outpath,args.inpfarfact,args.input)
	#print(acla_t.get_ascii(show_internal=True,attributes=["anno_cnt", "anno_tag"]))
	print("Files written to Output path "+ args.outpath)

