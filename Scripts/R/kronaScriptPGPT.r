suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(psadd))
library(argparse)
library(stringr)
source("phyloseq-extended/R/load-extra-functions.R")

parser <- ArgumentParser()
parser$add_argument("-i", "--InputFile" , default=TRUE,
    help="Path to input file in krona require format")
parser$add_argument("-o", "--OutputHTMLFile" , default=TRUE,
    help="Path to out file")
parser$add_argument("-m", "--mode" , default=TRUE,
    help="for img kegg mapper enter IMK while for Blast+Hmmer use BH")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

options(browser = "false")


parsing_to_require_formatIMG <- function(RequireInputFile) {
	InputFileForPlot<-read.table(RequireInputFile,sep="\t",header=TRUE)
	
	rownames(InputFileForPlot) <- paste0("OTU", 1:nrow(InputFileForPlot))
	level_name<-cbind(InputFileForPlot[,1:6])
	InputFileForPlot[,1:6]<-NULL
	LevelAbundace = otu_table(InputFileForPlot, taxa_are_rows = TRUE)
	LevelName=tax_table(as.matrix(level_name))
	Mapping <- data.frame (Name  = c("freq"))
	rownames(Mapping)<-Mapping[,1]
	metadata<-sample_data(Mapping)
	print(sample_data(metadata))
	ObjectPhyseq = merge_phyloseq(LevelAbundace, LevelName,metadata)
	count_to_prop <- function(x) {return( x / sum(x) )}
	physeqTrans <- transform_sample_counts(ObjectPhyseq, count_to_prop)
	return(physeqTrans)

}

parsing_to_require_format <- function(RequireInputFile) {
	InputFileForPlot<-read.table(RequireInputFile,sep="\t",header=TRUE)
	
	rownames(InputFileForPlot) <- paste0("OTU", 1:nrow(InputFileForPlot))
	level_name<-cbind(InputFileForPlot[,1:3])
	InputFileForPlot[,1:3]<-NULL
	LevelAbundace = otu_table(InputFileForPlot, taxa_are_rows = TRUE)
	LevelName=tax_table(as.matrix(level_name))
	Mapping <- data.frame (Name  = c("freq"))
	rownames(Mapping)<-Mapping[,1]
	metadata<-sample_data(Mapping)
	print(sample_data(metadata))
	ObjectPhyseq = merge_phyloseq(LevelAbundace, LevelName,metadata)
	count_to_prop <- function(x) {return( x / sum(x) )}
	physeqTrans <- transform_sample_counts(ObjectPhyseq, count_to_prop)
	return(physeqTrans)

}


ForKronaPlot <- function(TakeInputFile,TakeOutFileName,TakeMode) {
	if(TakeMode=='IMK'){
		plot_krona(parsing_to_require_formatIMG(TakeInputFile),TakeOutFileName,"Name",trim=T)
	}
	else {
		plot_krona(parsing_to_require_formatIMG(TakeInputFile),TakeOutFileName,"Name",trim=T)
	}
	
	# plot_krona(parsing_to_require_format(str_replace(TakeInputFile, "FileForKrona","FileForKronaBlast")),str_replace(TakeOutFileName,"_FileForKrona","_FileForKronaBlast"),"Name",trim=T)

}

ForKronaPlot(args$InputFile,args$OutputHTMLFile,args$mode)