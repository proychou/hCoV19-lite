# Generate consensus sequence from bam file
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(Biostrings);
library(RCurl);

#Get latest stable version of wgs_functions.R from github
# source('./wgs_functions.R'); #or locally
script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R',
							 ssl.verifypeer=FALSE)
eval(parse(text=script));

#Get args from command line: needs bamfname and sampname
args<-(commandArgs(TRUE));
if(length(args)==0){
	print("No arguments supplied.")
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
		print(args[[i]])
	}
}


#Generate consensus sequence -- this function is within wgs_functions.R (common to several viral WGS pipelines)
conseq<-generate_consensus(bamfname)

if(!is.na(conseq)){
	if(!dir.exists('./consensus_seqs')) dir.create('./consensus_seqs');
	writeXStringSet(conseq,file=paste('./consensus_seqs/',sampname,'.fasta',sep=''),format='fasta');
	
}else{
	print('Failed to generate consensus sequences.')
}