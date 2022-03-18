#!/usr/bin/env python
import argparse,time
from Bio import SeqIO
from deepRloopPre._version import __version__
def GetParser():
	RequiredArgs = GetRequiredArgs()
	OptionalArgs = GetOptionalArgs()
	parser = argparse.ArgumentParser(description="A tool for predicting R-loop profiles and R-loop locations.",formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='$ deepRloopPredict.py --h5 model.hdf5 --fasta genome.fasta --prefix test',parents=[RequiredArgs,OptionalArgs],add_help=False)
	return parser
def GetRequiredArgs():
	parser = argparse.ArgumentParser(add_help=False)
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--h5',dest="H5",required=True,help="Saved model")
	required.add_argument('--fasta',dest="Fasta",required=True,help='FASTA format file, usually a genome FASTA file.')
	required.add_argument("--prefix",dest="Prefix",required=True)
	return parser
def GetOptionalArgs():
	parser = argparse.ArgumentParser(add_help=False)
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument("--help", "-h", action="help",help="show this help message and exit")
	optional.add_argument('--chromsize',dest="ChromeSize",default=None,help='If this file is not provided, all the sequences in FASTA file will be predicted. If this file is provided, only the sequences contained in this file will be predicted. chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.')
	optional.add_argument("--strand",dest="Strand",default="fwd",choices=["fwd","rev"],help="fwd:Predict current sequence. rev:Predict the reverse complement sequence.")
	optional.add_argument("--threshold",dest="Threshold",default=0.5,help="The threshold score is used to judge whether an region is an R-loop region, the value range is from 0 to 1.",type=float)
	optional.add_argument("--mem",dest="Memery",default="high",choices=["high","low"],help="low option can reduce memory, but the program will run slower.")
	optional.add_argument('--version', action='version',version='deepRloopPre {}'.format(__version__))
	return parser
if __name__ == '__main__':
	args=GetParser().parse_args()
	Genome=args.Fasta
	Prefix=args.Prefix
	if not args.ChromeSize:
		temp=str(time.time()).split(".")[0]
		ChromeSize="%s_chrom_size_%s.txt"%(Prefix,temp)
		Fr=open(ChromeSize,"w")
		for record in SeqIO.parse(Genome,"fasta"):
			Fr.write(record.id+"\t"+str(len(record.seq))+"\n")
		Fr.close()
	else:
		ChromeSize=args.ChromeSize
	Strand=args.Strand
	Threshold=float(args.Threshold)
	Mem=args.Memery
	Win=128000
	Scale=128
	from deepRloopPre import RloopModel
	from tensorflow.keras.models import load_model
	Model=load_model(args.H5,{'Rsquare':RloopModel.Rsquare})
	RloopModel.Predict(Model,Win,Genome,ChromeSize,Strand,Prefix,Threshold,Scale,Mem)