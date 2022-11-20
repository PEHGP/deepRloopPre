#!/usr/bin/env python
import argparse,time
import numpy as np
from deepRloopPre._version import __version__
from deepRloopPre import RloopDeal
from Bio import SeqIO
def GetParser():
	RequiredArgs = GetRequiredArgs()
	OptionalArgs = GetOptionalArgs()
	parser = argparse.ArgumentParser(description="format the data for training.",formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='$ deepRloopData.py --fasta genome.fasta --drip test_drip.bdg --qpois test_qpois.bdg --prefix test',parents=[RequiredArgs,OptionalArgs],add_help=False)
	return parser
def GetRequiredArgs():
	parser = argparse.ArgumentParser(add_help=False)
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--fasta',dest="Fasta",required=True,help='FASTA format file, usually a genome FASTA file.')
	required.add_argument('--drip',dest="Drip",required=True,help='R-loop profiles in sorted bdg format obtained from ssDRIP-seq.')
	required.add_argument('--qpois',dest="Qpois",required=True,help='-np.log10(qvalue) in sorted bdg format obtained from macs2.')
	required.add_argument("--prefix",dest="Prefix",required=True)
	return parser
def GetOptionalArgs():
	parser = argparse.ArgumentParser(add_help=False)
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument("--help", "-h", action="help",help="show this help message and exit")
	optional.add_argument('--chromsize',dest="ChromeSize",default=None,help='If this file is not provided, all the sequences in FASTA file will be formated. If this file is provided, only the sequences contained in this file will be formated. chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.')
	optional.add_argument("--strand",dest="Strand",default="fwd",choices=["fwd","rev"],help="fwd:format current sequence. rev:format the reverse complement sequence.")
	optional.add_argument("--logq",dest="Qscore",default=2,help="The threshold score (-np.log10(qvalue)) is used to judge whether an region is an R-loop region.",type=float)
	optional.add_argument("--chrlist",dest="ChrList",default=None,nargs='+',help="A list of space-delimited sequence names containing those sequence that should be used as validation set. If this parameter is not provided, 20%% of the data in the training set will be randomly selected as the validation set.")
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
	DripBdgFile=args.Drip
	QpoisBdgFile=args.Qpois
	Qscore=float(args.Qscore)
	TestChrList=args.ChrList
	Win=128000
	Scale=128
	RloopDeal.GetData(Genome,ChromeSize,Win,Scale,Strand,DripBdgFile,QpoisBdgFile,Qscore,Prefix,TestChrList)