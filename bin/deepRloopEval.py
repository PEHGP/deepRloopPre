#!/usr/bin/env python
from deepRloopPre._version import __version__
import argparse,time,sys
import numpy as np
from deepRloopPre import RloopDeal
def GetParser():
	parser = argparse.ArgumentParser(description="""

Comprehensive evaluation of prediction results.
A detailed sub-commands help is available by typing:
	deepRloopEval.py getvalue -h
	deepRloopEval.py plotpr -h
	deepRloopEval.py plotpeak -h
	deepRloopEval.py plotgene -h
	deepRloopEval.py plotcorr -h
""",
		formatter_class=argparse.RawDescriptionHelpFormatter,epilog="""
example usage:
	deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak predicted_peaks.bed --prefix test
	deepRloopEval.py plotpr --target target --prefix test
	deepRloopEval.py plotpeak --truepeak observed_peaks.bed --predpeak predicted_all.bed --truebw true.bw --predbw pred.bw --prefix test
	deepRloopEval.py plotgene --genebed gene.bed --truefwdbw truefwd.bw --truerevbw truerev.bw --predfwdbw predfwd.bw --predrevbw predrev.bw --prefix test
	deepRloopEval.py plotcorr --genebed gene.bed --truefwdbw truefwd.bw --truerevbw truerev.bw --predfwdbw predfwd.bw --predrevbw predrev.bw --genebody --tssextend 300 --ttsextend 300 --prefix test
""",
		conflict_handler='resolve',)
	parser.add_argument('--version', action='version',version='deepRloopPre {}'.format(__version__))
	subparsers = parser.add_subparsers(title="commands",dest='Command',metavar='')
	subparsers.add_parser("getvalue",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("getvalue")],help="getting the F1-score, Precision, Recall, Jaccard, AP, Spearman.",add_help=False,usage="deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak predicted_peaks.bed --prefix test")
	subparsers.add_parser("plotpr",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("plotpr")],help="plot precision-recall curve for many samples.",add_help=False,usage="deepRloopEval.py plotpr --target target --prefix test")
	subparsers.add_parser("plotpeak",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("plotpeak")],help="plot R-loop peak region metaplot.",add_help=False,usage="deepRloopEval.py plotpeak --truepeak observed_peaks.bed --predpeak predicted_all.bed --truebw true.bw --predbw pred.bw --prefix test")
	subparsers.add_parser("plotgene",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("plotgene")],help="plot sense/antisense metaplot or region metaplot.",add_help=False,usage="deepRloopEval.py plotgene --genebed gene.bed --truefwdbw truefwd.bw --truerevbw truerev.bw --predfwdbw predfwd.bw --predrevbw predrev.bw --prefix test")
	subparsers.add_parser("plotcorr",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("plotcorr")],help="plot observed and predicted sense/antisense correlation.",add_help=False,usage="deepRloopEval.py plotcorr --genebed gene.bed --truefwdbw truefwd.bw --truerevbw truerev.bw --predfwdbw predfwd.bw --predrevbw predrev.bw --genebody --tssextend 300 --ttsextend 300 --prefix test")
	return parser
def GetCaseArgs(Case):
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument("--help", "-h", action="help",help="show this help message and exit")
	parser.add_argument("--prefix",dest="Prefix",required=True)
	if Case=="plotpr":
		#parser.add_argument('--truepeak',dest="TruePeak",help='observed R-loop location in bed format. Only the first three columns in bed format are required.',required=True)
		#parser.add_argument("--predpeak",dest="PredPeak",default=None,help="The bed file contains all the predicted scores of all sequences, not just the scores of the R-loop region. The bed file needs to provide the first six columns of information, of which the fifth column is the score. If the target file is provided, this parameter is not required.",required=True)
		#parser.add_argument("--havescale",dest="HaveScale",default=1,choices=[0,1],help="whether the predpeak bed file in the second column has been divided according to the 128bp window. If it has been divided and selected 1, there is no division selection 0.")
		parser.add_argument("--target",dest="Target",help="The target file contains six columns. The first column is the sample name. The second column is a bed file containing the observed R-loop location. The third column is a bed file containing the prediction scores of all sequences. The fourth column indicates whether the bed file in the second column has been divided according to the 128bp window. If it has been divided and selected 1, there is no division selection 0. The fifth column is chrome size file.The sixth column is color, using hexadecimal color code.",required=True)
		#parser.add_argument('--chromsize',dest="ChromeSize",help='chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.',required=True)
	elif Case=="getvalue":
		parser.add_argument('--truepeak',dest="TruePeak",default=None,help='observed R-loop location in bed format. Only the first three columns in bed format are required.')
		parser.add_argument("--predpeak",dest="PredPeak",default=None,help="predicted R-loop location in bed format. Only the first three columns in bed format are required.")
		parser.add_argument("--truebw",dest="TrueBw",default=None,help="observed R-loop abundance in bigwig format.")
		parser.add_argument("--predbw",dest="PredBw",default=None,help="predicted R-loop abundance in bigwig format.")
		parser.add_argument('--chromsize',dest="ChromeSize",default=None,help='chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.')
		parser.add_argument("--havescale",dest="HaveScale",default="1",choices=["0","1"],help="whether the predpeak bed file in the second column has been divided according to the 128bp window. If it has been divided and selected 1, there is no division selection 0.")
		parser.add_argument("--predall",dest="PredAllPeak",default=None,help="The bed file contains all the predicted scores of all sequences, not just the scores of the R-loop region. The bed file needs to provide the first six columns of information, of which the fifth column is the score.")
	elif Case=="plotpeak":
		parser.add_argument('--truepeak',dest="TruePeak",help='observed R-loop location in bed format. Only the first three columns in bed format are required.',required=True)
		parser.add_argument("--predpeak",dest="PredPeak",help="predicted R-loop location in bed format. Only the first three columns in bed format are required.",required=True)
		parser.add_argument("--truebw",dest="TrueBw",help="observed R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--predbw",dest="PredBw",help="predicted R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--extend",dest="Extend",default=1000,help="distance upstream of the start site of the regions and distance downstream of the end site of the regions.",type=int)
		parser.add_argument("--thread",dest="Thread",default=1,help="number of processors to use.",type=int)
	elif Case=="plotgene":
		parser.add_argument('--genebed',dest="GeneBed",help='gene location in bed format. If the bed file provides six columns, plot sense/antisense metaplot. If the bed file provides three columns, just plot gene region metaplot',required=True)
		parser.add_argument("--truefwdbw",dest="TrueFwdBw",help="observed Watson R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--truerevbw",dest="TrueRevBw",help="observed Crick R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--predfwdbw",dest="PredFwdBw",help="predicted Watson R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--predrevbw",dest="PredRevBw",help="predicted Crick R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--extend",dest="Extend",default=1000,help="distance upstream of the start site of the regions and distance downstream of the end site of the regions.",type=int)
		parser.add_argument("--thread",dest="Thread",default=1,help="number of processors to use.",type=int)
	elif Case=="plotcorr":
		parser.add_argument('--genebed',dest="GeneBed",help='gene location in bed format. Only the first six columns in bed format are required.',required=True)
		parser.add_argument("--truefwdbw",dest="TrueFwdBw",help="observed Watson R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--truerevbw",dest="TrueRevBw",help="observed Crick R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--predfwdbw",dest="PredFwdBw",help="predicted Watson R-loop abundance in bigwig format.",required=True)
		parser.add_argument("--predrevbw",dest="PredRevBw",help="predicted Crick R-loop abundance in bigwig format.",required=True)
		parser.add_argument('--chromsize',dest="ChromeSize",help='chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.',required=True)
		parser.add_argument("--genebody",dest="GeneBody",action='store_true',help="plot gene body sense/antisense correlation.")
		parser.add_argument("--tssextend",dest="TssExtend",default=None,help="distance upstream and downstream of the start site of the regions.",type=int)
		parser.add_argument("--ttsextend",dest="TtsExtend",default=None,help="distance upstream and downstream of the end site of the regions.",type=int)
		parser.add_argument("--thread",dest="Thread",default=1,help="number of processors to use.",type=int)
	return parser
if __name__ == '__main__':
	args=GetParser().parse_args()
	if not args.Command:
		GetParser().print_help()
		sys.exit()
	Scale=128 #as parameter?
	#print(args.Command)
	if args.Command=="plotpr":
		Prefix=args.Prefix
		#TruePeak=args.TruePeak
		#ChromeSize=args.ChromeSize
		#HaveScale=int(args.HaveScale)
		#PredictAllBed=args.PredPeak
		Target=args.Target
		#if HaveScale==1:
		#	HaveScale=True
		#else:
		#	HaveScale=False
		RloopDeal.PlotPrList(Scale,Prefix,Target)
	elif args.Command=="getvalue":
		Prefix=args.Prefix
		TruePeak=args.TruePeak
		PredPeak=args.PredPeak
		TrueBw=args.TrueBw
		PredBw=args.PredBw
		PredAllBed=args.PredAllPeak
		ChromeSize=args.ChromeSize
		HaveScale=int(args.HaveScale)
		if HaveScale==1:
			HaveScale=True
		else:
			HaveScale=False
		RloopDeal.Evaluation(Prefix,HaveScale,PredAllBed,ChromeSize,TruePeak,PredPeak,TrueBw,PredBw)
	elif args.Command=="plotpeak":
		Prefix=args.Prefix
		TruePeak=args.TruePeak
		PredPeak=args.PredPeak
		TrueBw=args.TrueBw
		PredBw=args.PredBw
		Extend=args.Extend
		Thread=args.Thread
		RloopDeal.PlotDist(Prefix,TruePeak,PredPeak,TrueBw,PredBw,Extend,Thread)
	elif args.Command=="plotgene":
		Prefix=args.Prefix
		GeneBed=args.GeneBed
		TrueFwdBw=args.TrueFwdBw
		TrueRevBw=args.TrueRevBw
		PredFwdBw=args.PredFwdBw
		PredRevBw=args.PredRevBw
		Extend=args.Extend
		Thread=args.Thread
		RloopDeal.PlotGeneDist(Prefix,GeneBed,TrueFwdBw,TrueRevBw,PredFwdBw,PredRevBw,Extend,Thread)
	elif args.Command=="plotcorr":
		Prefix=args.Prefix
		GeneBed=args.GeneBed
		TrueFwdBw=args.TrueFwdBw
		TrueRevBw=args.TrueRevBw
		PredFwdBw=args.PredFwdBw
		PredRevBw=args.PredRevBw
		ChromeSize=args.ChromeSize
		GeneBody=args.GeneBody
		TssExtend=args.TssExtend
		TtsExtend=args.TtsExtend
		Thread=args.Thread
		RloopDeal.PlotRloopCorr(Prefix,GeneBed,TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,ChromeSize,GeneBody,TssExtend,TtsExtend,Thread)
