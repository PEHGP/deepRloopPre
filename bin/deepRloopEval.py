#!/usr/bin/env python
from deepRloopPre._version import __version__
import argparse,time,sys
import numpy as np
from deepRloopPre import RloopDeal
def GetParser():
	parser = argparse.ArgumentParser(description="""

Evaluation of prediction results.This program can output F1-score, Precision, Recall for R-loop location and Spearman for R-loop abundance or draw a precision-recall curve.
A detailed sub-commands help is available by typing:
	deepRloopEval.py getvalue -h
	deepRloopEval.py plotpr -h
""",
		formatter_class=argparse.RawDescriptionHelpFormatter,epilog="""
example usage:
	deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak predicted_peaks.bed --prefix test
	deepRloopEval.py plotpr --truepeak observed_peaks.bed --predpeak predicted_all.bed --prefix test
""",
		conflict_handler='resolve',)
	parser.add_argument('--version', action='version',version='deepRloopPre {}'.format(__version__))
	subparsers = parser.add_subparsers(title="commands",dest='Command',metavar='')
	subparsers.add_parser("getvalue",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("getvalue")],help="getting the F1-score, Precision, Recall, Spearman.",add_help=False,usage="deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak predicted_peaks.bed --prefix test")
	subparsers.add_parser("plotpr",formatter_class=argparse.ArgumentDefaultsHelpFormatter,parents=[GetCaseArgs("plotpr")],help="plot precision-recall curve.",add_help=False,usage="deepRloopEval.py plotpr --truepeak observed_peaks.bed --predpeak predicted_all.bed --prefix test")
	return parser
def GetCaseArgs(Case):
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument("--prefix",dest="Prefix",required=True)
	if Case=="plotpr":
		parser.add_argument("--help", "-h", action="help",help="show this help message and exit")
		parser.add_argument('--truepeak',dest="TruePeak",default=None,help='observed R-loop location in bed format. Only the first three columns in bed format are required.')
		parser.add_argument("--predpeak",dest="PredPeak",default=None,help="The bed file contains all the predicted scores of all sequences, not just the scores of the R-loop region. The bed file needs to provide the first six columns of information, of which the fifth column is the score. If the target file is provided, this parameter is not required.")
		parser.add_argument("--havescale",dest="HaveScale",default=1,choices=["0","1"],help="whether the predpeak bed file in the second column has been divided according to the 128bp window. If it has been divided and selected 1, there is no division selection 0.")
		parser.add_argument("--target",dest="Target",default=None,help="If the figure contains multiple samples, a target file needs to be provided. The target file contains four columns. The first column is the sample name, the second column is a bed file containing the prediction scores of all sequences, and the third column indicates whether the bed file in the second column has been divided according to the 128bp window. If it has been divided and selected 1, there is no division selection 0. The fourth column is color, using hexadecimal color code.")
		parser.add_argument('--chromsize',dest="ChromeSize",default=None,help='chrom_size file has two columns separated by tabs. The first column is the sequence name and the second column is the sequence length.',required=True)
	else:
		parser.add_argument("--help", "-h", action="help",help="show this help message and exit")
		parser.add_argument('--truepeak',dest="TruePeak",default=None,help='observed R-loop location in bed format. Only the first three columns in bed format are required.')
		parser.add_argument("--predpeak",dest="PredPeak",default=None,help="predicted R-loop location in bed format. Only the first three columns in bed format are required.")
		parser.add_argument("--truebw",dest="TrueBw",default=None,help="observed R-loop abundance in bigwig format.")
		parser.add_argument("--predbw",dest="PredBw",default=None,help="predicted R-loop abundance in bigwig format.")
	return parser
if __name__ == '__main__':
	args=GetParser().parse_args()
	if not args.Command:
		GetParser().print_help()
		sys.exit()
	Scale=128
	print(args.Command)
	if args.Command=="plotpr":
		Prefix=args.Prefix
		TruePeak=args.TruePeak
		ChromeSize=args.ChromeSize
		HaveScale=int(args.HaveScale)
		PredictAllBed=args.PredPeak
		Target=args.Target
		RloopDeal.PlotPrList(TruePeak,ChromeSize,Scale,Prefix,HaveScale,PredictAllBed,Target)
	else:
		Prefix=args.Prefix
		TruePeak=args.TruePeak
		PredictPeak=args.PredPeak
		TrueBw=args.TrueBw
		PredictBw=args.PredBw
		RloopDeal.Evaluation(Prefix,TruePeak,PredictPeak,TrueBw,PredictBw)