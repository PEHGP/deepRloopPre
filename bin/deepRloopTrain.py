#!/usr/bin/env python
import argparse,time
import numpy as np
from deepRloopPre._version import __version__
def GetParser():
	RequiredArgs = GetRequiredArgs()
	OptionalArgs = GetOptionalArgs()
	parser = argparse.ArgumentParser(description="Building your own deepRloopPre model.",formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='$ deepRloopTrain.py --npz dataset.npz --prefix test',parents=[RequiredArgs,OptionalArgs],add_help=False)
	return parser
def GetRequiredArgs():
	parser = argparse.ArgumentParser(add_help=False)
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--npz',dest="Npz",required=True,help="Training set in npz format.")
	required.add_argument("--prefix",dest="Prefix",required=True)
	return parser
def GetOptionalArgs():
	parser = argparse.ArgumentParser(add_help=False)
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument("--help", "-h", action="help",help="show this help message and exit")
	optional.add_argument('--epoch',dest="Epoch",default=1000,help="Epoch",type=int)
	optional.add_argument("--bs",dest="BatchSize",default=20,help="batch size",type=int)
	optional.add_argument('--version', action='version',version='deepRloopPre {}'.format(__version__))
	return parser
if __name__ == '__main__':
	args=GetParser().parse_args()
	Data=np.load(args.Npz)
	Prefix=args.Prefix
	Epoch=int(args.Epoch)
	BatchSize=int(args.BatchSize)
	from deepRloopPre import RloopModel
	RloopModel.Train(Data,Prefix,Epoch,BatchSize)