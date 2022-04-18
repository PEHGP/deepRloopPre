#!/usr/bin/env python
import sys,os,ig,argparse
from Bio import SeqIO
import numpy as np
from tensorflow.keras.models import load_model,Model
def get_parser():
	parser = argparse.ArgumentParser(description="ig_seq.py",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--fasta',dest="fasta")
	parser.add_argument('--win',dest="win",type=int) #choices=["fwd","rev"]
	parser.add_argument('--chrome_size',dest="chrome_size")
	parser.add_argument('--h5',dest="h5")#help=""
	parser.add_argument('--select',dest="select",choices=["random","shuffle"],default="random")
	parser.add_argument("--num_steps",dest="num_steps",type=int,default="50")
	parser.add_argument("--num_runs",dest="num_runs",type=int,default="2")
	parser.add_argument("--prefix",dest="prefix")
	return parser
if __name__ == '__main__':
	parser=get_parser()
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit()
	else:
		Args = parser.parse_args()
	fasta=Args.fasta
	win=Args.win
	chrome_size=Args.chrome_size
	h5=Args.h5
	select=Args.select
	num_steps=Args.num_steps
	num_runs=Args.num_runs
	prefix=Args.prefix
	Fr_log=open(prefix+".log","w")
	Fr_log.write(" ".join(sys.argv)+"\n")
	Fr_bdg=open(prefix+"_imscore.bdg","w")
	one_hot={"A":[1,0,0,0],"T":[0,1,0,0],"C":[0,0,1,0],"G":[0,0,0,1],"N":[0,0,0,0]}
	dseq=SeqIO.index(fasta,"fasta")
	model = load_model(h5,{'Rsquare':ig.Rsquare})
	inputs=model.layers[0].input
	outputs=model.get_layer("rloop_regression").output
	fModel=Model(inputs=inputs,outputs=outputs)
	lines=os.popen("bedtools makewindows -g %s -w %s "%(chrome_size,win)).readlines()
	igrads_matrix=[] #for heatmap
	imscore_matrix=[] #for seqlogo
	for x in lines:
		x=x.rstrip()
		l=x.split("\t")
		seq=dseq[l[0]][int(l[1]):int(l[2])]
		ss=[]
		for s in seq.upper():
			if s in "ATCG":
				ss.append(one_hot[s])
			else:
				ss.append(one_hot["N"])
		if len(ss)<int(win):
			for i in range(int(win)-len(ss)):
				ss.append(one_hot["N"])
		X_test=np.array([ss])
		print("X_test shape",X_test.shape)
		igrads=ig.baseline_integrated_gradients(X_test,fModel,num_steps=num_steps,num_runs=num_runs,select=select) #num_step,num_run need change
		#print(igrads) #(1, 128000, 4)
		#igrads_matrix.append(igrads[0])
		score=np.sum(np.abs(igrads),axis=2) #right?
		print("score shape",score.shape) #(1, 128000)
		#print(score)
		score_final=X_test*score[:,:,None]
		print("score_final shape",score_final.shape) #(1, 128000, 4)
		#print(score_final)
		#imscore_matrix.append(score_final)
		for ni,pi,gi,si in zip(np.around(np.reshape(np.sum(score_final,axis=2),-1),4),range(int(l[1]),int(l[2])),igrads[0],score_final[0]):
			fm=l[0]+"\t"+str(pi)+"\t"+str(pi+1)+"\t"+str(ni)
			Fr_bdg.write(fm+"\n")
			igrads_matrix.append(gi)
			imscore_matrix.append(si)
	Fr_bdg.close()
	Fr_log.close()
	np.savez_compressed("%s_matrix.npz"%prefix,igrads_matrix=np.around(igrads_matrix,4),imscore_matrix=np.around(imscore_matrix,4))
	os.system("bedGraphToBigWig %s_imscore.bdg %s %s_imscore.bw"%(prefix,chrome_size,prefix))