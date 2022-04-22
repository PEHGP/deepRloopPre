#!/usr/bin/env python
#coding=utf-8
import numpy as np
from Bio import SeqIO
import sys,collections
def get_one_hot_format(seq,one_hot):
	ss=[]
	for s in seq.upper():
		if s in "ATCG":
			ss.append(one_hot[s])
		else:
			ss.append(one_hot["N"])
	return np.array(ss)
def get_bdg_d(bdg):
	d=collections.defaultdict(dict)
	for x in open(bdg):
		x=x.rstrip()
		l=x.split("\t")
		for i in range(int(l[1]),int(l[2])):
			d[l[0]][i]=float(l[-1])
	return d
def get_param(region,genome_fasta_file,win_length,position,prefix,extend): #region=Chr3:14191758-14205372，根据region得到mask序列，然后把序列左右延申得到128000bp的模型输入序列，mask就是那一段不变的序列
	Chr,po=region.split(":")
	Start,End=po.split("-")
	Start=int(Start)
	End=int(End)
	dseq=SeqIO.index(genome_fasta_file,"fasta")
	seq_length=len(dseq[Chr].seq)
	mask_seq=str(dseq[Chr].seq[Start:End])
	left=Start-extend
	right=End+extend
	if position=="center":
		flag=0
		while End-Start!=win_length:
			if flag%2==0:
				Start=Start-1
				if Start<0:
					Start=0
			else:
				End=End+1
				if End>seq_length:
					End=seq_length
			flag+=1
	elif position=="left":
		End=End+win_length-(End-Start)
		if End>seq_length:
			End=seq_length
			Start=Start-(win_length-(End-Start))
	elif position=="right":
		Start=Start-(win_length-(End-Start))
		if Start<0:
			Start=0
			End=win_length
	Fr=open(prefix+"_realseq.fasta","w") #mask及周边真实序列
	#Fr.write(">%s:%s-%s\n"%(Chr,Start,End))
	Fr.write(">%s\n"%(prefix))
	Fr.write(str(dseq[Chr].seq[Start:End])+"\n")
	Fr.close()
	mask_seq_left=str(dseq[Chr].seq[Start:left])
	mask_seq_right=str(dseq[Chr].seq[right:End])
	return mask_seq,Start,End,mask_seq_left,mask_seq_right
if __name__ == '__main__':
	mask_region=sys.argv[1] #Chr3:100-120,rloop peak区
	genome_fasta_file=sys.argv[2] #tair10_Chr3.fasta
	position=sys.argv[3] #left right center
	drip_bdg=sys.argv[4]
	pred_bdg=sys.argv[5] #模型预测的结果
	extend=int(sys.argv[6])
	scale=float(sys.argv[6]) #up or down rloop level
	prefix=sys.argv[7]
	one_hot={"A":[1,0,0,0],"T":[0,1,0,0],"C":[0,0,1,0],"G":[0,0,0,1],"N":[0,0,0,0]}
	dim=4
	scale_win=128 #need change
	win_length=128000 #need change
	mask_seq,win_start,win_end,mask_seq_left,mask_seq_right=get_param(mask_region,genome_fasta_file,win_length,position,prefix,extend)
	ch,po=mask_region.split(":")
	mask_start,mask_end=[int(i) for i in po.split("-")]
	mask_rela_start=mask_start-win_start
	mask_rela_end=mask_end-win_start
	mask_rela_left=mask_rela_start-extend
	mask_rela_right=mask_rela_end+extend
	print(win_start,win_end)
	print(mask_rela_start,mask_rela_end)
	mask=np.ones((win_length,dim))
	mask[mask_rela_start:mask_rela_end,:]=[0,0,0,0]
	mask[0:mask_rela_left,:]=[0,0,0,0]
	mask[mask_rela_right:,:]=[0,0,0,0]
	one_hot_seq=get_one_hot_format(mask_seq,one_hot)
	one_hot_seq_left=get_one_hot_format(mask_seq_left,one_hot)
	one_hot_seq_right=get_one_hot_format(mask_seq_right,one_hot)
	template=np.zeros((win_length,dim))
	template[mask_rela_start:mask_rela_end,:]=one_hot_seq
	template[0:mask_rela_left,:]=one_hot_seq_left
	template[mask_rela_right:,:]=one_hot_seq_right
	#Fr=open(prefix+"_rela.bed","w")
	#Fr.write(prefix+"\t"+str(mask_rela_start)+"\t"+str(mask_rela_end)+"\n")
	#Fr.close()
	drip_bdgd=get_bdg_d(drip_bdg)
	pred_bdgd=get_bdg_d(pred_bdg)
	y_drip=[]
	y_pred=[]
	for flag_t,i in enumerate(range(win_start,win_end,scale_win)):
		r_start=i
		r_end=i+scale_win
		if r_end>win_end:
			r_end=win_end
		#print(r_start,r_end)
		if mask_start >=r_start and mask_start<=r_end:
			y_mask_rela_start=flag_t
		if mask_end >=r_start and mask_end<=r_end:
			y_mask_rela_end=flag_t
		t_drip_sum=0.
		t_pred_sum=0.
		for p in range(r_start,r_end):
			t_drip_sum+=drip_bdgd[ch][p]
			t_pred_sum+=pred_bdgd[ch][p]
		t_drip_mean=t_drip_sum/float(scale_win)
		t_pred_mean=t_pred_sum/float(scale_win)
		y_drip.append(t_drip_mean)
		y_pred.append(t_pred_mean)
	y_drip=np.array(y_drip)
	y_pred=np.array(y_pred)
	y_scale=np.zeros(len(y_pred)) #除了拟合区，其他区域值都为0
	print(y_mask_rela_start,y_mask_rela_end)
	y_scale[y_mask_rela_start:y_mask_rela_end]=y_pred[y_mask_rela_start:y_mask_rela_end]*scale #pred,使用预测的结果作为拟合目标，因为预测的结果本身与实际结果丰度高低有差距，去拟合实际丰度用于生成序列不合理。
	np.savez_compressed("%s_yscale%s.npz"%(prefix,scale),mask=mask,template=template,y_scale=y_scale,y_real=y_drip,y_pred=y_pred,y_mask_start=y_mask_rela_start,y_mask_end=y_mask_rela_end,X_mask_start=mask_rela_start,X_mask_end=mask_rela_end)