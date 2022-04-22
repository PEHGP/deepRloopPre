#!/usr/bin/env python
#coding=utf-8
import sys,ig,generate_net
import tensorflow as tf
from tensorflow.keras.models import load_model,Model
import numpy as np
def ont_hot_to_seq(pwm):
	seq=""
	score=[]
	for x in pwm:
		#print(x)
		p=np.argwhere(x==np.max(x))
		#print(p)
		real_p=p[0][0]
		score.append(x[real_p])
		if real_p==0:
			seq+="A"
		elif real_p==1:
			seq+="T"
		elif real_p==2:
			seq+="C"
		elif real_p==3:
			seq+="G"
	return seq,score
if __name__ == '__main__':
	h5=sys.argv[1]
	npz=sys.argv[2]
	prefix=sys.argv[3]
	scale_win=128
	data=np.load(npz)
	mask=np.expand_dims(data["mask"],axis=0)
	template=np.expand_dims(data["template"],axis=0)
	fModel = load_model(h5,custom_objects={'Rsquare':ig.Rsquare,'GenerLayer':generate_net.GenerLayer,'my_loss_fn':generate_net_v2.my_loss_fn})
	gModel=Model(inputs=fModel.input,outputs=fModel.get_layer("gener_layer").output)
	pwm=gModel.predict([mask,template])
	seq,score=ont_hot_to_seq(pwm[0])
	Fr=open(prefix+"_generateseq.fasta","w")
	Fr.write(">"+prefix+"\n")
	Fr.write(seq+"\n")
	Fr.close()
	Fr=open(prefix+"_generateseq_score.bdg","w") #每个碱基的支持率或者频率
	for i,s in enumerate(score):
		Fr.write(prefix+"\t"+str(i)+"\t"+str(i+1)+"\t"+str(s)+"\n")
	Fr.close()
	r=fModel.predict([mask,template])
	Fr=open(prefix+"_generateseq_Rlooplevel.bdg","w") #生成序列的Rloop水平
	for i,s in enumerate(r[0]):
		st=i*scale_win
		ed=st+scale_win
		Fr.write(prefix+"\t"+str(st)+"\t"+str(ed)+"\t"+str(s[0])+"\n")
	Fr.close()
	#sys.exit()
	Fr=open(prefix+"_generateseq_targetRlooplevel.bdg","w") #生成序列去拟合的目标target
	for i,s in enumerate(data["y_scale"]):
		st=i*scale_win
		ed=st+scale_win
		Fr.write(prefix+"\t"+str(st)+"\t"+str(ed)+"\t"+str(s)+"\n")
	Fr.close()
	Fr=open(prefix+"_targetseq_realRlooplevel.bdg","w") #target目标序列的真实R-loop水平
	for i,s in enumerate(data["y_real"]):
		st=i*scale_win
		ed=st+scale_win
		Fr.write(prefix+"\t"+str(st)+"\t"+str(ed)+"\t"+str(s)+"\n")
	Fr.close()
	Fr=open(prefix+"_targetseq_predRlooplevel.bdg","w") #target目标序列使用模型预测的R-loop水平
	for i,s in enumerate(data["y_pred"]):
		st=i*scale_win
		ed=st+scale_win
		Fr.write(prefix+"\t"+str(st)+"\t"+str(ed)+"\t"+str(s)+"\n")
	Fr.close()
	np.savez_compressed(prefix+"_generateseq_pwm.npz",pwm=pwm[0])
	Fr=open(prefix+"_targetseq.bed","w")
	Fr.write(prefix+"\t"+str(data["X_mask_start"])+"\t"+str(data["X_mask_end"])+"\n")
	Fr.close()
