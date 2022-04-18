#!/usr/bin/env python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
from tensorflow.keras.models import load_model,Model
import sys,ig,generate_net,get_generate_seq
import numpy as np
def my_loss_fn(y_pred,filter_num):
	return 1-tf.reduce_mean(y_pred[:,:,filter_num])
def get_layer_model(layer_num,model,input_list):
	x=generate_net.GenerLayer()(input_list)
	for li in range(layer_num+1):
		for layer_name in ["conv1d","batch_normalization","activation","max_pooling1d","dropout"]:
			if li!=0:
				layer_name=layer_name+"_"+str(li)
			layer_item=model.get_layer(layer_name)
			layer_item.trainable=False
			if "dropout" in layer_name and li==layer_num:
				outputs=layer_item(x)
			else:
				x=layer_item(x)
	fModel=Model(inputs=input_list,outputs=outputs)
	fModel.summary()
	return fModel
def train_model(fModel,n_epoch,mask,template,filter_num):
	optimizer = tf.keras.optimizers.Adam(learning_rate=0.01)
	for i in range(n_epoch):
		with tf.GradientTape() as tape:
			y_pred=fModel([mask,template])
			#print(i,"relu",tf.reduce_mean(y_pred[:,:,filter_num]).numpy(),y_pred.shape)
			loss_value=my_loss_fn(y_pred,filter_num)
			#print(loss_value)
			#print("entropy",fModel.losses)
			loss_value +=fModel.losses[0]
			#print(loss_value)
		grads = tape.gradient(loss_value,fModel.trainable_weights)
		optimizer.apply_gradients(zip(grads,fModel.trainable_weights))
	y_pred=fModel([mask,template])
	max_act_value=tf.reduce_mean(y_pred[:,:,filter_num]).numpy()
	return max_act_value
if __name__ == '__main__':
	h5=sys.argv[1]
	prefix=sys.argv[2]
	win_length=128000
	dim=4
	n_epoch=2000
	scale_win=128
	model = load_model(h5,{'Rsquare':ig.Rsquare})
	mask=np.ones((1,win_length,dim))
	template=np.zeros((1,win_length,dim))
	input1=keras.layers.Input(shape=mask.shape[1:],name="input1")
	input2=keras.layers.Input(shape=template.shape[1:],name="input2")
	input_list=[input1,input2]
	Fr_seq=open(prefix+"_maxacseq.fasta","w")
	Fr_score=open(prefix+"_maxacseq_score.bdg","w")
	Fr_maxac=open(prefix+"_maxac.txt","w")
	#Fr_rloop=open(prefix+"_maxacseq_Rlooplevel.bdg","w")
	pwmd={}
	for layer_num in range(5): #top five layers
		if layer_num==0:
			filter_total_num=model.get_layer("conv1d").filters
		else:
			filter_total_num=model.get_layer("conv1d_%s"%layer_num).filters
		fModel=get_layer_model(layer_num,model,input_list)
		aModel=Model(inputs=fModel.input,outputs=fModel.get_layer("gener_layer").output)
		weights =fModel.get_layer("gener_layer").get_weights()
		for filter_num in range(filter_total_num):
			max_act_value=train_model(fModel,n_epoch,mask,template,filter_num)
			print("layer_num,filter_num,max_act_value:",layer_num,filter_num,max_act_value)
			pwm=aModel.predict([mask,template])
			fModel.get_layer("gener_layer").set_weights(weights)
			#print(pwm[0])
			seq,score=get_generate_seq.ont_hot_to_seq(pwm[0])
			#r=model.predict(pwm)
			seq_name=prefix+"_conv1d_%s_filter_%s"%(layer_num,filter_num)
			Fr_seq.write(">"+seq_name+"\n")
			Fr_seq.write(seq+"\n")
			for i,s in enumerate(score):
				Fr_score.write(seq_name+"\t"+str(i)+"\t"+str(i+1)+"\t"+str(s)+"\n")
			pwmd[seq_name]=pwm[0]
			Fr_maxac.write(seq_name+"_maxac"+"\t"+str(max_act_value)+"\n")
			Fr_maxac.flush()
			#for i,s in enumerate(r[0][0]):
			#	st=i*scale_win
			#	ed=st+scale_win
			#	Fr_rloop.write(prefix+"\t"+str(st)+"\t"+str(ed)+"\t"+str(s[0])+"\n")
		K.clear_session()
		#break
	Fr_seq.close()
	Fr_score.close()
	Fr_maxac.close()
	#Fr_rloop.close()
	np.savez_compressed(prefix+"_maxacseq_pwm.npz",**pwmd)	
