#!/usr/bin/env python
import os,collections
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import ReduceLROnPlateau,ModelCheckpoint,CSVLogger,EarlyStopping
from sklearn.metrics import precision_recall_curve,average_precision_score
import numpy as np
import sys,glob,pandas
from tensorflow.keras import layers
from tensorflow.keras import backend as K
import tensorflow as tf
import scipy,gc,random
from Bio import SeqIO
from deepRloopPre import RloopDeal
def Rsquare(y_true, y_pred):
	SS_res =  K.sum(K.square( y_true-y_pred ))
	SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
	return 1 - SS_res/(SS_tot + K.epsilon())
def Predict(ModelIn,Win,Genome,ChromeSize,Strand,Prefix,Threshold,Scale,Mem="high"):
	Lines=os.popen("bedtools makewindows -g %s -w %s"%(ChromeSize,Win)).readlines()
	Wind=collections.defaultdict(list)
	for x in Lines:
		x=x.rstrip()
		l=x.split("\t")
		Wind[l[0]].append(x)
	Cl=list(Wind.keys())
	print(Cl)
	FrR=open("%s_predict.bdg"%Prefix,"w")
	FrC=open("%s_predict.bed"%Prefix,"w")
	for record in SeqIO.parse(Genome,"fasta"):
		if not record.id in Cl:
			continue
		ChromSeq=str(record.seq)
		Xpredict=[]
		for x in Wind[record.id]:
			l=x.split("\t")
			WinSeq=ChromSeq[int(l[1]):int(l[2])]
			Ss=RloopDeal.GetOnehotSeq(WinSeq,Win,Strand)
			if Mem=="low":
				r=ModelIn(np.array([Ss]),training=False)
				GetResults(r,FrR,FrC,Strand,[x],Scale)
			else:
				Xpredict.append(Ss)
		if Mem!="low":
			r=ModelIn.predict(Xpredict,batch_size=20,verbose=1)
			GetResults(r,FrR,FrC,Strand,Wind[record.id],Scale)
	FrR.close()
	FrC.close()
	if Strand=="fwd":
		st="+"
	elif Strand=="rev":
		st="-"
	os.system("bedGraphToBigWig %s_predict.bdg %s %s_predict.bw"%(Prefix,ChromeSize,Prefix))
	if Threshold>0:
		os.system("awk '$5>=%s{print $0}' %s_predict.bed|bedtools merge -i - -d 300 -c 5 -o mean|awk -F'\t' '$3-$2>=100{num+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_peak_\"num\"\t\"$4\"\t%s\"}' >%s_final_predict.bed"%(Threshold,Prefix,Prefix,st,Prefix))
def GetResults(Output,FrR,FrC,Strand,WindList,Scale):
	RegressionList=Output[0]
	ClassificationList=Output[1]
	if Strand=="fwd":
		st="+"
	elif Strand=="rev":
		st="-"
	for rl,cl,x in zip(RegressionList,ClassificationList,WindList):
		l=x.split("\t")
		lr=[ri[0] for ri in rl] #regression
		lc=[ri[0] for ri in cl] #classification
		if Strand=="rev":
			lr=lr[::-1] #right?
			lc=lc[::-1]
		t=0
		for i in range(int(l[1]),int(l[2]),Scale):
			start=i
			end=i+Scale
			if end>int(l[2]):
				end=int(l[2])
			fm_r=l[0]+"\t"+str(start)+"\t"+str(end)+"\t"+str(np.array(lr[t]))
			fm_c=l[0]+"\t"+str(start)+"\t"+str(end)+"\t.\t"+str(np.array(lc[t]))+"\t"+st
			FrR.write(fm_r+"\n")
			FrC.write(fm_c+"\n")
			t+=1
def ConvBlock(x,filters,kernel_size,pool_size,drop_out,name=None):
	if name:
		x=layers.Conv1D(filters=filters,kernel_size=kernel_size,strides=1,padding="same",use_bias=False,name="conv1d_"+name)(x)
		x=layers.BatchNormalization(axis=2,epsilon=1.001e-5,name="batch_normalization_"+name)(x)
	else:
		x=layers.Conv1D(filters=filters,kernel_size=kernel_size,strides=1,padding="same",use_bias=False)(x) #64,7
		x=layers.BatchNormalization(axis=2,epsilon=1.001e-5)(x)
	x=layers.Activation('relu')(x)
	x=layers.MaxPooling1D(pool_size=pool_size,strides=None, padding='same',data_format='channels_last')(x)
	x=layers.Dropout(drop_out)(x) #0.2
	return x
def CreateModel(Inputs,Transfer=False):
	x=ConvBlock(Inputs,256,3,2,0.01,name="1")
	x=ConvBlock(x,64,3,4,0.15,name="2")
	x=ConvBlock(x,1024,3,4,0.19,name="3")
	x=ConvBlock(x,1024,21,4,0.168,name="4")
	x=ConvBlock(x,1024,19,1,0.01,name="5")
	x=layers.Bidirectional(layers.LSTM(5,return_sequences=True,stateful=False,dropout=0.14),name="BLSTM_rloop_regression")(x)
	RloopRegression=layers.Conv1D(filters=1,kernel_size=7,strides=1,padding="same",use_bias=False,activation="relu",name='rloop_regression')(x)
	if Transfer:
		RloopClassification=layers.Bidirectional(layers.LSTM(5,return_sequences=True,stateful=False),name="BLSTM_rloop_classification")(RloopRegression)
		RloopClassification=layers.Conv1D(filters=1,kernel_size=7,strides=1,padding="same",use_bias=False,activation="sigmoid",name='rloop_classification')(RloopClassification)
	with tf.device('/cpu:0'):
		if Transfer:
			Model = tf.keras.models.Model(Inputs,outputs=[RloopRegression,RloopClassification])
		else:
			Model = tf.keras.models.Model(Inputs,RloopRegression)
	return Model
def Train(Data,Prefix,Epoch,BatchSize):
	timesteps=128000
	data_dim=4
	batch_size=BatchSize #20
	n_epoch=Epoch #1000
	X_train=Data["X_train"]
	y_rlooptrain_regression=Data["y_rlooptrain_regression"]
	y_rlooptrain_classification=Data['y_rlooptrain_classification']
	X_test=Data['X_test']
	y_rlooptest_regression=Data['y_rlooptest_regression']
	y_rlooptest_classification=Data['y_rlooptest_classification']
	y_rlooptrain_regression=np.expand_dims(np.array(y_rlooptrain_regression),axis=2)
	y_rlooptrain_classification=np.expand_dims(np.array(y_rlooptrain_classification),axis=2)
	y_rlooptest_regression=np.expand_dims(np.array(y_rlooptest_regression),axis=2)
	y_rlooptest_classification=np.expand_dims(np.array(y_rlooptest_classification),axis=2)
	Inputs=layers.Input(shape=(timesteps,data_dim))
	ModelRegression=CreateModel(Inputs)
	ModelRegression.summary()
	CheckPoint=ModelCheckpoint(filepath="%s_{epoch:02d}_{val_Rsquare:.5f}_regression.hdf5"%Prefix,monitor='val_Rsquare',verbose=1,save_best_only=True,mode='max') #
	CsvLogger=CSVLogger('%s_regression_training.log'%Prefix,separator='\t')
	EarlyStop=EarlyStopping(monitor='val_Rsquare',patience=300,mode='max')#monitor need change
	CallBacks=[CheckPoint,CsvLogger,EarlyStop]
	ModelRegression.compile(optimizer='Adam',loss={'rloop_regression':'poisson'},metrics={"rloop_regression":Rsquare})
	ModelRegression.fit(X_train,{"rloop_regression":y_rlooptrain_regression},batch_size=batch_size,epochs=n_epoch,shuffle=True,validation_data=(X_test,y_rlooptest_regression),callbacks=CallBacks)
	RegressionH5fList=[]
	for h5f in glob.glob("%s_*_*_regression.hdf5"%Prefix):
		RegressionH5fList.append((os.path.getctime(h5f),h5f))
	FinalRegressionH5=sorted(RegressionH5fList)[-1][1]
	print("Optimal regression model:%s"%FinalRegressionH5)
	ModelRegression = load_model(FinalRegressionH5,{'Rsquare':Rsquare})
	ModelClassification=CreateModel(Inputs,Transfer=True)
	for layer in ModelClassification.layers:
		#print(layer.name)
		if layer.name.startswith("conv1d") or layer.name=="BLSTM_rloop_regression" or layer.name.startswith("batch_normalization") or layer.name=="rloop_regression": #bidirectional
			#ModelClassification.get_layer(layer.name).set_weights(layer.get_weights())
			layer.set_weights(ModelRegression.get_layer(layer.name).get_weights())
			#ModelClassification.get_layer(layer.name).trainable = False
			layer.trainable = False
	ModelClassification.summary()
	CheckPoint=ModelCheckpoint(filepath="%s_{epoch:02d}_{val_rloop_classification_loss:.5f}_classification.hdf5"%Prefix,monitor='val_rloop_classification_loss',verbose=1,save_best_only=True,mode='min')
	CsvLogger=CSVLogger('%s_classification_training.log'%Prefix,separator='\t')
	EarlyStop=EarlyStopping(monitor='val_rloop_classification_loss',patience=300,mode='min')#monitor need change
	CallBacks=[CheckPoint,CsvLogger,EarlyStop]
	ModelClassification.compile(optimizer='Adam',loss={'rloop_regression':'poisson','rloop_classification':'binary_crossentropy'},metrics={"rloop_regression":Rsquare,"rloop_classification":['Precision','Recall']})
	ModelClassification.fit(X_train,{"rloop_regression":y_rlooptrain_regression,"rloop_classification":y_rlooptrain_classification},batch_size=batch_size,epochs=n_epoch,shuffle=True,validation_data=(X_test,[y_rlooptest_regression,y_rlooptest_classification]),callbacks=CallBacks)
	ClassificationH5fList=[]
	for h5f in glob.glob("%s_*_*_classification.hdf5"%Prefix):
		ClassificationH5fList.append((os.path.getctime(h5f),h5f))
	FinalClassificationH5=sorted(ClassificationH5fList)[-1][1]
	print("Optimal classification model:%s"%FinalClassificationH5)
