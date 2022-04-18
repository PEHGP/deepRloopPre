#!/usr/bin/env python
#coding=utf-8
import sys
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.applications import xception
from tensorflow.keras.models import load_model,Model
from tensorflow.keras import backend as K
def Rsquare(y_true, y_pred):
	SS_res=K.sum(K.square( y_true-y_pred ))
	SS_tot=K.sum(K.square( y_true - K.mean(y_true) ) )
	return 1 - SS_res/(SS_tot + K.epsilon())
def get_gradients(X_test,model):
	X_test=tf.convert_to_tensor(X_test,dtype=tf.float32)
	with tf.GradientTape() as tape:
		tape.watch(X_test)
		preds = model(X_test)
	grads = tape.gradient(preds,X_test)
	return grads
def get_integrated_gradients(X_test, model, baseline=None, num_steps=50):
	if baseline is None:
		baseline = np.zeros(X_test.shape).astype(np.float32)
	else:
		baseline = baseline.astype(np.float32)

	# 1. Do interpolation.
	X_test = X_test.astype(np.float32)
	interpolated = [baseline + (step / num_steps) * (X_test - baseline) for step in range(num_steps + 1)]
	interpolated = np.array(interpolated).astype(np.float32)

	# 3. Get the gradients
	grads = []
	for i,ipd in enumerate(interpolated):
		#print("ipd shape",ipd.shape) #ipd shape (1, 128000, 4)
		grad = get_gradients(ipd,model)
		#print(grad)
		grads.append(grad[0])
	grads = tf.convert_to_tensor(grads, dtype=tf.float32)

	# 4. Approximate the integral using the trapezoidal rule
	grads = (grads[:-1] + grads[1:]) / 2.0
	avg_grads = tf.reduce_mean(grads, axis=0)

	# 5. Calculate integrated gradients and return
	integrated_grads = (X_test - baseline) * avg_grads
	return integrated_grads

def baseline_integrated_gradients(X_test, model,num_steps=50, num_runs=2,select="random"):
	# 1. List to keep track of Integrated Gradients (IG) for all the images
	integrated_grads = []

	# 2. Get the integrated gradients for all the baselines
	for run in range(num_runs):
		if select=="random":
			#np.random.seed(10)
			baseline = np.random.random(X_test.shape)
			total=np.sum(baseline,axis=2) #(1,12800),total[:,:,None].shape==(1,12800,1)
			baseline=baseline/total[:,:,None] #(1,12800,4)
			#print(np.sum(baseline,axis=2))
		elif select=="shuffle":
			baseline=X_test.copy()
			list(map(np.random.shuffle, baseline)) #Multi-dimensional arrays are only shuffled along the first axis，直接shuffle原始序列，不加list，结果不shuffle
			#print("baseline shape",baseline.shape) #(1, 128000, 4)
		#baseline=None
		igrads = get_integrated_gradients(X_test,model,baseline=baseline,num_steps=num_steps)
		integrated_grads.append(igrads)
	# 3. Return the average integrated gradients for the image
	integrated_grads = tf.convert_to_tensor(integrated_grads)
	#print("integrated_grads shape",integrated_grads.shape) # (num_runs, 1, 128000, 4)
	return tf.reduce_mean(integrated_grads, axis=0)
if __name__ == '__main__':
	h5=sys.argv[1]
	npz=sys.argv[2]
	prefix=sys.argv[3]
	data=np.load(npz)
	X_test=data["X_test"][0:10,:,:]
	model = load_model(h5,{'Rsquare':Rsquare})
	inputs=model.layers[0].input
	outputs=model.get_layer("rloop_regression").output
	fModel=Model(inputs=inputs,outputs=outputs)
	Fr_step=open("%s_step.xls"%prefix,"w")
	for xi in X_test:
		xi=xi[None,:,:]
		print("xi shape",xi.shape)
		igrads=baseline_integrated_gradients(xi,fModel,num_steps=50,num_runs=10,select="random")
		#print(igrads) #(1, 128000, 4)
		score=np.sum(np.abs(igrads),axis=2) #right?
		print("score shape",score.shape) #(1, 128000)
		#print(score)
		score_final=xi*score[:,:,None]
		print("score_final shape",score_final.shape) #(1, 128000, 4)
		#print(score_final)
		for ni in np.reshape(np.sum(score_final,axis=2),-1):
			Fr_step.write(str(ni)+"\n")
	Fr_step.close()
	"""
	Fr_batch=open("%s_batch.xls"%prefix,"w")
	igrads=baseline_integrated_gradients(X_test,fModel,num_steps=50,num_runs=10,select="random")
	score=np.sum(np.abs(igrads),axis=2)
	print("score shape",score.shape)
	score_final=X_test*score[:,:,None]
	print("score_final shape",score_final.shape)
	for timesteps in np.sum(score_final,axis=2):
		for ni in timesteps:
			Fr_batch.write(str(ni)+"\n")
	Fr_batch.close()
	"""
