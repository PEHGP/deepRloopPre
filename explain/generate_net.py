#!/usr/bin/env python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
from tensorflow.keras.models import load_model,Model
from Bio import SeqIO
import sys,ig
import numpy as np
from tensorflow.keras.callbacks import ModelCheckpoint,CSVLogger,EarlyStopping
class GenerLayer(keras.layers.Layer):
	def __init__(self,**kwargs):
		super(GenerLayer, self).__init__(**kwargs)
	def build(self, input_shape):
		#assert isinstance(input_shape, list)
		#print(input_shape)
		self.w = self.add_weight(name='weight',shape=input_shape[0][1:],initializer='glorot_uniform',trainable=True)
		super(GenerLayer, self).build(input_shape)
	def call(self, x):
		assert isinstance(x, list)
		mask,template=x
		target_bits=2.0
		norm_weights=keras.layers.Softmax()(self.w)
		#print(norm_weights)
		target_entropy_mse_loss=target_entropy_mse(norm_weights,target_bits)
		self.add_loss(target_entropy_mse_loss) #https://keras.io/api/losses/
		return tf.math.multiply(norm_weights,mask)+template #add(+) is right?
	#def compute_output_shape(self, input_shape):
	#	return input_shape[0]
def target_entropy_mse(pwm,target_bits): #seqprop_optimizer.py
	entropy = pwm * -tf.math.log(tf.clip_by_value(pwm, K.epsilon(), 1. - K.epsilon()))/tf.math.log(2.0)
	#print(entropy)
	entropy = tf.reduce_sum(entropy,axis=1)
	conservation = 2.0 - entropy
	target_entropy_mse_loss=tf.reduce_mean((conservation - target_bits)**2, axis=-1)
	#target_entropy_sme_loss=(tf.reduce_mean(conservation, axis=-1) - target_bits)**2
	return target_entropy_mse_loss
def my_loss_fn(y_true, y_pred):
	#loss=keras.losses.poisson(y_true[:,y_mask_start:y_mask_end,:], y_pred[:,y_mask_start:y_mask_end,:])
	mse=keras.losses.MeanSquaredError()
	loss=mse(y_true[:,y_mask_start:y_mask_end,:], y_pred[:,y_mask_start:y_mask_end,:])
	return tf.reduce_mean(loss) 
if __name__ == '__main__':
	h5=sys.argv[1]
	npz=sys.argv[2]
	prefix=sys.argv[3]
	model = load_model(h5,{'Rsquare':ig.Rsquare})
	data=np.load(npz)
	mask=np.expand_dims(data["mask"],axis=0)
	template=np.expand_dims(data["template"],axis=0)
	y_scale=np.expand_dims(data["y_scale"],axis=(0,2))
	global y_mask_start,y_mask_end
	y_mask_start=data["y_mask_start"]
	y_mask_end=data["y_mask_end"]
	input1=keras.layers.Input(shape=mask.shape[1:],name="input1")
	input2=keras.layers.Input(shape=template.shape[1:],name="input2")
	x=GenerLayer()([input1,input2])
	for layer_item in model.layers:
		print(layer_item.name)
		if layer_item.name=="input_1":
			continue
		layer_item.trainable=False
		if layer_item.name=="rloop_regression":
			outputs=layer_item(x)
			break
		else:
			x=layer_item(x)
	fModel=Model(inputs=[input1,input2],outputs=outputs)
	fModel.summary()
	#sys.exit()
	#checkpoint=ModelCheckpoint(filepath="generate_model_trained/%s_{epoch:02d}_{loss:.4f}.hdf5"%prefix,monitor='loss',verbose=1,save_best_only=True,mode='min')
	checkpoint=ModelCheckpoint(filepath="generate_model_trained/%s.hdf5"%prefix,monitor='loss',verbose=1,save_best_only=True,mode='min')
	csv_logger=CSVLogger('%s_training.log'%prefix,separator='\t')
	callbacks=[csv_logger,checkpoint]
	fModel.compile(optimizer="Adam",loss={'rloop_regression':my_loss_fn},metrics={"rloop_regression":ig.Rsquare})
	fModel.fit([mask,template],y_scale,epochs=10000,callbacks=callbacks)
	#for layer_item in fModel.layers:
	#	print(layer_item.name)
	gModel=Model(inputs=fModel.input,outputs=fModel.get_layer("gener_layer").output)
	pwm=gModel.predict([mask,template])
	print(pwm)
