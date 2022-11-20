#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import os,collections
from sklearn.metrics import precision_recall_curve,average_precision_score
import numpy as np
import sys,glob,pandas
import scipy,gc,random
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn import datasets, linear_model
def GetOnehot():
	Onehot={"A":[1,0,0,0],"T":[0,1,0,0],"C":[0,0,1,0],"G":[0,0,0,1],"N":[0,0,0,0]}
	return Onehot
def GetOnehotSeq(Sequence,Win,Strand):
	Onehot=GetOnehot()
	if Strand=="rev":
		Sequence=str(Seq(Sequence).reverse_complement())
	Ss=[]
	for s in Sequence.upper():
		if s in "ATCG":
			Ss.append(Onehot[s])
		else:
			Ss.append(Onehot["N"])
	if len(Sequence)<int(Win):
		if Strand=="rev":
			for i in range(int(Win)-len(Ss)):
				Ss.insert(0,Onehot["N"])
		else:
			for i in range(int(Win)-len(Ss)):
				Ss.append(Onehot["N"])
	return Ss

def GetBdgd(bdg): #bdg must sort
	d=collections.defaultdict(list)
	for x in open(bdg):
		x=x.rstrip()
		l=x.split("\t")
		for i in range(int(l[1]),int(l[2])):
			d[l[0]].append(float(l[-1])) #good
	return d
def GetySet(winlines,win,scale_win,drip_d,qpois_d,test_winlines,strand,qscore):
	y_train_regression=[]
	y_train_classification=[]
	y_test_regression=[]
	y_test_classification=[]
	for x in winlines:
		x=x.rstrip()
		l=x.split("\t")
		sc_regression=[]
		sc_classification=[]
		for i in range(int(l[1]),int(l[2]),scale_win):
			start=i
			end=i+scale_win
			if end>int(l[2]):
				end=int(l[2])
			t_sum=0
			for p in range(start,end):
				try:
					t_sum+=drip_d[l[0]][p]
				except:
					print("Waring:%s is not in drip.bdg"%l[0],l[0],p)
					continue
			t_mean=t_sum/float(scale_win)
			for p in range(start,end): 
				try:
					if qpois_d[l[0]][p]>=qscore: #qscore=2 -np.log10(0.01)
						sc_classification.append(1)
						break
				except:
					print(l[0],p)
					sc_classification.append(0) #pass is right
					break
			else:
				sc_classification.append(0)
			sc_regression.append(t_mean)
		if len(sc_regression)<int(win/scale_win):
			for i in range(int(win/scale_win)-len(sc_regression)):
				sc_regression.append(0)
				sc_classification.append(0)
		if strand=="rev":
			sc_regression=sc_regression[::-1]
			sc_classification=sc_classification[::-1]
		if x in test_winlines:
			y_test_regression.append(sc_regression)
			y_test_classification.append(sc_classification)
		else:
			y_train_regression.append(sc_regression)
			y_train_classification.append(sc_classification)
	return y_train_regression,y_train_classification,y_test_regression,y_test_classification
def GetXset(winlines,win,dseq,strand,test_chr_list=[]):
	X_train=[]
	X_test=[]
	if test_chr_list:
		test_winlines=[]
		for x in winlines:
			x=x.rstrip()
			l=x.split("\t")
			if l[0] in test_chr_list:
				test_winlines.append(x)
	else:
		test_winlines=np.random.choice(winlines,int(len(winlines)*0.2),replace=False)
	if len(test_winlines)<1:
		print("test set error")
		sys.exit()
	for x in winlines:
		x=x.rstrip()
		l=x.split("\t")
		seq=dseq[l[0]][int(l[1]):int(l[2])]
		ss=GetOnehotSeq(seq,win,strand)
		if x in test_winlines:
			X_test.append(ss)
		else:
			X_train.append(ss)
	return X_train,X_test,test_winlines
def GetData(Genome,ChromeSize,Win,Scale,Strand,DripBdgFile,QpoisBdgFile,Qscore,Prefix,TestChrList=[]):
	Dseq={}
	for record in SeqIO.parse(Genome,"fasta"):
		Dseq[record.id]=str(record.seq)
	WinLines=os.popen("bedtools makewindows -g %s -w %s"%(ChromeSize,Win)).readlines()
	random.shuffle(WinLines)
	WinLines=[x.rstrip() for x in WinLines]
	X_train,X_test,TestWinlines=GetXset(WinLines,Win,Dseq,Strand,TestChrList)
	DripD=GetBdgd(DripBdgFile)
	QpoisD=GetBdgd(QpoisBdgFile)
	y_train_regression,y_train_classification,y_test_regression,y_test_classification=GetySet(WinLines,Win,Scale,DripD,QpoisD,TestWinlines,Strand,Qscore)
	del DripD,QpoisD
	gc.collect()
	np.savez_compressed("%s.npz"%Prefix,X_train=X_train,y_rlooptrain_classification=y_train_classification,y_rlooptrain_regression=y_train_regression,X_test=X_test,y_rlooptest_classification=y_test_classification,y_rlooptest_regression=y_test_regression)
def OnebaseEval(TruePeak,PredictPeak):
	TP=os.popen("bedtools intersect -a %s -b %s -nonamecheck|bedtools sort -i -|bedtools merge -i -|awk '{sum+=$3-$2}END{print sum}'"%(TruePeak,PredictPeak)).readlines()[0].rstrip()
	TP=int(TP)
	TotalPredict=os.popen("awk '{sum+=$3-$2}END{print sum}' %s"%(PredictPeak)).readlines()[0].rstrip()
	TotalTrue=os.popen("awk '{sum+=$3-$2}END{print sum}' %s"%(TruePeak)).readlines()[0].rstrip()
	FP=int(TotalPredict)-TP
	FN=int(TotalTrue)-TP
	Precision=TP/float(TP+FP)
	Recall=TP/float(TP+FN) #total_true
	F1=2*Precision*Recall/(Precision+Recall)
	Jaccard=os.popen("bedtools sort -i %s|bedtools jaccard -a %s -b - -nonamecheck"%(PredictPeak,TruePeak)).readlines()[1].rstrip().split()[2]
	return F1,Precision,Recall,Jaccard
def OverlapEval(TruePeak,PredictPeak):
	TP_predict=os.popen("bedtools intersect -wa -a %s -b %s -nonamecheck|sort -u|wc -l"%(PredictPeak,TruePeak)).readlines()[0].rstrip()
	TP_true=os.popen("bedtools intersect -wa -a %s -b %s -nonamecheck|sort -u|wc -l"%(TruePeak,PredictPeak)).readlines()[0].rstrip()
	TP_predict=int(TP_predict)
	TP_true=int(TP_true)
	TotalPredict=os.popen("wc -l %s"%PredictPeak).readlines()[0].rstrip().split()[0]
	TotalTrue=os.popen("wc -l %s"%TruePeak).readlines()[0].rstrip().split()[0]
	#FP=int(total_predict)-TP
	#FN=int(total_true)-TP wrong
	Precision=TP_predict/float(TotalPredict)
	Recall=TP_true/float(TotalTrue) #total_true
	F1=2*Precision*Recall/(Precision+Recall)
	return F1,Precision,Recall
def AbundanceEval(TrueBw,PredictBw,Prefix):
	os.system("multiBigwigSummary bins -p 30 -bs 128 -b %s %s -o %s_results.npz --outRawCounts %s_results.xls"%(TrueBw,PredictBw,Prefix,Prefix)) #128 need as a parameter
	os.system("plotCorrelation -in %s_results.npz --corMethod spearman --skipZeros --whatToPlot scatterplot -o %s_spearman.png --outFileCorMatrix %s_spearman.tab"%(Prefix,Prefix,Prefix))
	df=pandas.read_csv("%s_spearman.tab"%Prefix,sep="\t",header=0,index_col=0,comment="#")
	Spearman=df.iloc[0,1]
	return Spearman
def GetyScore(scale,chrom_size,predict_bed,have_scaled=True):
	y_scored={}
	if have_scaled:
		lines=open(predict_bed).readlines()
	else:
		lines=os.popen("bedtools makewindows -g %s -w %s|bedtools intersect -nonamecheck -a - -b %s -loj|bedtools groupby -i - -g 1,2,3 -c 8 -o mean"%(chrom_size,scale,predict_bed)).readlines() #-g -c need change?
	for x in lines:
		x=x.rstrip()
		l=x.split("\t")
		if have_scaled:
			y_scored[l[0]+"\t"+l[1]+"\t"+l[2]]=float(l[4])
		else:
			y_scored[l[0]+"\t"+l[1]+"\t"+l[2]]=float(l[3])
	return y_scored
def GetyTrue(scale,chrom_size,true_bed):
	y_trued={}
	lines=os.popen("bedtools makewindows -g %s -w %s|bedtools intersect -nonamecheck -a - -b %s -loj"%(chrom_size,scale,true_bed)).readlines()
	for x in lines:
		x=x.rstrip()
		l=x.split("\t")
		if l[-1]=="." or l[-1]=="-1":
			y_trued[l[0]+"\t"+l[1]+"\t"+l[2]]=0
		else:
			y_trued[l[0]+"\t"+l[1]+"\t"+l[2]]=1
	return y_trued
def GetPr(true_bed,predict_bed,chrom_size,scale,have_scaled=True):
	lines=os.popen("bedtools makewindows -g %s -w %s"%(chrom_size,scale)).readlines()
	y_trued=GetyTrue(scale,chrom_size,true_bed)
	y_scored=GetyScore(scale,chrom_size,predict_bed,have_scaled)
	y_true=[]
	y_score=[]
	for x in lines:
		x=x.rstrip()
		l=x.split("\t")
		y_true.append(y_trued[l[0]+"\t"+l[1]+"\t"+l[2]])
		y_score.append(y_scored[l[0]+"\t"+l[1]+"\t"+l[2]])
	#fpr,tpr,thresholds=roc_curve(y_true,y_score,pos_label=1)
	precision,recall,thresholds = precision_recall_curve(y_true,y_score)
	#auc=roc_auc_score(y_true,y_score)
	AP=average_precision_score(y_true,y_score)
	#np.save("thre.npy",thresholds)
	#print(auc)
	return precision,recall,thresholds,AP
def PlotPr(precisiond,recalld,APd,colord,sample_list,prefix):
	plt.rcParams['svg.fonttype'] = 'none'
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(1,1,1)
	ax.grid(True,linestyle='dotted')
	#color_list=["#C00000","#2E75B6","#548235","#BF9000"]
	#sample_list=["single_scale128","multi_scale128","single_base","multi_base"]
	for s in sample_list:
		ax.plot(recalld[s],precisiond[s],label='%s (AP=%s)'%(s,format(APd[s],'.4f')),linewidth=2,color=colord[s])
	#ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
	#ax.set_xlim([0.0, 1.0])
	ax.set_ylim(top=1.05)
	ax.set_xlabel('Recall')
	ax.set_ylabel('Precision')
	ax.set_title('Precision-Recall curve')
	ax.legend(loc="lower left")
	fig.savefig("%s_pr.png"%prefix,format='png')
	fig.savefig("%s_pr.svg"%prefix,format='svg')
	fig.clf()
	plt.close(fig)
def Evaluation(Prefix,HaveScale=True,PredAllBed=None,ChromeSize=None,TruePeak=None,PredPeak=None,TrueBw=None,PredBw=None):
	Fr=open(Prefix+"_eval.xls","w")
	if TruePeak and PredPeak:
		if PredAllBed:
			if not ChromeSize:
				print("ChromSize is not provided.")
				AP="-"
			else:
				Scale=128 #as parameter?
				a,b,c,AP=GetPr(TruePeak,PredAllBed,ChromeSize,Scale,HaveScale)
		else:
			AP="-"
		Fr.write("Method\tF1\tPrecision\tRecall\tJaccard\tAP\n")
		F1,Precision,Recall,Jaccard=OnebaseEval(TruePeak,PredPeak)
		Fr.write("onebase\t"+str(F1)+"\t"+str(Precision)+"\t"+str(Recall)+"\t"+str(Jaccard)+"\t-"+"\n")
		F1,Precision,Recall=OverlapEval(TruePeak,PredPeak)
		Fr.write("overlap\t"+str(F1)+"\t"+str(Precision)+"\t"+str(Recall)+"\t-\t"+str(AP)+"\n")
	if TrueBw and PredBw:
		Spearman=AbundanceEval(TrueBw,PredBw,Prefix)
		Fr.write("\n\nspearman:%s\n"%Spearman)
	Fr.close()
def PlotPrList(Scale,Prefix,Target): #HaveScale=True PredictAllBed=None
	Precisiond={}
	Recalld={}
	APd={}
	Colord={}
	Fr=open(Prefix+"_pr.xls","w")
	SampleList=[]
	for x in open(Target): #Name\tture_peak.bed\tpredict_all.bed\tHaveScale(0/1)\tChromeSize\tcolor
		x=x.rstrip()
		l=x.split("\t")
		Name=l[0]
		TruePeak=l[1]
		PredictBed=l[2]
		HaveScale=bool(int(l[3]))
		ChromeSize=l[4]
		SampleList.append(Name)
		Colord[Name]=l[-1]
		Precision,Recall,Thresholds,AP=GetPr(TruePeak,PredictBed,ChromeSize,Scale,HaveScale)
		Precisiond[Name]=Precision.tolist()
		Recalld[Name]=Recall.tolist()
		APd[Name]=AP
		for pi,ri,ti in zip(Precision,Recall,Thresholds):
			Fr.write(Name+"\t"+str(pi)+"\t"+str(ri)+"\t"+str(ti)+"\n")
	Fr.close()
	PlotPr(Precisiond,Recalld,APd,Colord,SampleList,Prefix)
	#elif PredictAllBed:
	#	precision,recall,thresholds,AP=GetPr(TruePeak,PredictAllBed,ChromeSize,Scale,HaveScale)
	#	precisiond[Prefix]=precision.tolist()
	#	recalld[Prefix]=recall.tolist()
	#	APd[Prefix]=AP
	#	sample_list.append(Prefix)
	#	colord[Prefix]="#C00000"
	#	for pi,ri,ti in zip(precision,recall,thresholds):
	#		Fr.write(Prefix+"\t"+str(pi)+"\t"+str(ri)+"\t"+str(ti)+"\n")
	#	Fr.close()
	#	PlotPr(precisiond,recalld,APd,colord,sample_list,Prefix)
def PlotDist(Prefix,TruePeak,PredPeak,TrueBw,PredBw,Extend=1000,Thread=12):
	os.system("bedtools intersect -nonamecheck -wa -a %s -b %s|sort -u|bedtools sort -i - >%s_overlap_predict.bed"%(PredPeak,TruePeak,Prefix))
	os.system("bedtools intersect -nonamecheck -wa -a %s -b %s|sort -u|bedtools sort -i - >%s_overlap_true.bed"%(TruePeak,PredPeak,Prefix))
	os.system("bedtools intersect -nonamecheck -v -a %s -b %s|sort -u|bedtools sort -i - >%s_nooverlap_predict.bed"%(PredPeak,TruePeak,Prefix))
	os.system("bedtools intersect -nonamecheck -v -a %s -b %s|sort -u|bedtools sort -i - >%s_nooverlap_true.bed"%(TruePeak,PredPeak,Prefix))
	os.system("computeMatrix scale-regions --missingDataAsZero -S %s -R %s_overlap_predict.bed %s_nooverlap_predict.bed %s_overlap_true.bed %s_nooverlap_true.bed -o %s_true_matrix.gz --startLabel start --endLabel end -b %s -a %s -m %s -p %s"%(TrueBw,Prefix,Prefix,Prefix,Prefix,Prefix,Extend,Extend,Extend,Thread))
	os.system("computeMatrix scale-regions --missingDataAsZero -S %s -R %s_overlap_predict.bed %s_nooverlap_predict.bed %s_overlap_true.bed %s_nooverlap_true.bed -o %s_predict_matrix.gz --startLabel start --endLabel end -b %s -a %s -m %s -p %s"%(PredBw,Prefix,Prefix,Prefix,Prefix,Prefix,Extend,Extend,Extend,Thread))
	os.system("plotHeatmap -m %s_true_matrix.gz -o %s_true_metaplot.svg --outFileNameMatrix %s_true_matrix.gz --startLabel start --endLabel end --perGroup -y 'R-loop level' --dpi 300 --samplesLabel 'R-loop level' --regionsLabel OverP NoverP OverT NoverT --colorMap RdBu_r --plotFileFormat svg --legendLocation none "%(Prefix,Prefix,Prefix))
	os.system("plotHeatmap -m %s_predict_matrix.gz -o %s_predict_metaplot.svg --outFileNameMatrix %s_predict_matrix.gz --startLabel start --endLabel end --perGroup -y 'R-loop level' --dpi 300 --samplesLabel 'R-loop level' --regionsLabel OverP NoverP OverT NoverT --colorMap RdBu_r --plotFileFormat svg --legendLocation none"%(Prefix,Prefix,Prefix))
def PlotGeneDist(Prefix,GeneBed,TrueFwdBw,TrueRevBw,PredFwdBw,PredRevBw,Extend=1000,Thread=12):
	sl=os.popen("awk -F'\t' '{print $6}' %s|sort -u"%(GeneBed)).readlines()
	d={}
	d["true"]={"fwd":TrueFwdBw,"rev":TrueRevBw}
	d["pred"]={"fwd":PredFwdBw,"rev":PredRevBw}
	if len(sl)==2:
		sl=[s.rstrip() for s in sl]
		if "+" in sl and "-" in sl:
			os.system("grep '+' %s >%s_positive.bed"%(GeneBed,Prefix))
			os.system("grep '\\-' %s >%s_negative.bed"%(GeneBed,Prefix))
			for p in d:
				for dr in ["fwd","rev"]:
					for zf in ["positive","negative"]:
						os.system("computeMatrix scale-regions -p %s --missingDataAsZero -S %s -R %s_%s.bed -bs 5 -b %s -a %s -m %s --skipZeros -o %s_%s_%s_%s.gz"%(Thread,d[p][dr],Prefix,zf,Extend,Extend,Extend,Prefix,p,zf,dr))
				os.system("computeMatrixOperations rbind -m %s_%s_positive_rev.gz %s_%s_negative_fwd.gz -o %s_%s_antisense.gz"%(Prefix,p,Prefix,p,Prefix,p))
				os.system("computeMatrixOperations relabel -m %s_%s_antisense.gz -o %s_%s_antisense_deal.gz --sampleLabels \"R-loop level\" --groupLabels \"%s_antisense\""%(Prefix,p,Prefix,p,p))
				os.system("computeMatrixOperations rbind -m %s_%s_negative_rev.gz %s_%s_positive_fwd.gz -o %s_%s_sense.gz"%(Prefix,p,Prefix,p,Prefix,p))
				os.system("computeMatrixOperations relabel -m %s_%s_sense.gz -o %s_%s_sense_deal.gz --sampleLabels \"R-loop level\" --groupLabels \"%s_sense\""%(Prefix,p,Prefix,p,p))
			os.system("computeMatrixOperations rbind -m %s_true_sense_deal.gz %s_true_antisense_deal.gz %s_pred_sense_deal.gz %s_pred_antisense_deal.gz -o %s_deal.gz"%(Prefix,Prefix,Prefix,Prefix,Prefix))
			os.system("plotHeatmap -m %s_deal.gz -o %s_sense_antisense_metaplot.svg --perGroup -y 'R-loop level' --dpi 300 --samplesLabel 'R-loop level' --regionsLabel TrueSense TrueAntisense PredictSense PredictAntisense --colorMap RdBu_r --legendLocation none --yMin 0 --startLabel TSS --endLabel TTS"%(Prefix,Prefix))
		else:
			print(sl)
			print("strand do not include + or -.")
	else:
		for p in d:
			for dr in ["fwd","rev"]:
				os.system("computeMatrix scale-regions -p %s --missingDataAsZero -S %s -R %s -bs 5 -b %s -a %s -m %s --skipZeros -o %s_%s_%s_matrix.gz"%(Thread,d[p][dr],GeneBed,Extend,Extend,Extend,Prefix,p,dr))
			os.system("computeMatrixOperations rbind -m %s_%s_fwd_matrix.gz %s_%s_rev_matrix.gz -o %s_%s_matrix.gz"%(Prefix,p,Prefix,p,Prefix,p))
			os.system("computeMatrixOperations relabel -m %s_%s_matrix.gz -o %s_%s_matrix_deal.gz --sampleLabels \"R-loop level\" --groupLabels \"%s\""%(Prefix,p,Prefix,p,p))
		os.system("computeMatrixOperations rbind -m %s_true_matrix_deal.gz %s_pred_matrix_deal.gz -o %s_deal.gz"%(Prefix,Prefix,Prefix))
		os.system("plotHeatmap -m %s_deal.gz -o %s_metaplot.svg --perGroup -y 'R-loop level' --dpi 300 --samplesLabel 'R-loop level' --regionsLabel True Predict --colorMap RdBu_r --legendLocation none --yMin 0 --startLabel TSS --endLabel TTS"%(Prefix,Prefix))
def GetSenseAntisenseCounts(BedFile,TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,Thread,Prefix):
	Sense=pandas.DataFrame()
	Antisense=pandas.DataFrame()
	DfBed=pandas.read_table(BedFile,sep="\t",header=None,index_col=[0,1,2],names=["chr","start","end","gene","score","strand"])
	for s in ["fwd","rev"]:
		if s=="fwd":
			TrueBwFile=TrueFwdBw
			PredBwFile=PredFwdBw
		elif s=="rev":
			TrueBwFile=TrueRevBw
			PredBwFile=PredRevBw
		cmd="multiBigwigSummary BED-file --BED %s -p %s --bwfiles %s %s --label observation prediction --outRawCounts %s_%s_counts.xls -o %s_%s.npz"%(BedFile,Thread,TrueBwFile,PredBwFile,Prefix,s,Prefix,s)
		os.system(cmd)
		os.system("sed \"1s/[#|']//g\" %s_%s_counts.xls >%s_%s_counts_deal.xls"%(Prefix,s,Prefix,s))
		Df=pandas.read_table("%s_%s_counts_deal.xls"%(Prefix,s),sep="\t",header=0,index_col=[0,1,2])
		DfUni=Df[~Df.index.duplicated(keep=False)]
		DfBedUni=DfBed[~DfBed.index.duplicated(keep=False)]
		#df_con=pandas.concat([df_uni,df_bed_uni],axis=1)
		print(Df.shape)
		print(DfUni.shape)
		print(DfBed.shape)
		print(DfBedUni.shape)
		DfR=pandas.concat([DfUni,DfBedUni],axis=1)
		print(DfR)
		labels=["gene"]+list(Df.columns)
		if s=="fwd":
			Sense=Sense.append(DfR[DfR["strand"]=="+"].loc[:,labels],ignore_index=True)
			Antisense=Antisense.append(DfR[DfR["strand"]=="-"].loc[:,labels],ignore_index=True)
		elif s=="rev":
			Sense=Sense.append(DfR[DfR["strand"]=="-"].loc[:,labels],ignore_index=True)
			Antisense=Antisense.append(DfR[DfR["strand"]=="+"].loc[:,labels],ignore_index=True)
	Sense.to_csv("%s_sense_counts.xls"%Prefix,sep="\t",index=False)
	Antisense.to_csv("%s_antisense_counts.xls"%Prefix,sep="\t",index=False)
	TrueX=np.log2(Sense["observation"]+1)
	PredictY=np.log2(Sense["prediction"]+1)
	PlotCorr(TrueX,PredictY,Prefix+"_sense")
	TrueX=np.log2(Antisense["observation"]+1)
	PredictY=np.log2(Antisense["prediction"]+1)
	PlotCorr(TrueX,PredictY,Prefix+"_antisense")
	#return Sense,Antisense
def PlotRloopCorr(Prefix,GeneBed,TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,ChromSize,GeneBody=None,TssExtend=None,TtsExtend=None,Thread=12):
	if GeneBody:
		GetSenseAntisenseCounts(GeneBed,TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,Thread,Prefix+"_genebody")
	if TssExtend:
		os.system("awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >%s_tss_%s.bed"%(GeneBed,ChromSize,TssExtend,Prefix,TssExtend))
		GetSenseAntisenseCounts("%s_tss_%s.bed"%(Prefix,TssExtend),TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,Thread,Prefix+"_tss_"+str(TssExtend))
	if TtsExtend:
		os.system("awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >%s_tts_%s.bed"%(GeneBed,ChromSize,TtsExtend,Prefix,TtsExtend))
		GetSenseAntisenseCounts("%s_tts_%s.bed"%(Prefix,TtsExtend),TrueFwdBw,PredFwdBw,TrueRevBw,PredRevBw,Thread,Prefix+"_tts_"+str(TtsExtend))
def PlotCorr(TrueX,PredictY,Name):
	xr=np.array([[i] for i in TrueX])
	yr=np.array([[i] for i in PredictY])
	regr = linear_model.LinearRegression()
	regr.fit(xr,yr)
	s=scipy.stats.spearmanr(TrueX,PredictY)
	plt.rcParams['svg.fonttype'] = 'none'
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(1,1,1)
	ax.grid(True,linestyle='dotted')
	bp=ax.scatter(TrueX,PredictY,alpha=0.3,marker=".",linewidths=0)
	py=np.array([i[0] for i in regr.predict(xr)])
	ax.plot(TrueX,py,color='red',linewidth=1)
	xi=max(TrueX)*1/5.
	yi=max(PredictY)*5/6.
	ax.text(xi,yi,"Spearman=%.4f"%s[0])
	ax.set_ylabel("log2(Prediction+1)")
	ax.set_xlabel("log2(Observation+1)")
	ax.set_title(Name)
	fig.savefig("%s_corr.svg"%Name,format="svg")
	fig.savefig("%s_corr.png"%Name,format="png")
	fig.clf()
	plt.close(fig)
