<p align="center">
  <img src="https://github.com/PEHGP/deepRloopPre/blob/main/logo.png?raw=true" alt="deepRloopPre"/>
</p>

## Introduction
R-loops are three-stranded nucleic acid structures composed of a DNA:RNA hybrid strand and a single-stranded DNA. R-loops are widespread in different species and participate in a variety of biological processes. The ssDRIP-seq technology developed in our laboratory can efficiently and strand-specifically detect the whole genome R-loops and was widely used in multiple species. Here, we develop a prediction tool deepRloopPre based on deep learning from the ssDRIP-seq data, which could predict the locations and profiles of strand-specific genome-wide R-loops. The deepRloopPre is only dependent on DNA sequences and its performance is better than other R-loop prediction tools.  
  
deepRloopPre contains the following files:
- `deepRloopPredict.py` (main executable script for predicting R-loop)
- `deepRloopData.py` (main executable script for getting your own training set)
- `deepRloopTrain.py` (main executable script for training your own model)
- `deepRloopEval.py` (main executable script for evaluating model)

In addition, we also provide a online [tool](http://bioinfor.kib.ac.cn/R-loopAtlas/deepRloopPre/), but only one sequence of 128kb can be predicted at a time.
## Installation
Create a new environment:
```
$ conda create -n deepRloopPre python=3.7
```
Install deepRloopPre:
```
$ conda activate deepRloopPre
$ git clone https://github.com/PEHGP/deepRloopPre
$ cd deepRloopPre
$ conda env update -n deepRloopPre -f ml.yaml
$ python setup.py install --user
```
## Quick start
### Just prediction
```
$ deepRloopPredict.py --h5 Arabidopsis_thaliana.hdf5 --fasta your_genome.fasta --prefix test
```
### Training your own model step by step
#### 1. Format your training data
Getting the input file `test_drip.bdg`
```
$ bamCoverage -v -p 30 -b test_fwd.bam -o test_drip.bdg --binSize 1 --effectiveGenomeSize 2131846805 --normalizeUsing RPGC --outFileFormat bedgraph
```
Getting the input file `test_qpois.bdg`
```
$ macs2 callpeak -f BAMPE --trackline -B -t test_fwd.bam -g 2131846805 -n test_fwd --keep-dup all
$ macs2 bdgcmp -t test_fwd_treat_pileup.bdg -c test_fwd_control_lambda.bdg -m qpois --o-prefix test
```
Getting the training data
```
$ deepRloopData.py --fasta genome.fasta --drip test_drip.bdg --qpois test_qpois.bdg --prefix test
```
#### 2. Training your model
```
$ deepRloopTrain.py --npz dataset.npz --prefix test
```
#### 3. Using your model
Getting `test_final_predict.bed` and `test_predict.bw`
```
$ deepRloopPredict.py --h5 your_model.hdf5 --fasta your_genome.fasta --prefix test
```
Getting `test_all_predict.bed`
```
$ deepRloopPredict.py --h5 your_model.hdf5 --fasta your_genome.fasta --threshold 0 --prefix test_all
```
#### 4. Evaluation your model
```
$ deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak test_final_predict.bed --truebw observed_drip.bw --predbw test_predict.bw --prefix test 
$ deepRloopEval.py plotpr --truepeak observed_peaks.bed --predpeak test_all_predict.bed --prefix test
```
## Prediction  
You can download the trained model [here](http://bioinfor.kib.ac.cn/R-loopAtlas/deepRloopPre/downloads.php).
```
$ deepRloopPredict.py -h
usage: $ deepRloopPredict.py --h5 model.hdf5 --fasta genome.fasta --prefix test

A tool for predicting R-loop profiles and R-loop locations.

Required arguments:
  --h5 H5               Saved model (default: None)
  --fasta FASTA         FASTA format file, usually a genome FASTA file. (default: None)
  --prefix PREFIX

Optional arguments:
  --help, -h            show this help message and exit
  --chromsize CHROMESIZE
                        If this file is not provided, all the sequences in FASTA file will be predicted. If this file is provided, only the sequences
                        contained in this file will be predicted. chrom_size.txt file has two columns separated by tabs. The first column is the sequence name
                        and the second column is the sequence length. (default: None)
  --strand {fwd,rev}    fwd:Predict current sequence. rev:Predict the reverse complement sequence. (default: fwd)
  --threshold THRESHOLD
                        The threshold score is used to judge whether an region is an R-loop region, the value range is from 0 to 1. (default: 0.5)
  --mem {high,low}      low option can reduce memory, but the program will run slower. (default: high)
  --version             show program's version number and exit
```
chrom_size.txt file has two columns separated by tabs. The file format is as follows:
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
```
Prefix is set to test and the `deepRloopPredict.py` output as follows:
- `test_predict.bw` R-loop profile prediction results of bigWig format
- `test_predict.bdg` R-loop profile prediction results of bedGraph format
- `test_predict.bed` original R-loop location prediction results of bed format
- `test_final_predict.bed` final R-loop location prediction results of bed format. Merging R-loop region with a distance of less than 300bp and then deleting region with a distance of less than 100bp.
## Training your own data
### 1. Format your training data
```
$ deepRloopData.py -h
usage: $ deepRloopData.py --fasta genome.fasta --drip test_drip.bdg --qpois test_qpois.bdg --prefix test

format the data for training.

Required arguments:
  --fasta FASTA         FASTA format file, usually a genome FASTA file. (default: None)
  --drip DRIP           R-loop profiles in bdg format obtained from ssDRIP-seq. (default: None)
  --qpois QPOIS         -np.log10(qvalue) in bdg format obtained from macs2. (default: None)
  --prefix PREFIX

Optional arguments:
  --help, -h            show this help message and exit
  --chromsize CHROMESIZE
                        If this file is not provided, all the sequences in FASTA file will be formated. If this file is provided, only the sequences
                        contained in this file will be formated. chrom_size file has two columns separated by tabs. The first column is the sequence name
                        and the second column is the sequence length. (default: None)
  --strand {fwd,rev}    fwd:format current sequence. rev:format the reverse complement sequence. (default: fwd)
  --logq QSCORE         The threshold score (-np.log10(qvalue)) is used to judge whether an region is an R-loop region. (default: 2)
  --chrlist CHRLIST [CHRLIST ...]
                        A list of space-delimited sequence names containing those sequence that should be used as validation set. If this parameter is not
                        provided, 20% of the data in the training set will be randomly selected as the validation set. (default: None)
  --version             show program's version number and exit
```
The input `test_drip.bdg` can be obtained from BAM file:
```
$ bamCoverage -v -p 30 -b test_fwd.bam -o test_drip.bdg --binSize 1 --effectiveGenomeSize 2131846805 --normalizeUsing RPGC --outFileFormat bedgraph
```
`--effectiveGenomeSize` parameter needs to be set to your own.  
`test_fwd.bam` BAM needs to be obtained by aligning raw data to the genome. If you use ssDRIP-seq, you can refer to this [wiki](https://github.com/PEHGP/ssDripPipeline/wiki).  
Of course, you can also generate `test_drip.bdg` based on the [bedGraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) format by yourself.  
The input `test_qpois.bdg` can be obtained as follows:
```
$ macs2 callpeak -f BAMPE --trackline -B -t test_fwd.bam -g 2131846805 -n test_fwd --keep-dup all
$ macs2 bdgcmp -t test_fwd_treat_pileup.bdg -c test_fwd_control_lambda.bdg -m qpois --o-prefix test
```
`-g` means genome size and needs to be set to your own.  
`-f` means format of tag file and needs to be set to your own. If it is BAM for paired-end reads, **BAMPE** is used.  
Prefix is set to test and the `deepRloopData.py` output as follows:
- `test.npz`
### 2. Training your model
```
$ deepRloopTrain.py -h
usage: $ deepRloopTrain.py --npz dataset.npz --prefix test

Building your own deepRloopPre model.

Required arguments:
  --npz NPZ        Training set in npz format. (default: None)
  --prefix PREFIX

Optional arguments:
  --help, -h       show this help message and exit
  --epoch EPOCH    Epoch (default: 1000)
  --bs BATCHSIZE   batch size (default: 20)
  --version        show program's version number and exit
```
The input `dataset.npz` can be obtained from `deepRloopData.py`  
The output is a series of files in hdf5 format, among which the newly generated hdf5 file is the final trained model.  
### 3. Evaluation your model
```
$ deepRloopEval.py -h
usage: deepRloopEval.py [-h] [--version]  ...

Evaluation of prediction results.This program can output F1-score, Precision, Recall for R-loop location and Spearman for R-loop abundance or draw a precision-recall curve.
A detailed sub-commands help is available by typing:
        deepRloopEval.py getvalue -h
        deepRloopEval.py plotpr -h

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

commands:
  
    getvalue  getting the F1-score, Precision, Recall, Spearman.
    plotpr    plot precision-recall curve.

example usage:
        deepRloopEval.py getvalue --truepeak observed_peaks.bed --predpeak predicted_peaks.bed --prefix test
        deepRloopEval.py plotpr --truepeak observed_peaks.bed --predpeak predicted_all.bed --prefix test
```
The `observed_peaks.bed` is the R-loop location file obtained by ssDRIP-seq, or the true R-loop location file prepared by yourself. The `predicted_peaks.bed` is the predicted R-loop locations.  
`getvalue` subcommand can calculate the spearman of the R-loop profiles, and you need to use `--truebw` and `--predbw`.
## Third-party software
- [BEDTools](https://bedtools.readthedocs.io/en/latest/index.html)
- [deepTools](https://deeptools.readthedocs.io/en/develop/)
- [bedGraphToBigWig](http://rohsdb.cmb.usc.edu/goldenPath/help/bigWig.html)
- [MACS2](http://github.com/taoliu/MACS/)
## Citing this work
