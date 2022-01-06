# deepRloopPre
The R-loop is a three-stranded nucleic acid structure composed of a DNA:RNA hybrid strand and single-stranded DNA. R-loops are widespread in different species and participate in a variety of biological processes. The ssDRIP-seq technology developed in our laboratory can efficiently and strand-specifically detect the whole genome R-loop, and is widely used in multiple species. We used ssDRIP-seq data based on deep learning to develop the prediction tool deepRloopPre, which can predict the locations and profiles of strand-specific R-loops in the whole genome using only sequences, and is better than other R-loop prediction tools.  
deepRloopPre contains the following files:
- `deepRloopPredict.py` (main executable script for predicting R-loop)
- `deepRloopData.py` (main executable script for getting your own training set)
- `deepRloopTrain.py` (main executable script for training your own model)
- `deepRloopEval.py` (main executable script for evaluating models)
## Installation
Create a new environment:
```
$ conda create -n deepRloopPre python=3.7
```
Import configuration:
```
$ conda env export -n deepRloopPre --file ml.yaml
```
Install deepRloopPre:
```
$ git clone git@github.com:PEHGP/deepRloopPre.git
$ python setup.py install --user
```
## Commands
**prediction of R-loop profiles and R-loop locations**
The trained model is in the model folder.There are three training sets of models that can be selected for prediction.
```
deepRloopPredict.py -h
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
prefix is test and the `deepRloopPredict.py` output as follows:
`test_predict.bw` R-loop profile prediction results of bigWig format
`test_predict.bdg` R-loop profile prediction results of bedGraph format
`test_predict.bed` original R-loop location prediction results of bed format
`test_final_predict.bed` final R-loop location prediction results of bed format. Merging R-loop region with a distance of less than 300bp and then deleting region with a distance of less than 100bp.
**Training your own data**
1. Format training data
```
```
## Third-party software
- bedtools
- deeptools
- bedGraphToBigWig