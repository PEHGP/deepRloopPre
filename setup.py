import re
from setuptools import setup, find_packages
from deepRloopPre._version import __version__
setup(
	name='deepRloopPre',
	version=__version__,
	author='LiKuan',
	author_email="396777306@qq.com",
	packages=find_packages(),
	scripts=['bin/deepRloopData.py','bin/deepRloopTrain.py','bin/deepRloopPredict.py','bin/deepRloopEval.py'],
	#include_package_data=True,
	url='https://github.com/PEHGP/deepRloopPre',
	license='GPL-3.0',
	description='A tool for predicting R-loop profiles and R-loop locations.',
	long_description='',
	classifiers=[
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
	install_requires=[
		"pandas >= 1.1.3",
		"biopython >=1.78",
		"matplotlib==3.2.2",
		"tensorflow==2.4.1",
		"scikit-learn>=0.23.2",

	],
	zip_safe=False,
)
