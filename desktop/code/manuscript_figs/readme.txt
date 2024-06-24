G Oldford, March 2024
g.oldford@dfo-mpo.gc.ca; greig.oldford@gmail.com

Python notebooks (Python 3.9x) have .ipynb extension
Any files with .ksh extension are meant to be run on unix server
The .py files were developed using PyCharm
.qgz file extension is QGIS map file

please contact greig.oldford@dfo-mpo.gc.ca if you would like source data

Figs 3 and 4 generated using scripts on unix server hosting model outputs
Requires custom python library currently referred to as 'pyap', configured for the server
Please contact greig.oldford@dfo-mpo.gc.ca or michael.dunphy@dfo-mpo.gc.ca for assistance

Figure 4 - evaluation using CTDs for each model for 2016 - 2018 experimental eval period
	 - must be run on server
	 - to run (unix): 
	 -      go to: project/6006412/goldford
	 - 	module load python/3.8.10
	 -      module load StdEnv/2020
	 - 	module load scipy-stack/2020b
	 - 	edit /project/6006412/goldford/gitlab_pyap_fork20230331/analysispkg/pkg_plot_cast.py
	 - 	python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/compare_several_GO.yaml -i CTD -j 3
	 -  	output to: project/6006412/goldford/ANALYSIS/SalishSea1500-ALLRUNS_test/CTD/hindcast/ProfilePDF_comparisons



