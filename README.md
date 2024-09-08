Files to support manuscript submitted to Geoscientific Model Development related to the HOTSSea v1 model. GMD-2024-58.  
  
 

repo structure
<code>
desktop  
[desktop code and resources for analysis and plots]  
|  
|--> code  
|   |  
|   |--> manuscript_figs  
|   |	[desktop PC code for most manuscript plots; ORAS5.. code and Fig 3 code (see pypkg/etc) are meant to run on server]  
|   |  
|   |--> data_prep  
|   |	[notebooks etc for prepping data]  
|  
|--> data  
|   |  
|   |--> bathymetry  
|   |  
|	|--> eval  
|	|  
|	|--> grid  
|	|  
|	|--> mesh mask  
|	|  
|	|--> outputs  
|	|	[due to space limits, contains only mod-obs extraction from running pypkg]  
|	|  
|	|--> reference  
|  
|--> figs  
|	[output folder for manuscript plots]  
|  
|--> venv  
|	[python libraries and virtual environment files]  
|  
serverside  
[server code (computecanada, graham, unix) and resurces for running model, mod-obs extraction, analysis, and plots]  
|  
|-->NEMO  
|	|[NEMO v3.6 code, config, namelists, and some smaller forcings, from mdunphy/ folder on server]  
|	|  
|	|-->control  
|	|	[control files for each model experiment and mapping of experiment coding]  
|	|  
|	|-->Forcing  
|	|	[bdy, surface, runoff; most had to be stored in separate Zenodo repos due to size constraints]  
|	|  
|	|-->Software  
|		[NEMO v3.6 software]  
|  
|-->other scripts  
|	|[custom scripts (other than pypkg) used for data prep, extraction, and analysis]  
|   |
|   |--> 'ORAS5..' - temperature bias correction (Fig 3)
|
|-->pypkg  
|	|[NEMO-model python 'analysis package' used for mod-obs and statistics, under development, not for redist]
|	|
|	|-->etc
|	   [batch job submission scripts. Look to plot_CTD_comparison.sh for script to generate Fig 3]
|  
|-->ANALYSIS  
	|[outputs from running pypkg scripts, mainly model-obs stats] 
	|
	|--> SalishSea1500-ALLRUNS_test/CTD/hindcast
	|     [outputs of running the analysis used to generate Fig 3]
	|
	|-->
</code>
