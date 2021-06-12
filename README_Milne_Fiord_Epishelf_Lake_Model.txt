README 
Milne Fiord Epishelf Lake Model

Last updated 2021-05-13

Please cite this dataset as follow:
Jérémie Bonneau, Bernard E. Laval, Derek Mueller, Andrew K. Hamilton, Drew M. Friedrichs and  Alexander L. Forrest,. (2021). Winter Dynamics in an Epishelf Lake: Quantitative Mixing Estimates and Ice Shelf Basal Channel Considerations. Journal of Geophysical Research: Oceans. Submitted 2021-03-11. Manuscript #2021JC017324. In review. 


-----
DESCRIPTION:
This is the model used by Bonneau et al. (see above). This archive version reproduces the results. This model uses the mooring data in Milne Fiord from 2011 to 2019 to estimate the amount of mixing happening in the epishelf lake and below, and the minimum draft of Milne Ice Shelf. 
This model was created using matlab.
This model uses the gsw toolbox; it will not run without it --> http://www.teos-10.org/software.htm#1 
For details on the model please refer to the above publication. 
For a workflow schematic of the model, please refer to the model_workflow.png 

SCRIPTS
	model_11_12  : script to run the model for the 2011-2012 winter
	model_12_13  : script to run the model for the 2012-2013 winter
	model_13_14  : script to run the model for the 2013-2014 winter
	model_14_15  : script to run the model for the 2014-2015 winter
	model_15_16  : script to run the model for the 2015-2016 winter
	model_16_17  : script to run the model for the 2016-2017 winter
	model_17_18  : script to run the model for the 2017-2018 winter
	model_18_19  : script to run the model for the 2018-2019 winter
	
FUNCTIONS
	Pre-process
		meshing          : creates a mesh for the model
		IC_mooring_ctd   : writes the salinity and temperature initial conditions on the mesh
		TBC              : writes the temperature boundary conditions

	Process
		find_daily_K     : computes the salinty and temperature for the specified time interval looping through the possible values of mixing coefficients. Uses K_N2_custom_fast, CN_Nstep_fast, Outflow_weir and model_score
		model_score      : computes the coefficient of agreement (C_a) between the model and the data from the mooring. Uses CT_Eval_fast and SA_Eval_fast
		CN_Nstep_fast    : solves N time steps dt of the 1D transport eddy diffusivity equation, no advection term, using a Crank-Nicolson scheme.
		K_N2_custom_fast : computes the combined (molecular+turbulent) mixing for heat and salt using the N2 frequency profile. Uses K_SA_from_K_CT
		N2_profile       : computes the N2 frequency of the water column using the gsw toolbox
		Outflow_weir     : simulates the outflow of water from the epishelf lake using an inverse rectangular equation
		K_SA_from_K_CT   : computes the diffusivity (molecular+turb) for salinity from the turbulent eddy diffusivity for heat using the parameterization proposed by Jackson and Rehman (2014); Experiments on differential scalar mixing in turbulence in a sheared, stratified flow. Journal of Physical Oceanography, 44(10), 2661-2680.

	Post-process
		CT_Eval_fast     : computes the temperature difference between the model and the mooring thermistors data
		SA_Eval_fast     : computes the salinity difference between the model and the mooring data

DATA
	K_SA_from_K_CT_struct2interp: data to interpolate K_SA from K_CT
	SA_day                      : salinity data from the mooring, daily averaged, to evaluate the model
	SA_IC_ctd                   : salinity data to write the initial conditons
	Tmatrix_0_30_30m_10cm       : temperature data from the mooring to write boundary conditions; too big for Github, download at https://drive.google.com/file/d/1uWGEAiKqjb5h5NCoJFM_BDw-OVdUIgCm/view?usp=sharing 
	Tmatrix_0_30_day_10cm       : temperature data from the mooring, daily averaged, to write initial conditions and evaluate the model

RESULTS
	01sept2011_01jun2012_h067_Cb25 : structure with the output of the model for winter 2011-2012
	02sept2012_10jun2013_h075_Cb49 : structure with the output of the model for winter 2012-2013
	11aug2013_12jun2014_h075_Cb74  : structure with the output of the model for winter 2013-2014
	10sept2014_10jun2015_h070_Cb76 : structure with the output of the model for winter 2014-2015
	06sept2015_04jun2016_h053_Cb34 : structure with the output of the model for winter 2015-2016
	25sept2016_01jun2017_h079_Cb55 : structure with the output of the model for winter 2016-2017
	25aug2017_01jun2018_h073_Cb72  : structure with the output of the model for winter 2017-2018
	05sept2018_01jun2019_h078_Cb55 : structure with the output of the model for winter 2018-2019

----
ACKNOWLEDGEMENTS

Funding from:
-The Natural Sciences and Engineering Research Council of Canada
-Canada Foundation for Innovation
-ArcticNet Network of Centres of Excellence of Canada 
-The University of British Columbia
-Carleton University

Logistical support was provided by:
-Polar Continental Shelf Program

----
ARCHIVE:
The Milne Fiord CTD data is permanently archived and updated as necessary at the Polar Data Catalogue (CCIN 12101 and 12102; https://polardata.ca/).

----
CONTACT:

Please direct inquiries to:
Jérémie Bonneau, The University of British Columbia, Vancouver, Canada (jeremie.bonneau@alumni.ubc.ca)
