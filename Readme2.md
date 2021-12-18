#Python scripts for comparing A2010 and PAM models
Martin Montes, SSAI-NASA
Lanham, April 30, 2020

Extraction of AERONET parameters
The main script is located in parseAEROglobe.py and can extract all AERONET aerosol properties excluding phase functions. Main calls the method parse_aero and is calling an iterative function (parse_line) to parse line-by-line each AERONET file with *.all extension. The result is saved as a csv file. Notice that an additional aerosol propertie is computed for each record based on AODa measured at two wavelengths (the Angstrom exponent for absorption within 870-1020 nm).

Creation of Mie input files 
These routines are organized in make_MieInput.py and include code to produce formatted Mie input files for MODTRAN and FORTRAN. make_MieInput.py calls inp_miefile() to initiate the processing of files and includes two main methods (rv_to_rm2 and nd_frac2) to convert median radius from volume to number density space and calculate number density fractions, respectively. There are two directories to input cluster’s ids, aerosol properties and one directory to output Mie input files for MODTRAN and FORTRAN simulations. Standard AERONET wavelengths are interpolated to 17 wavelengths extracted from A2010 models. Volume fractions levels and integration intervals for PSD are obtained from A2010 as well. Lastly, Mie output files are created for each PAM cluster.

Extraction of Mie output files
This task is subdivided between two scripts devoted to parse Mie results from PAM (extract_MIEPAM.py) and A2010 (extract_MIEZIA.py) models. In all cases, the extracted information corresponds to Mie simulations of radiative aerosol properties as computed from Jean-Claude Roger code in FORTRAN. For PAM, there is an input file as txt for each cluster. For A2010, there is only one big output file with Mie simulations results. The parsed results are saved as csv files.

Spectral interpolation of A2010 to customized wavelengths
The script interpZIA.py adjust the aerosol properties with spectral dependency of A2010 models to those corresponding to 6s (default configuration) or AERONET (0.44, 0.675,0.870 and 1.02 microns). The function performing the cubic interpolation is interp_ref. Beyond 1.02 microns, the interpolation results are not reliable and should be ignored.

Similarity metrics
Aerosol properties of PAM and A2010 models are compared for each number density fraction (NDF) using the script simPAMZIA.py. The range of comparisons must be constrained a posteriori based on mean NDF values corresponding to each PAM cluster. simPAMZIA.py has one method (cosim()) for calculating the cosine similarity angle also known as spectral angle. The remaining metrics for comparing spectral and absolute values of aerosol attributes are computed based on pointers to python libraries. There are two input (unsorted Mie files for PAM and sorted for A2010) and two output (sorted Mie files for PAM and metrics results) directories. After reading the PAM and A2010 models, the following step is to select which aerosol property will be compared (ssa,g, Kext or Ka?). Notice that these comparisons are spectral containing 17 wavelengths per aerosol property (e.g., the spectral curve of SSA from PAM and A2010 is compared). Five metrics are computed for each inter-model comparison: cosine similarity index, root mean square error, median absolute error, normalized root mean square error and normalized median absolute error (COSIM, RMSE, MAE, NRMSE and NMAE, respectively). Spectral differences are quantified with COSIM. Conversely, RMSE, MAE, NRMSE and NMAE metrics are used to quantify absolute differences.

Clusters statistics
These methods are described in clus_stat.py and allow the calculation of number of clusters, number of records, dominant cluster, number of records per dominant cluster and dominant cluster fraction per month. There is one input path for two files needed to obtain AERONET sites and clusters id along site names and aerosol properties.
