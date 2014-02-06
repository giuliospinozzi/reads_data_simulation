###Header################################################
header = """

+------------------------------------------------------+
 Module: simulation_module
 Author: Stefano Brasca
 Date:  January 31th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions used in reads data
  	simulations
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import sys
import random
import scipy
from scipy import stats
###############################

###Import Module(s)#####################
import classes_for_data_simulation
########################################






def input_controls (source_distribution, distribution_parameters, span, n_of_events):
	'''
	*** Type control on arguments given in input ***
	          [for GENERATE_IS_HISTOGRAM]
	'''
	# distribution_parameters must be a dictionary
	if (type(distribution_parameters) is not dict):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'distribution_parameters' argument must be a dictionary."
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
	# distribution_parameters[source_distribution] must be a dictionary too
	if (type(distribution_parameters[source_distribution]) is not dict):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'distribution_parameters[{0}]' must be a dictionary.".format(source_distribution)
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

	# span must be an integer
	if (type(span) is not int):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'span' argument must be an integer, not '{0}'.".format(str(type(span)))
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
	# span must be >= 1
	if (span < 1):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'span' argument must be an integer >= 1, not '{0}'.".format(str(span))
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

	# n_of_events must be an integer
	if (type(n_of_events) is not int):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'n_of_events' argument must be an integer, not '{0}'.".format(str(type(n_of_events)))
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
	# n_of_events must be >= 1
	if (n_of_events < 1):
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'n_of_events' argument must be an integer >= 1, not '{0}'.".format(str(n_of_events))
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

	return 0

def distribution_parameters_controls(distribution_parameters, source_distribution, parameter):
	'''
	*** Perform controls on distribution_parameters dictionary given in input ***
	          [for GENERATE_IS_HISTOGRAM]
	'''
	try:
		distribution_parameters[source_distribution][parameter]
	except:
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'distribution_parameters[{0}][{1}]' doesn't exist.".format(source_distribution, parameter)
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")







def GENERATE_IS_HISTOGRAM (source_distribution, distribution_parameters, span, n_of_events):
	
	'''
	*** Brief description ***

	INPUT: - source_distribution: the string defining the distribution called to generate 
	                              the histogram; in case of incorrect type or unavailable
	                              choice sys.exit('some-explanation') will be called

	       - distribution_parameters: nested dictionary of the parameters defining
	                                  source_distribution;
	                                  external key = source_distribution
	                                  internal key = 'name-of-desired-parameter'
	                                  
	                                  e.g. gauss:
	                                  st_dev = distribution_parameters['gauss']['st_dev']

	       - span: integer representing the maximum number of bin allowed for the histogram;
	               span must be >= 1

	       - n_of_events: integer representing the number of discretized_realizations used to generate
	       				  the histogram (sum of occurrences); n_of_events must be >= 1

	OUTPUT: - Histogram object (see classes_for_data_simulation module).
			  Such an object has following attributes:

			  	- source_distribution, distribution_parameters, n_of_events like data
			  	  given in input

				- bins: list of integers, from 0 to span-1; likewise len(bins) == span.
						It represents the list of 'bin centers' for the just generated histogram
						
						[Note for symmetric 'source_distribution']:
						bins goes from 0 to span-1 (if span was odd) or to span (if 
			            span was even); likewise len(bins) == span (odd) or span+1 (even).
						As an even span requires corrections, a warning message will be 
						produced

				- discrete_realizations: list of realization from source_distribution, discretized 
										 and cleaned (realizations beyond bin edges are discarded)
										 [Note: from cleaned_discretized_realizations variable of
										  this function]

				- occurrencies: list of oredered occurrencies for each item in bins

				- frequencies: list of oredered frequencies for each item in bins
						
				- expected_value: integer, representing the expected value underlaying the
				                  histogram computation

				- discrete_realizations_beyond_edges: dictionary: 'string, +/- distance from expected_value'
											as key, 'integer, number of occurrencies found' as item

				- p-value: float, result of GOF test

				- another_GOF_indicator: None for general purpose; KS_test if 'gauss' is the
										 source_distribution

			[further data are available if you need them:
			 - span: like the one given in input, sometimes modified to fit chosen 'source_distribution'
			 - occurrencies: list of oredered occurrecies for each item in bins]

	NOTE: 1) about (source_distribution == 'gauss') case:
			 A Kolmogorov-Smirnov test for goodness of fit is performed
			 [see e.g. - http://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
			 		   - http://ocw.mit.edu/courses/mathematics/18-443-statistics-for-applications-fall-2006/lecture-notes/lecture14.pdf]
			 Remember that the null-hypothesis is that the two sample tested are coming
			 from the same distribution; we can reject this hypotesis if p-value is lower than a threshold (e.g. <0.05).
			 "The higher the p-value, the more 'gaussian' is the shape"

	'''

	### GAUSS DISTRIBUTION CASE #######################################################################################################################################################################
	if (source_distribution == 'gauss'):

		# Input arguments controls
		input_controls (source_distribution, distribution_parameters, span, n_of_events)
		distribution_parameters_controls(distribution_parameters, source_distribution, 'st_dev')

		# span (number of bins) should be odd or will be corrected adding '1'; [WARNING] will be produced but  the function goes on
		# Symmetrize span, if needed
		if (span % 2 == 0):
			span = span + 1
			message = "[WARNING - GENERATE_IS_HISTOGRAM call]\t'span' value has been symmetrized to fit 'gauss' shape: from '{0}'' to '{1}'".format(str(span-1), str(span))
			print "\n\n" + message + "\n\t[CONTINUE]\n\n"

		# Get bins
		bins = range(0, span) ### VARIABLE TO RETURN

		# Get expected_value
		expected_value = (span-1) / 2 ### VARIABLE TO RETURN
		if (expected_value != float((span-1))/2):
			error_message = "[ERROR]\tTroubles in GENERATE_IS_HISTOGRAM: can't retrieve the correct 'expected_value' for 'gauss' distribution"
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
		# Get expected_value_index
		expected_value_index = bins.index(expected_value)

		# Get st_dev
		st_dev = distribution_parameters[source_distribution]['st_dev']
		##################################################################################################################################
		### Perhaps further controls might be useful about st_dev type, value and value-effect over outcomes (with respect to span/bins) #
		##################################################################################################################################

		# Generate discretized_realizations
		realizations = []
		discretized_realizations = []
		for i in range(0, n_of_events):
			r = random.gauss(expected_value,st_dev)
			realizations.append(r)
			discretized_realizations.append(int(r))

		###################################################################################################################
		### From now on, discretized_realizations falling beyond bin boundaries (bins[0] and bins[-1]) are not accounted! #
		###################################################################################################################

		# Compute occurrencies
		occurrencies = [] ### VARIABLE TO RETURN
		total_occurrencies = 0
		for bin in bins:
			bin_occurrence = discretized_realizations.count(bin)
			occurrencies.append(bin_occurrence)
			total_occurrencies = total_occurrencies + bin_occurrence

		# Accounting for discretized_realizations falling beyond bin boundaries, then clean discretized_realizations
		discrete_realizations_beyond_edges = {}   # key = string, +/- distance from expected_value
												 # item = integer, number of occurrencies found
		cleaned_discretized_realizations = [] ### VARIABLE TO RETURN
											  # The same as discretized_realizations but without realizations falling beyond bin boundaries
		for discrete_realization in discretized_realizations:
			if (discrete_realization not in bins):
				# Update discrete_realizations_beyond_edges
				key = str(int(discrete_realization - expected_value)) # int should be redundant
				if (discrete_realizations_beyond_edges.has_key(key)):
					item = discrete_realizations_beyond_edges[key] + 1
					discrete_realizations_beyond_edges[key] = item
				else:
					item = 1
					discrete_realizations_beyond_edges.update({key:item})
			else:
				# Update cleaned_discretized_realizations
				cleaned_discretized_realizations.append(discrete_realization)

		# Compute frequencies
		frequencies = [] ### VARIABLE TO RETURN
		for occurrence in occurrencies:
			frequencies.append(float(occurrence)/total_occurrencies)

		# 'Fit' data and test goodness-of-fit
		realizations_standardized = [ float(r - expected_value)/float(st_dev) for r in realizations]
		KS_test, p_value = stats.kstest(realizations_standardized, 'norm') ###VARIABLES TO RETURN
		###############################################################################
		### Kolmogorov-Smirnov test for goodness of fit is done over 'realizations'   #
		### (not discrete)                                                            #
		###############################################################################

		# CREATE FINAL OBJECT TO RETURN #
		return classes_for_data_simulation.Histogram(source_distribution, distribution_parameters, n_of_events, expected_value, bins, cleaned_discretized_realizations, frequencies, occurrencies, discrete_realizations_beyond_edges, p_value, KS_test)

		###############################################################################################################################################################################################




	### WHATEVER DISTRIBUTION CASE ####################################################################################################################################################################
	elif (source_distribution == 'whatever'):

		# Input arguments controls
		input_controls (distribution_parameters, span, n_of_events)
		distribution_parameters_controls(distribution_parameters, source_distribution, 'somewhat-parameter-name')
		
		return 0
	###################################################################################################################################################################################################



	### SOME TROUBLES OCCURRED IN source_distribution CHOICE ###
	else:

		# source_distribution is not a string -> quit
		try:
			str(source_distribution)
		except:
			error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: 'source_distribution' argument must be a string."
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

		# source_distribution is not a valid string -> quit
		error_message = "[ERROR]\tGENERATE_IS_HISTOGRAM bad calling: '{0}' does not correspond to any distribution available.".format(source_distribution)
		sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")