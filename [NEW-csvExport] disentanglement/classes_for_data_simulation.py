###Header################################################
header = """

+------------------------------------------------------+
 Module: classes_for_data_simulation
 Author: Stefano Brasca
 Date:  January 31th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains classes used in reads data
  	simulations (simulation_module.py)
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import os
import copy
import random
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
###############################

###Import Module(s)#####################
#import yyyy
########################################






class IS_Histogram:

	def __init__ (self, source_distribution, distribution_parameters, n_of_events, expected_value, bins, discrete_realizations, frequencies, occurrencies, discrete_realizations_beyond_edges, p_value, another_GOF_indicator = None):

		self.source_distribution = source_distribution # String
		self.distribution_parameters = distribution_parameters[source_distribution] # Dictionary
		self.n_of_events = n_of_events # Integer
		self.expected_value = expected_value # Integer
		self.bins = bins # List of integers, starting from 0
		self.discrete_realizations = discrete_realizations # list of realization from source_distribution, discretized 
														   # and then cleaned (realizations beyond bin edges are discarded)
		self.occurrencies = occurrencies # List of integers, paired with bins
		self.frequencies = frequencies # List of float, paired with bins
		self.discrete_realizations_beyond_edges = discrete_realizations_beyond_edges # Dictionary 	# key = string, +/- distance from expected_value
																									# item = integer, number of occurrencies found
		self.p_value = p_value
		self.another_GOF_indicator = another_GOF_indicator #KS_test if source_distribution='gauss'
		
		# Post Amplification
		self.amplified = False
		self.amplification_bias = None
		self.slippage_bias = None
		self.n_of_mapped_sequencies = None # Integer; 'Amplified analogous' of n_of_events
		self.sequence_count = [None]*len(occurrencies) # List of integers, 'amplified analogous' of occurrencies, paired with bins
		self.abundances = [None]*len(frequencies) # List of float, 'amplified analogous' of frequencies, paired with bins
		self.p_value_post_amplification = None
		self.another_GOF_indicator_post_amplification = None #KS_test if source_distribution='gauss'




	def amplify (self, minimum_amplification_factor=1, maximum_amplification_factor=100, amplification_bias=True, slippage_bias=True, minimum_splippage_percentage=0, maximum_splippage_percentage=1):

		# This method works only with unamplified IS_Histograms
		if (self.amplified == False):

			### Control and Cast
			if ((type(minimum_amplification_factor) is not int) or (type(maximum_amplification_factor) is not int)):
				minimum_amplification_factor = int(minimum_amplification_factor)
				maximum_amplification_factor = int(maximum_amplification_factor)
				print "\n[WARNING] 'Amplify()' method of IS_Histogram Class supports only integer arguments: cast has been done."

			### Create a copy of IS_Histogram object
			amplified_IS_Histogram = copy.deepcopy(self) ### OBJECT TO RETURN

			# Tick amplification
			amplified_IS_Histogram.amplified = True

			# Inizialize results collectors
			amplified_discrete_realizations = copy.deepcopy(self.discrete_realizations)
			total_amplified_n_of_events = copy.deepcopy(self.n_of_events) ### VARIABLE TO RETURN
			amplified_occurrencies = copy.deepcopy(self.occurrencies) ### VARIABLE TO RETURN

			### Introduce PCR-AMPLIFICATION BIAS for each realization
			if (amplification_bias == True): 

				# Tick amplification bias
				amplified_IS_Histogram.amplification_bias = True

				# Set total_amplified_n_of_events
				total_amplified_n_of_events = 0
				amplified_discrete_realizations = []
				amplified_occurrencies = [0]*len(self.occurrencies)

				# Get a copy of self.discrete_realizations (to work on it)
				discrete_realizations = copy.deepcopy(self.discrete_realizations)

				# Loop: until discrete_realizations local variable is []  # Comupte: - amplified_discrete_realizations (for GOF test)
																		  #          - total_amplified_n_of_events -> amplified_IS_Histogram.n_of_mapped_sequencies
																		  #          - amplified_occurrencies -> amplified_IS_Histogram.sequence_count (then amplified_frequencies -> amplified_IS_Histogram.abundances)
				while (discrete_realizations != []):

					# Randomly choose a realization
					selected_realization = random.choice(discrete_realizations)
					# remove it (like pop)
					discrete_realizations.remove(selected_realization)
					# Amplify: relate a number to selected_realization
					current_amplification_factor = np.random.random_integers(minimum_amplification_factor,maximum_amplification_factor)
					# Update amplified_discrete_realizations
					amplified_discrete_realizations = amplified_discrete_realizations + [selected_realization]*current_amplification_factor
					# Update total_amplified_n_of_events
					total_amplified_n_of_events = total_amplified_n_of_events + current_amplification_factor
					# Sum the result in the related amplified_occurrency
					amplified_occurrencies_index = self.bins.index(selected_realization)
					amplified_occurrencies[amplified_occurrencies_index] = amplified_occurrencies[amplified_occurrencies_index] + current_amplification_factor

			else: #Perform a simple 'rescaling' with a random-selected factor (between 1 and maximum_amplification_factor)
				# Tick amplification bias
				amplified_IS_Histogram.amplification_bias = False

				# Choose a superimposed amplification_factor
				amplification_factor = np.random.random_integers(minimum_amplification_factor,maximum_amplification_factor)

				# Apply amplification_factor
				total_amplified_n_of_events = self.n_of_events*amplification_factor
				amplified_discrete_realizations = []
				for realization in self.discrete_realizations:
					amplified_discrete_realizations = amplified_discrete_realizations + [realization]*amplification_factor
				amplified_occurrencies = []
				for occurrence in self.occurrencies:
					amplified_occurrencies.append(occurrence*amplification_factor)


			### Introduce PCR-SLIPPAGE BIAS for each occurrence (bin)
			if (slippage_bias == True):

				# Tick slippage bias
				amplified_IS_Histogram.slippage_bias = True

				# Inizialize Slippage Biases Dictionary: #Key: 'bin value'; #Item: related_change_in_occurrencies
				slippage_biases_dictionary = {}
				for bin in self.bins:
					slippage_biases_dictionary.update({str(bin):0})

				# State gain and loss in occurrencies for each bin; store infos in Slippage Biases Dictionary
				for bin, occurrency in zip(self.bins, amplified_occurrencies):

					# State how much occurrencies go away:
					# RATIONALE: only few sequences will change their mapping location because of a stutter (e.g. the stutter occurs in firsts nucleotides... 0to1%...)
					occurrencies_variation = np.random.random_integers(int(occurrency*minimum_splippage_percentage*0.01),int(occurrency*maximum_splippage_percentage*0.01))

					# State 'where this stuff goes'
					#key = str(random.choice(self.bins)) # Old code, working but conceptually poor, even more if you increase span!
					# New code - RATIONALE: 'A stochastic model of the processes in PCR based amplification [...]', Jos Weusten, Jos Herbergs 2012
					#                        doi:10.1016/j.fsigen.2011.01.003
					min_p_two_stutters = 0.0
					max_p_two_stutters = 1.0
					min_p_one_stutter = 2.5
					max_p_one_stutter = 25.0
					p_two_stutter = random.triangular(min_p_two_stutters, max_p_two_stutters) / 100.0
					p_one_stutter = random.triangular(min_p_one_stutter, max_p_one_stutter) /100.0
					p_no_stutter = 1.0 - p_one_stutter - p_two_stutter
					direction = random.choice([-1,1])
					sample = [bin]*int(round(p_no_stutter,3)*1000) + [bin+direction]*int(round(p_one_stutter,3)*1000) + [bin+(2*direction)]*int(round(p_two_stutter,3)*1000)
					choice = random.choice(sample)
					key = str(choice)
					if (choice not in self.bins):
						key = str(bin)

					# Update Slippage Biases Dictionary
					slippage_biases_dictionary[str(bin)] = slippage_biases_dictionary[str(bin)] - occurrencies_variation # Current bin loose...
					slippage_biases_dictionary[key] = slippage_biases_dictionary[key] + occurrencies_variation # ...what a random-selected bin gains!

				# Update amplified_occurrencies taking advantage of slippage_biases_dictionary (total_amplified_n_of_events is conserved) and amplified_discrete_realizations (needed by GOF Test)
				amplified_discrete_realizations = []
				for bin in self.bins:
					index = self.bins.index(bin)
					amplified_occurrencies[index] = amplified_occurrencies[index] + slippage_biases_dictionary[str(bin)]
					amplified_discrete_realizations = amplified_discrete_realizations + [bin]*amplified_occurrencies[index]

			else:
				# Tick slippage bias
				amplified_IS_Histogram.slippage_bias = False


			### Calculate amplified_frequencies (abundances)
			amplified_frequencies = [] ### VARIABLE TO RETURN
			for amplified_occurrence in amplified_occurrencies:
				amplified_frequencies.append(float(amplified_occurrence)/float(total_amplified_n_of_events))


			### Calculate p_value_post_amplification and another_GOF_indicator_post_amplification
			p_value_post_amplification = None ### VARIABLE TO RETURN
			another_GOF_indicator_post_amplification = None ### VARIABLE TO RETURN
			if (self.source_distribution == 'gauss'):
				# Calculate p-value of Kolmogorov-Smirnov test (with respect to initial distribution settings)
				expected_value = self.expected_value
				st_dev = self.distribution_parameters['st_dev']
				standardized_data = [ float(r - expected_value+0.5)/float(st_dev) for r in amplified_discrete_realizations]
				KS_test, p_value = stats.kstest(standardized_data, 'norm')
				p_value_post_amplification = p_value
				another_GOF_indicator_post_amplification = KS_test
			else:
				print "\n\n[WARNING] 'Amplify()'' method of IS_Histogram Class, up to now, can't provide either 'p_value_post_amplification' or 'another_GOF_indicator_post_amplification' for {0} distribution...!".format(self.source_distribution)
				print "          By now, 'None'(s) have been returned; conversely, all the other attributes are correct.\n"

			### Update 'amplified' attribute
			amplified_IS_Histogram.n_of_mapped_sequencies = total_amplified_n_of_events
			amplified_IS_Histogram.sequence_count = amplified_occurrencies
			amplified_IS_Histogram.abundances = amplified_frequencies
			amplified_IS_Histogram.p_value_post_amplification = p_value_post_amplification
			amplified_IS_Histogram.another_GOF_indicator_post_amplification = another_GOF_indicator_post_amplification #KS_test if source_distribution='gauss'

			### Return the amplified IS_Histogram
			return amplified_IS_Histogram

		else:
			print "\n\n[WARNING] You are trying to amplify() an IS_Histogram object that seems to be already amplified...!"
			print "          Because of this, the same object has been returned by amplify() method.\n"

			### Resturn itself
			return self




	def bar_plot (self, bar_width = 0.8, color = 'blue', title = 'auto', show = False, save = False, name = 'auto', path = 'current_location', id_num = 0):

		### Open figure
		plt.figure()


		### Set figure details

		# Choose Title
		figure_title = ""
		if (title == 'auto'):
			if (self.source_distribution == 'gauss'):
				figure_title = self.source_distribution + " (N={0}, mean={1}, std={2}, p={3})".format(str(self.n_of_events), str(self.expected_value), str(self.distribution_parameters['st_dev']), str(self.p_value)[:4])
			else:
				figure_title = self.source_distribution + " (N={0}, expectation={1}, next_moment={2}, p={3})".format(str(self.n_of_events), str(self.expected_value), str(self.distribution_parameters), str(self.p_value)[:4])
		else:
			figure_title = title
		# Case of Amplified IS_Histogram
		if (self.amplified == True):
			figure_title = figure_title + "\n[Amplified to MS = {0}, AmpBias = {2}, SlipBias = {3}; p = {1}]".format(str(self.n_of_mapped_sequencies), str(self.p_value_post_amplification)[:4], str(self.amplification_bias)[0], str(self.slippage_bias)[0])

		# Insert Title
		plt.title(figure_title)

		# Axes Labels
		plt.xlabel('bp')
		y_label = "#N of Occurrencies"
		if (self.amplified == True):
			y_label = "Sequence_Count"
		plt.ylabel(y_label)


		### Create histogram in figure
		if (self.amplified == True):
			plt.bar(self.bins, self.sequence_count, align='center', width=bar_width, hold=True, facecolor=color)
		else:
			plt.bar(self.bins, self.occurrencies, align='center', width=bar_width, hold=True, facecolor=color)


		### Show
		if (show == True):
			plt.show()

		### Save
		if (save == True):
			
			# Name
			file_name = ""
			if (name == 'auto'):
				if (self.source_distribution == 'gauss'):
					file_name =  "#" + str(id_num) + "_Hist_"+self.source_distribution+"_[N={0}-mu={1}-sigma={2}]_p={3}".format(str(self.n_of_events), str(self.expected_value), str(self.distribution_parameters['st_dev']), str(self.p_value)[:4])
				else:
					next_moment = ""
					if ((distribution_parameters is dict) and (distribution_parameters != {})):
						next_moment = "&"+str(self.distribution_parameters).split("': ")[1][:-1]
					file_name = "#" + str(id_num) + "_Hist_"+self.source_distribution+"_N{0}_[{1}-{2}]_p={3}".format(str(self.n_of_events), str(self.expected_value), next_moment, str(self.p_value)[:4])
				# Case of Amplified IS_Histogram
				if (self.amplified == True):
					file_name = file_name + "_[Amplified_MS{0}_ab{1}_sb{2}]".format(str(self.n_of_mapped_sequencies), str(self.amplification_bias)[0], str(self.slippage_bias)[0])
				file_name = file_name + ".pdf"
			else:
				file_name = name


			# Path
			file_path_and_name=""
			if (path == 'current_location'):
				file_path_and_name=os.path.normpath(os.path.join(os.getcwd(), file_name))
			else:
				file_path_and_name=os.path.normpath(path, file_name)

			# Create file
			plt.savefig(file_path_and_name)

		### Nothing to return ###
		return 0