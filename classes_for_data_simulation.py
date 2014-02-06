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
  	simulations
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import os
import copy
import random
from scipy import stats
import matplotlib.pyplot as plt
###############################

###Import Module(s)#####################
import simulation_module
########################################






class Histogram:

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
		self.n_of_mapped_sequencies = None # Integer; 'Amplified analogous' of n_of_events
		self.sequence_count = [None]*len(occurrencies) # List of integers, 'amplified analogous' of occurrencies, paired with bins
		self.abundances = [None]*len(frequencies) # List of float, 'amplified analogous' of frequencies, paired with bins
		self.p_value_post_amplification = None
		self.another_GOF_indicator_post_amplification = None #KS_test if source_distribution='gauss'




	def amplify (self, minimum_amplification_factor=1, maximum_amplification_factor=1000):

		# This method works only with unamplified Histograms
		if (self.amplified == False):

			### CAST
			if ((minimum_amplification_factor is not int) or (maximum_amplification_factor is not int)):
				minimum_amplification_factor = int(minimum_amplification_factor)
				maximum_amplification_factor = int(maximum_amplification_factor)
				print "\n[WARNING] 'Amplify()' method of Histogram Class supports only integer arguments. [Cast has been done]"

			### Create a copy of Histogram object
			amplified_Histogram = copy.deepcopy(self) ### OBJECT TO RETURN


			### Introduce PCR bias for each realization

			# Get a copy of self.discrete_realizations (to work on it)
			discrete_realizations = copy.deepcopy(self.discrete_realizations)

			# Inizialize results collectors
			amplified_occurrencies = [0]*len(self.occurrencies) ### VARIABLE TO RETURN
			amplified_discrete_realizations = []
			total_amplified_n_of_events = 0 ### VARIABLE TO RETURN


			# Loop: until discrete_realizations local variable is []
			while (discrete_realizations != []):

				# Randomly choose a realization
				selected_realization = random.choice(discrete_realizations)
				# remove it (like pop)
				discrete_realizations.remove(selected_realization)
				# Amplify: relate a number to selected_realization
				current_amplification_factor = random.randint(minimum_amplification_factor,maximum_amplification_factor+1)
				################################################################
				# A specific algorithm should be implemented soon              #
				# Something like:                                              #
				# amp_factor = simulation_module.get_amplification_factor(...) #
				################################################################
				# Update amplified_discrete_realizations
				amplified_discrete_realizations = amplified_discrete_realizations + [selected_realization]*current_amplification_factor
				# Update total_amplified_n_of_events
				total_amplified_n_of_events = total_amplified_n_of_events + current_amplification_factor
				# Sum the result in the related amplified_occurrency
				amplified_occurrencies_index = self.bins.index(selected_realization)
				amplified_occurrencies[amplified_occurrencies_index] = amplified_occurrencies[amplified_occurrencies_index] + current_amplification_factor

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
				standardized_data = [ float(r - expected_value)/float(st_dev) for r in amplified_discrete_realizations]
				KS_test, p_value = stats.kstest(standardized_data, 'norm')
				p_value_post_amplification = p_value
				another_GOF_indicator_post_amplification = KS_test
			else:
				print "\n[WARNING] 'Amplify()'' method of Histogram Class, up to now, can't provide either 'p_value_post_amplification' or 'another_GOF_indicator_post_amplification' for {0} distribution...!".format(self.source_distribution)
				print "          By now, 'None'(s) have been returned; conversely, all the other attributes are correct.\n"


			### Update 'amplified' attribute
			amplified_Histogram.amplified = True
			amplified_Histogram.n_of_mapped_sequencies = total_amplified_n_of_events
			amplified_Histogram.sequence_count = amplified_occurrencies
			amplified_Histogram.abundances = amplified_frequencies
			amplified_Histogram.p_value_post_amplification = p_value_post_amplification
			amplified_Histogram.another_GOF_indicator_post_amplification = another_GOF_indicator_post_amplification #KS_test if source_distribution='gauss'


			### Return the amplified Histogram
			return amplified_Histogram

		else:
			print "\n[WARNING] You are trying to amplify() an Histogram object that seems to be already amplified...!"
			print "          Because of this, the same object has been returned by amplify() method.\n"

			return self



### OBSOLETE ####################################################################################################################################################################
	def amplify_old (self, ):

		# This method works only with unamplified Histograms
		if (self.amplified == False):

			### Create a copy of Histogram object
			amplified_Histogram = copy.deepcopy(self)

			### Amplification factor
			amp_factor = 1000
			################################################################
			# A specific algorithm should be implemented soon              #
			# Something like:                                              #
			# amp_factor = simulation_module.get_amplification_factor(...) #
			################################################################

			### Amplify n_of_events -> Update n_of_mapped_sequencies
			amplified_Histogram.n_of_mapped_sequencies = (amplified_Histogram.n_of_events)*amp_factor


			### Introduce PCR biases bin by bin

			# Inizialize Biases Dictionary: #Key: 'bin value'; #Item: related_change_in_frequency
			biases_dictionary = {}
			for bin in amplified_Histogram.bins:
				biases_dictionary.update({str(bin):0})

			# State gain and loss in frequency for each bin; store infos in Biases Dictionary
			for bin, frequency in zip(amplified_Histogram.bins, amplified_Histogram.frequencies):

				### State 'how much frequency' goes away
				frequency_variation = random.triangular(0, (frequency*0.1))
				###################################################################
				# A specific algorithm should be implemented soon                 #
				# Something like:                                                 #
				# frequency_variation = simulation_module.get_frequency_loss(...) #
				###################################################################

				### State 'where this stuff goes'
				key = str(random.choice(amplified_Histogram.bins))
				###################################################
				# A specific algorithm should be implemented soon #
				# Something like:                                 #
				# key = simulation_module.get_bias_target (...)   #
				###################################################

				### Update Biases Dictionary
				biases_dictionary[str(bin)] = biases_dictionary[str(bin)] - frequency_variation # Current bin loose...
				biases_dictionary[key] = biases_dictionary[key] + frequency_variation # ...what a random-selected bin gains!


			### Set up sequence_count taking advantage of biases_dictionary
			for bin in amplified_Histogram.bins:
				index = amplified_Histogram.bins.index(bin)
				amplified_Histogram.sequence_count[index] = (amplified_Histogram.frequencies[index] + biases_dictionary[str(bin)]) * (amplified_Histogram.n_of_mapped_sequencies)

			### Update 'amplified' attribute
			amplified_Histogram.amplified = True

			### Return the amplified Histogram
			return amplified_Histogram

		else:
			print "\n[WARNING] You are trying to amplify() an Histogram object that seems to be already amplified...!"
			print " Because of this, the same object has been returned by amplify() method.\n"

			return self
####################################################################################################################################################################################



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
		# Case of Amplified Histogram
		if (self.amplified == True):
			figure_title = figure_title + "\n[Amplified to MS = {0}; p = {1}]".format(str(self.n_of_mapped_sequencies), str(self.p_value_post_amplification))

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
					file_name = "#" + str(id_num) + "_Hist_"+self.source_distribution+"_N{0}_[{1}{2}]_p={3}".format(str(self.n_of_events), str(self.expected_value), next_moment, str(self.p_value)[:4])
				# Case of Amplified Histogram
				if (self.amplified == True):
					file_name = file_name + "_[Amplified_MS{0}]".format(str(self.n_of_mapped_sequencies))
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