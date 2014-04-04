###Header################################################
header = """

+------------------------------------------------------+
 Module: classes_for_dataset_generation
 Author: Stefano Brasca
 Date:  February 10th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains classes used in dataset
    generation and DB export (dataset_module.py)
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import os
import uuid
from datetime import datetime
import copy
import sys
###############################

###Import Module(s)#####################
#import yyyy
########################################






class Simulation_RUN:

	def __init__ (self, N_IS, source_distribution, expected_value, distribution_parameters, span, n_of_events_per_IS, IS_list, n_IS_Histogram_draining_out, total_n_of_realization_out, average_p_value, p_value_st_dev, pre_amplification_diversity,
		Amplified_IS_list, amplification_bias, minimum_amplification_factor, maximum_amplification_factor, slippage_bias, minimum_splippage_percentage, maximum_splippage_percentage,
		total_n_of_sequencies, mean_amplification_factor, average_p_value_post_amplification, p_value_st_dev_post_amplification, post_amplification_diversity):


		### Input parameters, IS Histograms list and summary

		self.N_IS = N_IS
		self.source_distribution = source_distribution
		self.expected_value = expected_value
		self.distribution_parameters = distribution_parameters
		self.span = span
		self.n_of_events_per_IS = n_of_events_per_IS

		self.IS_list = IS_list

		self.n_IS_Histogram_draining_out = n_IS_Histogram_draining_out
		self.total_n_of_realization_out = total_n_of_realization_out
		self.average_p_value = average_p_value
		self.p_value_st_dev = p_value_st_dev
		self.pre_amplification_diversity = pre_amplification_diversity


		### Post PCR: Amplified IS Histograms list, amplification settings and summary

		self.Amplified_IS_list = Amplified_IS_list

		self.amplification_bias = amplification_bias
		self.minimum_amplification_factor = None
		self.maximum_amplification_factor = None
		if (self.amplification_bias == True):
			self.minimum_amplification_factor = minimum_amplification_factor
			self.maximum_amplification_factor = maximum_amplification_factor

		self.slippage_bias = slippage_bias
		self.minimum_splippage_percentage = None
		self.maximum_splippage_percentage = None
		if (self.amplification_bias == True):
			self.minimum_splippage_percentage = minimum_splippage_percentage
			self.maximum_splippage_percentage = maximum_splippage_percentage


		### Results Report
		self.total_n_of_sequencies = total_n_of_sequencies
		self.mean_amplification_factor = float(self.total_n_of_sequencies) / float((self.N_IS*self.n_of_events_per_IS) - self.total_n_of_realization_out)
		self.average_p_value_post_amplification = average_p_value_post_amplification
		self.p_value_st_dev_post_amplification = p_value_st_dev_post_amplification
		self.post_amplification_diversity = post_amplification_diversity


		### Related Files

		self.bedFile = False
		self.bedFile_path_and_name = None

		self.reference_bedFile = False
		self.reference_bedFile_path_and_name = None

		self.summaryFile = False
		self.summaryFile_path_and_name = None

		self.associationFile = False
		self.associationFile_path_and_name = None

		self.reference_associationFile = False
		self.reference_associationFile_path_and_name = None





	def generate_bedFile (self, chromosome, IS_distance = 100, group = 1, group_internal_distance = 2, bedFile_name = 'auto', bedFile_path = 'current_location', summary_file = True, summary_file_name = 'auto', summary_file_path = 'current_location'):

		# group = 1 (min) -> isolated IS, separated by a distance of IS_distance
		# group > 1 -> ensembles of n=group IS, internally separated by a distance of group_internal_distance, and mutually separated by a distance of IS_distance

		### Args control
		if ((IS_distance < 100) or (type(IS_distance) is not int)):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'IS_distance' argument should be an integer >= 100; your input: '{0}'.".format(str(IS_distance))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
		if ((group < 1) or (type(group) is not int)):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'group' argument must be an integer >= 1; your input: '{0}'.".format(str(group))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
		if ((group != 1) and ((group_internal_distance < 1) or (type(group_internal_distance) is not int))):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'group_internal_distance' argument must be an integer >= 1; your input: '{0}'.".format(str(group_internal_distance))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

		### Prepare summary file name and path
		summary_file_path_and_name = None

		if (summary_file == False):
			summary_file_name = None
		
		else:
			
			if (summary_file_name == 'auto'):
				summary_file_name = "bedFiles_summary.tsv"
			else:
				summary_file_name = str(summary_file_name)

			if (summary_file_path == 'current_location'):
				summary_file_path_and_name = os.path.normpath(os.path.join(os.getcwd(), summary_file_name))
			else:
				summary_file_path_and_name = os.path.normpath(str(summary_file_path), summary_file_name)

		### Prepare bedFile name and path
		bedFile_path_and_name = ""

		if (bedFile_name == 'auto'):
			if (group == 1):
				bedFile_name = "chr{0}_{7}x{1}IS_shift{6}_{2}reads_{3}_ab{4}_sb{5}.bed".format(str(chromosome), str(self.N_IS), str(self.total_n_of_sequencies), str(self.source_distribution), str(self.amplification_bias), str(self.slippage_bias), str(IS_distance), str(self.n_of_events_per_IS))
			else:
				bedFile_name = "chr{0}_{7}x{1}IS_group{8}_dist{9}_shift{6}_{2}reads_{3}_ab{4}_sb{5}.bed".format(str(chromosome), str(self.N_IS), str(self.total_n_of_sequencies), str(self.source_distribution), str(self.amplification_bias), str(self.slippage_bias), str(IS_distance), str(self.n_of_events_per_IS), str(group), str(group_internal_distance))
		else:
			bedFile_name = str(bedFile_name)

		if (bedFile_path == 'current_location'):
			bedFile_path_and_name = os.path.normpath(os.path.join(os.getcwd(), bedFile_name))
		else:
			bedFile_path_and_name = os.path.normpath(str(bedFile_path), bedFile_name)


		### Write summary file
		with open(summary_file_path_and_name, 'a') as summary_file_IO:

			summary_file_IO.write ('BedFile:\t' + bedFile_path_and_name + '\n\n')

			summary_file_IO.write ( 'N_IS:\t' + str(self.N_IS) + '\n')
			summary_file_IO.write ( 'source_distribution:\t' + str(self.source_distribution) + '\n')
			summary_file_IO.write ( 'expected_value:\t' + str(self.expected_value) + '\n')
			summary_file_IO.write ( 'distribution_parameters:\t' + str(self.distribution_parameters) + '\n')
			summary_file_IO.write ( 'span:\t' + str(self.span) + '\n')
			summary_file_IO.write ( 'n_of_events_per_IS:\t' + str(self.n_of_events_per_IS) + '\n')
			if (group > 1):
				summary_file_IO.write ( 'grouping:\t' + 'True' + '\n')
				summary_file_IO.write ( 'n_of_IS_per_group:\t' + str(group) + '\n')
				summary_file_IO.write ( 'distance_between_IS_in_each_group:\t' + str(group_internal_distance) + '\n')
				summary_file_IO.write ( 'distance_between_groups:\t' + str(IS_distance) + '\n')
			else:
				summary_file_IO.write ( 'grouping:\t' + 'False' + '\n')
				summary_file_IO.write ( 'distance_between_ISs:\t' + str(IS_distance) + '\n')
			summary_file_IO.write ( 'n_IS_Histogram_draining_out:\t' + str(self.n_IS_Histogram_draining_out) +'\n')
			summary_file_IO.write ( 'total_n_of_realization_out:\t' + str(self.total_n_of_realization_out) + '\n')
			summary_file_IO.write ( 'average_p_value:\t' + str(self.average_p_value) + '\t+/-\t'+ str(self.p_value_st_dev) + '\n')
			summary_file_IO.write ( 'pre_amplification_diversity (%):\t' + str(self.pre_amplification_diversity*100.0) + '\n\n')

			summary_file_IO.write ( 'amplification_bias:\t' + str(self.amplification_bias) + '\t(min={0}, max={1})'.format(str(self.minimum_amplification_factor), str(self.maximum_amplification_factor)) +'\n')
			summary_file_IO.write ( 'slippage_bias:\t' + str(self.slippage_bias) + '\t(min={0}%, max={1}%)'.format(str(self.minimum_splippage_percentage), str(self.maximum_splippage_percentage)) +'\n')
			summary_file_IO.write ( 'post_amplification_diversity (%):\t' + str(self.post_amplification_diversity*100.0) + '\n\n')

			summary_file_IO.write ( 'total_n_of_sequencies:\t' + str(self.total_n_of_sequencies) + '\n')
			summary_file_IO.write ( 'mean_amplification_factor:\t' + str(self.mean_amplification_factor) + '\n')
			summary_file_IO.write ( 'average_p_value_post_amplification:\t' + str(self.average_p_value_post_amplification) + '\t+/-\t'+ str(self.p_value_st_dev_post_amplification) + '\n\n')

			summary_file_IO.write ('\n***\t***\t***\t***\t***\t***\n\n\n')
			summary_file_IO.close()


		### Write bedFile
		with open(bedFile_path_and_name, 'w') as bedFile_IO:
			bedFile_rows = []
			shift = -1 * IS_distance
			row = None
			group_control = copy.deepcopy(group)

			for Amplified_IS in self.Amplified_IS_list:
				group_control += 1
				if (group_control > group):
					shift = shift + IS_distance
					group_control = 1
				else:
					shift = shift + group_internal_distance
				for bin, sequence_count in zip(Amplified_IS.bins, Amplified_IS.sequence_count):
					for i in range(0,sequence_count): # Rows with 'zero' SC are skipped!
						row = "chr{0}\t{1}\t{2}\t{3}\t60\t+\t00M".format(str(chromosome), str(bin+shift), str(bin+shift), str(uuid.uuid4()))
						bedFile_rows.append(row)

			bedFile_IO.write("\n".join(bedFile_rows))
			bedFile_IO.close()


		### Update attributes
		self.bedFile = True
		self.bedFile_path_and_name = bedFile_path_and_name
		self.summaryFile = summary_file #Boolean
		self.summaryFile_path_and_name = summary_file_path_and_name

		### Return complete (name and) path of created files (None if file has not be created)
		return bedFile_path_and_name, summary_file_path_and_name




	def generate_reference_bedFile (self, chromosome, IS_distance = 100, group = 1, group_internal_distance = 2, bedFile_name = 'auto', bedFile_path = 'current_location'):

		# group = 1 (min) -> isolated IS, separated by a distance of IS_distance
		# group > 1 -> ensembles of n=group IS, internally separated by a distance of group_internal_distance, and mutually separated by a distance of IS_distance

		### Args control
		if ((IS_distance < 100) or (type(IS_distance) is not int)):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'IS_distance' argument should be an integer >= 100; your input: '{0}'.".format(str(IS_distance))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
		if ((group < 1) or (type(group) is not int)):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'group' argument must be an integer >= 1; your input: '{0}'.".format(str(group))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")
		if  ((group != 1) and ((group_internal_distance < 1) or (type(group_internal_distance) is not int))):
			error_message = "[ERROR]\tgenerate_bedFile bad calling: 'group_internal_distance' argument must be an integer >= 1; your input: '{0}'.".format(str(group_internal_distance))
			sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

		### Prepare reference bedFile name and path
		bedFile_path_and_name = ""

		if (bedFile_name == 'auto'):
			if (group == 1):
				bedFile_name = "chr{0}_{7}x{1}IS_shift{6}_{2}reads_{3}_ab{4}_sb{5}_REFERENCE.bed".format(str(chromosome), str(self.N_IS), str(self.total_n_of_sequencies), str(self.source_distribution), str(self.amplification_bias), str(self.slippage_bias), str(IS_distance), str(self.n_of_events_per_IS))
			else:
				bedFile_name = "chr{0}_{7}x{1}IS_group{8}_dist{9}_shift{6}_{2}reads_{3}_ab{4}_sb{5}_REFERENCE.bed".format(str(chromosome), str(self.N_IS), str(self.total_n_of_sequencies), str(self.source_distribution), str(self.amplification_bias), str(self.slippage_bias), str(IS_distance), str(self.n_of_events_per_IS), str(group), str(group_internal_distance))
		else:
			bedFile_name = str(bedFile_name)

		if (bedFile_path == 'current_location'):
			bedFile_path_and_name = os.path.normpath(os.path.join(os.getcwd(), bedFile_name))
		else:
			bedFile_path_and_name = os.path.normpath(str(bedFile_path), bedFile_name)


		### Write bedFile
		with open(bedFile_path_and_name, 'w') as bedFile_IO:
			bedFile_rows = []
			shift = -1 * IS_distance
			row = None
			group_control = copy.deepcopy(group)

			for Amplified_IS in self.Amplified_IS_list:
				group_control += 1
				if (group_control > group):
					shift = shift + IS_distance
					group_control = 1
				else:
					shift = shift + group_internal_distance
				for i in range(0,Amplified_IS.n_of_mapped_sequencies):
					row = "chr{0}\t{1}\t{2}\t{3}\t60\t+\t00M".format(str(chromosome), str(Amplified_IS.expected_value+shift), str(Amplified_IS.expected_value+shift), str(uuid.uuid4()))
					bedFile_rows.append(row)

			bedFile_IO.write("\n".join(bedFile_rows))
			bedFile_IO.close()

		### Update attributes
		self.reference_bedFile = True
		self.reference_bedFile_path_and_name = bedFile_path_and_name

		### Return complete (name and) path of created files (None if file has not be created)
		return bedFile_path_and_name