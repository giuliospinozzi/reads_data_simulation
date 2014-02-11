###Header################################################
header = """

+------------------------------------------------------+
 Module: dataset_module
 Author: Stefano Brasca
 Date:  February 10th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions used in dataset
    generation and DB export
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import numpy as np
import os
###############################

###Import Module(s)#####################
import classes_for_dataset_generation
import simulation_module
########################################






def GENERATE_SIMULATION_RUN (n_of_events_per_IS, amplification_bias, slippage_bias, N_IS = 2000, source_distribution = 'gauss', distribution_parameters = {'gauss':{'st_dev':1.0}}, span = 7):

	'''
	*** The best way to get an Simulation_RUN Class Instance from few input parameters ***

	INPUT: - [...]

	OUTPUT: Simulation_RUN object

	'''

	### Set SIMULATION RUN DEFAULTS
	
	# Amplification Bias settings (if given True)
	minimum_amplification_factor=0
	maximum_amplification_factor=1000

	# Slippage Bias settings (if given True)
	minimum_splippage_percentage=0
	maximum_splippage_percentage=1


	#################################################
	# TO DO:                                        #
	# - At least a little piece of documentation!!! #
	# - Some controls (and related warnings!)       #
	#################################################

	### Initialize variables to store results

	IS_Histogram_list = []
	IS_Histogram_frequencies_set = set()

	Amplified_IS_Histogram_list = []
	Amplified_IS_Histogram_frequencies_set = set()

	total_n_of_sequencies = 0
	n_IS_Histogram_draining_out = 0
	total_n_of_realization_out = 0


	# Initialize temp variables

	p_value_list = []
	p_value_post_amplification_list = []
	realization_count = 0
	sequence_mapped_count = 0

	for i in range(1, N_IS+1):

		# Generate IS Hist -> IS_Histogram_object
		IS_Histogram_object = simulation_module.GENERATE_IS_HISTOGRAM(source_distribution, distribution_parameters, span, n_of_events_per_IS)
		IS_Histogram_list.append(IS_Histogram_object)
		IS_Histogram_frequencies_set.add(tuple(IS_Histogram_object.frequencies))

		# Other tasks on IS_Histogram_object
		if (IS_Histogram_object.discrete_realizations_beyond_edges != {}):
			n_IS_Histogram_draining_out += 1
			for distance, n_realization_out in IS_Histogram_object.discrete_realizations_beyond_edges.iteritems():
				total_n_of_realization_out += n_realization_out
		p_value_list.append(IS_Histogram_object.p_value)
		realization_count += len(IS_Histogram_object.discrete_realizations)


		# Amplify IS Hist -> Amplified_IS_Histogram_object
		Amplified_IS_Histogram_object = IS_Histogram_object.amplify(minimum_amplification_factor=minimum_amplification_factor, maximum_amplification_factor=maximum_amplification_factor, amplification_bias=amplification_bias, slippage_bias=slippage_bias, minimum_splippage_percentage=minimum_splippage_percentage, maximum_splippage_percentage=maximum_splippage_percentage)
		Amplified_IS_Histogram_list.append(Amplified_IS_Histogram_object)
		Amplified_IS_Histogram_frequencies_set.add(tuple(Amplified_IS_Histogram_object.abundances))

		# Other tasks on Amplified_IS_Histogram_object
		total_n_of_sequencies += Amplified_IS_Histogram_object.n_of_mapped_sequencies
		p_value_post_amplification_list.append(Amplified_IS_Histogram_object.p_value_post_amplification)
		sequence_mapped_count += Amplified_IS_Histogram_object.n_of_mapped_sequencies


	# Finalize results
	average_p_value = np.mean(p_value_list)
	p_value_st_dev = np.std(p_value_list)

	average_p_value_post_amplification = np.mean(p_value_post_amplification_list)
	p_value_st_dev_post_amplification = np.std(p_value_post_amplification_list)

	pre_amplification_diversity = float(len(IS_Histogram_frequencies_set)) / float(len(IS_Histogram_list))
	post_amplification_diversity = float(len(Amplified_IS_Histogram_frequencies_set)) / float(len(Amplified_IS_Histogram_list))

	mean_amplification_factor = float(sequence_mapped_count) / float(realization_count)


	# CREATE FINAL OBJECT TO RETURN #
	return classes_for_dataset_generation.Simulation_RUN(N_IS, source_distribution, distribution_parameters, span, n_of_events_per_IS, IS_Histogram_list, n_IS_Histogram_draining_out, total_n_of_realization_out, average_p_value, p_value_st_dev, pre_amplification_diversity,
		Amplified_IS_Histogram_list, amplification_bias, minimum_amplification_factor, maximum_amplification_factor, slippage_bias, minimum_splippage_percentage, maximum_splippage_percentage,
		total_n_of_sequencies, mean_amplification_factor, average_p_value_post_amplification, p_value_st_dev_post_amplification, post_amplification_diversity)




def generateAssociationFiles (List_of_Simulation_RUNs, assFile_name = 'auto', assFile_path = 'current_location'):

		### Check for *.bed files
		for Simulation_RUN in List_of_Simulation_RUNs:
			if (Simulation_RUN.bedFile == False):
				error_message = "[ERROR]\tgenerateAssociationFiles() bad calling: requested *.bed files are not available!"
				sys.exit("\n\n" + error_message + "\n\t[QUIT]\n\n")

		### Prepare assFile name and path
		assFile_path_and_name = ""

		if (assFile_name == 'auto'):
			assFile_name = "AssociationFile.tsv"
		else:
			assFile_name = str(assFile_name)

		if (assFile_path == 'current_location'):
			assFile_path_and_name = os.path.normpath(os.path.join(os.getcwd(), assFile_name))
		else:
			assFile_path_and_name = os.path.normpath(str(assFile_path), assFile_name)

		### Check for *_REFERENCE.bed files
		reference_files = True
		for Simulation_RUN in List_of_Simulation_RUNs:
			if (Simulation_RUN.reference_bedFile == False):
				reference_files = False
		if (reference_files == False):
			message = "[WARNING - generateAssociationFiles() call]\t*_REFERENCE.bed files are not available"
			print "\n\n" + message + "\n\t[CONTINUE without processing REFERENCE DATA]\n\n"

		### Prepare REFassFile name and path
		REFassFile_path_and_name = None
		if (reference_files == True):
			REFassFile_name = "REFERENCE_" + assFile_name
			if (assFile_path == 'current_location'):
				REFassFile_path_and_name = os.path.normpath(os.path.join(os.getcwd(), REFassFile_name))
			else:
				REFassFile_path_and_name = os.path.normpath(str(assFile_path), REFassFile_name)

		### LOOP OVER List_of_Simulation_RUNs: computing assFile_rows and REFassFile_rows
		assFile_rows = []
		REFassFile_rows = []
		void = "NULL"
		for Simulation_RUN in List_of_Simulation_RUNs:

			### assFile
			bedFile_name = os.path.basename(Simulation_RUN.bedFile_path_and_name) # bedFile_name in place of 'barcode join', 'LAM_id' and 'complete name of LAM'
			treatment = str(Simulation_RUN.n_of_events_per_IS) # Simulation_RUN.n_of_events_per_IS in place of treatment (time point)
			tissue = str(Simulation_RUN.amplification_bias) # Simulation_RUN.amplification_bias in place of tissue
			sample = str(Simulation_RUN.slippage_bias)
			ass_row = bedFile_name + '\t' + bedFile_name + '\t' + tissue + '\t' + void + '\t' + treatment + '\t' + bedFile_name + '\t' + bedFile_name + '\t' + sample +  '\t' + void + '\t' + void + '\t' + void
			assFile_rows.append(ass_row)

			### REFassFile
			if (reference_files == True):
				REFbedFile_name = os.path.basename(Simulation_RUN.reference_bedFile_path_and_name) # REFbedFile_name in place of 'barcode join', 'LAM_id' and 'complete name of LAM'
				REFass_row = REFbedFile_name + '\t' + REFbedFile_name + '\t' + tissue + '\t' + void + '\t' + treatment + '\t' + REFbedFile_name + '\t' + REFbedFile_name + '\t' + sample +  '\t' + void + '\t' + void + '\t' + void
				REFassFile_rows.append(REFass_row)


		### Write files

		### assFile
		with open(assFile_path_and_name, 'w') as assFile_IO:
			assFile_IO.write("\n".join(assFile_rows))
			assFile_IO.close()

		### REFassFile
		if (reference_files == True):
			with open(REFassFile_path_and_name, 'w') as REFassFile_IO:
				REFassFile_IO.write("\n".join(REFassFile_rows))
				REFassFile_IO.close()

		### Return complete (name and) path of created files (None if file has not be created)
		return assFile_path_and_name, REFassFile_path_and_name
