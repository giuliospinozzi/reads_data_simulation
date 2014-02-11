###Import Module(s)#####
import dataset_module
########################



### Set-up variables ##############################
source_distribution = 'gauss'
distribution_parameters = {'gauss':{'st_dev':1.0}}
span = 7

N_of_Simulation = 2
N_IS = 100
n_of_events_per_IS = 10

amplification_bias = True
slippage_bias = True
##################################################

### Initialize variables #########################
List_of_Simulation_RUNs = []
##################################################




### Starting print
print "\n\n[START]\tTest Environment!"

### LOOP over N_of_Simulation
for i in range(1,N_of_Simulation+1):

	# Run a simulation
	print "\n\n\t### Running simulation {0} ... ".format(str(i)),
	Simulation_RUN_object = dataset_module.GENERATE_SIMULATION_RUN (n_of_events_per_IS = n_of_events_per_IS, amplification_bias = amplification_bias, slippage_bias = slippage_bias, N_IS = N_IS, source_distribution = source_distribution, distribution_parameters = distribution_parameters, span = span)
	List_of_Simulation_RUNs.append(Simulation_RUN_object)
	print "Done!"

	# Create bed(s) and summary file related to the simulation just run
	print "\n\t\t* Creating related *.bed file and summary ... ",
	chromosome = str(i) ### REMEMBER TO CHANGE CHROMOSOME EACH LOOP!!
	bedFile_name_and_path, summary_file_path_and_name = Simulation_RUN_object.generate_bedFile(chromosome)
	print "Done!"
	print "\n\t\t* Creating a reference *.bed file ... ",
	reference_bedFile_name_and_path = Simulation_RUN_object.generate_reference_bedFile(chromosome)
	print "Done!"

	# Create Association File(s)
	print "\n\t\t* Creating Associations File(s) ... ",
	assFile_path_and_name, REFassFile_path_and_name = dataset_module.generateAssociationFiles(List_of_Simulation_RUNs)
	print "Done!"

	n_of_events_per_IS += 1

### End print
print "\n\n[QUIT]\tBye!\n\n"