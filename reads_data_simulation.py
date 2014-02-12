###Requested Package(s) Import###
from time import gmtime, strftime
#################################

###Import Module(s)#####
import dataset_module
########################



### Set-up Paremaeters ###########################################################################################


### Distribution features
source_distribution = 'gauss' # The only available up to now
distribution_parameters = {'gauss':{'st_dev':1.0}} # Sigma = 1bp
												   # Due to discretized construction, the distribution
												   # 'core' is centered in the 3bp in the middle of span
span = 7 # N of bin per IS; expected value in the middle (3) for symmetric distribution; realization
		 # beyond edges are discarded but tracked


### Simulation Runs Parameters
starting_n_of_events_per_IS = 10 # N of discrete realization per IS Histogram
ending_n_of_events_per_IS = 20 # N of discrete realization per IS Histogram
n_of_simulation_RUN = ending_n_of_events_per_IS - starting_n_of_events_per_IS + 1

n_IS_per_simulation_RUN = 1000 # N of IS Histogram built for each n in 
							   # range(starting_n_of_events_per_IS, ending_n_of_events_per_IS +1)

amplification_bias = True # Fine tuning of 'how to' is performable in code: amplify() method of IS_Histogram Class
slippage_bias = True # Fine tuning of 'how to' is performable in code: amplify() method of IS_Histogram Class


### Export Data to DB
export_data_to_DB = False
db = "local"

dbschema = "SimulatedData"

dbtable = "All"
reference_dbtable = "Reference"

patient = "Simulation"
pool = "Stefano"
tag = None ### This variable, each time, must host bedFile name!!!

### OTHER PARAMETER IN 'GENERATE_SIMULATION_RUN' FUNCTION (dataset_module.py)
### (see ###Set SIMULATION RUN DEFAULTS)

##################################################################################################################



### Initialize global variables ##################
List_of_Simulation_RUNs = []
##################################################






####################################################################################################################################################################
### BEGIN ##########################################################################################################################################################
####################################################################################################################################################################

### Starting print
print "\n\n{0}\t[START]\tVECTOR INTEGRATION SITES READS DATA SIMULATION".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
print "\n\t\t\t\tBrief Summary:"
print "\t\t\t\t* N of Simulation to Run: ", n_of_simulation_RUN
print "\t\t\t\t* IS per Simulation: ", n_IS_per_simulation_RUN
print "\t\t\t\t* N of Events per IS: from {0} to {1} with integer steps of 1 each Run".format(str(starting_n_of_events_per_IS), str(ending_n_of_events_per_IS))
print "\t\t\t\t* Events source Distribution: ", source_distribution
print "\t\t\t\t  (available span {0}, mu {1}, sigma {2})".format(str(span), str(span/2), str(distribution_parameters[source_distribution]['st_dev']))
print "\t\t\t\t* Add PCR Amplification Bias: ", amplification_bias
print "\t\t\t\t* Add Polymerase Slippage Bias: ", slippage_bias
print "\t\t\t\t* EXPORT DATA TO DB: ", export_data_to_DB
if (export_data_to_DB == True):
	print "\t\t\t\t  - DB choice: ", db
	print "\t\t\t\t  - Generated Data: {0}.{1}".format(dbschema, dbtable)
	print "\t\t\t\t  - Reference Data: {0}.{1}".format(dbschema, reference_dbtable)
print "\n"

### Setting up loop variables
amplification_bias_switcher = [False]
slippage_bias_switcher = [False]
if (amplification_bias == True):
	amplification_bias_switcher.append(True)
if (slippage_bias == True):
	slippage_bias_switcher.append(True)


### LOOP FOR SIMULATION RUNS
for i in range(starting_n_of_events_per_IS, ending_n_of_events_per_IS+1): # Change N of Events per IS each loop (+1)

	for aB in amplification_bias_switcher: # Switch Amplification Bias

		for sB in slippage_bias_switcher: # Switch Slippage Bias

			# Run a simulation
			print "\n\n{0}\t### Running simulation {1} ({2} event per IS, AB {3}, SB {4}) ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(i-starting_n_of_events_per_IS+1), str(i), str(aB), str(sB)),
			Simulation_RUN_object = dataset_module.GENERATE_SIMULATION_RUN (n_of_events_per_IS = i, amplification_bias = aB, slippage_bias = sB, N_IS = n_IS_per_simulation_RUN, source_distribution = source_distribution, distribution_parameters = distribution_parameters, span = span)
			List_of_Simulation_RUNs.append(Simulation_RUN_object)
			print "Done!"

			# Create bed(s) and summary file related to the simulation just run
			print "\n\t\t\t\t* Creating related *.bed file and summary ... ",
			chromosome = str(i) ### REMEMBER TO CHANGE CHROMOSOME EACH LOOP!!
			bedFile_name_and_path, summary_file_path_and_name = Simulation_RUN_object.generate_bedFile(chromosome)
			print "Done!"
			print "\t\t\t\t ", bedFile_name_and_path
			print "\n\t\t\t\t* Creating a reference *.bed file ... ",
			reference_bedFile_name_and_path = Simulation_RUN_object.generate_reference_bedFile(chromosome)
			print "Done!"
			print "\t\t\t\t ", reference_bedFile_name_and_path