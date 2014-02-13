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
span = 15 # N of bin per IS; expected value in the middle (7) for symmetric distribution; realization
		  # beyond edges are discarded but tracked


### Simulation Runs Parameters
starting_n_of_events_per_IS = 1 # N of discrete realization per IS Histogram
ending_n_of_events_per_IS = 20 # N of discrete realization per IS Histogram
n_of_simulation_RUN = ending_n_of_events_per_IS - starting_n_of_events_per_IS + 1

n_IS_per_simulation_RUN = 1000 # N of IS Histogram built for each n in 
							   # range(starting_n_of_events_per_IS, ending_n_of_events_per_IS +1)

amplification_bias = True # Fine tuning of 'how to' is performable in code: amplify() method of IS_Histogram Class
slippage_bias = True # Fine tuning of 'how to' is performable in code: amplify() method of IS_Histogram Class


### Export Data to DB
export_data_to_DB = True
PYTHON_TOOL = "dbimport_redundantiss_from_bed.v2.py"

DB = "local"

dbschema = "simulation_data"

dbtable = "alldata"
reference_dbtable = "reference"

patient = "Simulation"
pool = "Stefano"
tag = None ### This variable, each time, must host bedFile name!!!
		     # Let be 'None' to do it automatically

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
	print "\t\t\t\t  - DB choice: ", DB
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
j = 1
for i in range(starting_n_of_events_per_IS, ending_n_of_events_per_IS+1): # Change N of Events per IS each loop (+1)

	for aB in amplification_bias_switcher: # Switch Amplification Bias

		for sB in slippage_bias_switcher: # Switch Slippage Bias

			# Run a simulation
			print "\n\n{0}\t### Running simulation {1} of {5} ({2} event per IS, AB {3}, SB {4}) ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(j), str(i), str(aB), str(sB), str(n_of_simulation_RUN*len(amplification_bias_switcher)*len(slippage_bias_switcher))),
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
			j += 1
# Now results are in List_of_Simulation_RUNs


### Create Association Files (one for simulated data and another for reference data)
print "\n\n\n{0}\t[ASSOCIATION FILES GENERATION] Processing ... \n".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),
associationFile_path_and_name, REFassociationFile_path_and_name = dataset_module.generateAssociationFiles(List_of_Simulation_RUNs)
print "Done!"
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), associationFile_path_and_name
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), REFassociationFile_path_and_name



### Launch external app to write result on DB
print "\n\n\n{0}\t[EXPORT DATA TO DB] Processing ... \n".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),

# Loop over List_of_Simulation_RUNs
j = 1
for Simulation_RUN in List_of_Simulation_RUNs:
	print "\n{0}\t### Exporting simulation {1} of {2} ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(j), str(n_of_simulation_RUN*len(amplification_bias_switcher)*len(slippage_bias_switcher))),
	std_output, errors = dataset_module.exportDataToDB (PYTHON_TOOL, Simulation_RUN.bedFile_path_and_name, Simulation_RUN.associationFile_path_and_name, patient, pool, tag, DB, dbschema, dbtable)
	print "Done!"
	print "{0}\t### Exporting related reference data ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),
	std_output_REF, errors_REF = dataset_module.exportDataToDB (PYTHON_TOOL, Simulation_RUN.reference_bedFile_path_and_name, Simulation_RUN.reference_associationFile_path_and_name, patient, pool, tag, DB, dbschema, reference_dbtable)
	print "Done!"
	j += 1

### Final
print "\n\n[QUIT]\tAll tasks performed. Bye!\n\n"

####################################################################################################################################################################
### END ############################################################################################################################################################
####################################################################################################################################################################


### Tip # Integration Analysis command template


# Classic
# python Integration_Analysis.py --dbDataset "simulation_data.alldata,simulation_data.reference" --columns treatment,tissue,sample --IS_method classic --collision --set_radius 0 --tsv --statistics

# Gauss
# python Integration_Analysis.py --dbDataset "simulation_data.alldata,simulation_data.reference" --columns treatment,tissue,sample --IS_method gauss --interaction_limit 3 --alpha 0.5  --collision --set_radius 0 --tsv --statistics

# Comments
# 1) treatment,tissue and sample  are the only discriminant in our datasets (respectively 'original N of gaussian Events per IS' (int), 'PCR Amplification Bias' (bool) and 'Polymerase Slippage Bias' (bool))
# 2) --collision --set_radius 0 underlines perfect matching with 'reference' data
# 3) default value for bushman_bp_rule in classic (3) and --interaction_limit 3 (=bushman_bp_rule) in gauss guarantee the same ensamble construction rule
# 4) --alpha 0.5 in gauss states '1sigma=1bp', that is the simulation criteria
# 5) If a phenomenon was regulated by a gaussian distribution like the one used by my algorithm, the 99.9 % of events would be under the comparison histogram used during the analysis

# Obscure point: why we modeled an integration event as coming from such a gaussian?!


### Further improvements

# better if amplification factor was extracted from 0 (now 1) to max (now 100), in order to account for sequences lost
# this change should be selective, and active only until at least a sequence remain in histogram!



