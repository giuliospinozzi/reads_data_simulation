###Requested Package(s) Import###
from time import gmtime, strftime
import os
#################################

###Import Module(s)#####
import NEW_dataset_module
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

amplification_bias = True # Fine tuning of default 'how to' is performable in code: amplify() method of IS_Histogram Class
						  # ACTUALLY IS 'GENERATE_SIMULATION_RUN' THAT FIX DEFAULT VALUE FOR SIMULATION (in code)
slippage_bias = True # Fine tuning of default 'how to' is performable in code: amplify() method of IS_Histogram Class
					 # ACTUALLY IS 'GENERATE_SIMULATION_RUN' THAT FIX DEFAULT VALUE FOR SIMULATION (in code)
### *.bed files features
IS_distance = 100 # 100 or more, distance between IS / IS-group
group = 2 # 1 or more, N of IS per group: 1 for standard simulations, more for disentanglement
group_internal_distance = 2 # 1 or more, distance between IS in the same group



### Export Data to DB
export_data_to_DB = True
PYTHON_TOOL = "dbimport_redundantiss_from_bed.v2.py"

DB = "local"
dbschema = "sequence_simulation"
dbtable = "disentanglement_allbiases"
### Setting up names ###
if (group > 1):
	dbtable = dbtable + "_group{0}_dist{1}".format(str(group), str(group_internal_distance))

###      #####       ###
reference_dbtable = dbtable + "_reference"

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

# ### Setting up loop variables
# amplification_bias_switcher = [False]
# slippage_bias_switcher = [False]
# if (amplification_bias == True):
# 	amplification_bias_switcher.append(True)
# if (slippage_bias == True):
# 	slippage_bias_switcher.append(True)

### TEMP MOD, to manually split parameters combination in different tables
amplification_bias_switcher = [amplification_bias]
slippage_bias_switcher = [slippage_bias]


### LOOP FOR SIMULATION RUNS
j = 1
for i in range(starting_n_of_events_per_IS, ending_n_of_events_per_IS+1): # Change N of Events per IS each loop (+1)

	for aB in amplification_bias_switcher: # Switch Amplification Bias

		for sB in slippage_bias_switcher: # Switch Slippage Bias

			# Run a simulation
			print "\n\n{0}\t### Running simulation {1} of {5} ({2} event per IS, AB {3}, SB {4}) ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(j), str(i), str(aB), str(sB), str(n_of_simulation_RUN*len(amplification_bias_switcher)*len(slippage_bias_switcher))),
			Simulation_RUN_object = NEW_dataset_module.GENERATE_SIMULATION_RUN (n_of_events_per_IS = i, amplification_bias = aB, slippage_bias = sB, N_IS = n_IS_per_simulation_RUN, source_distribution = source_distribution, distribution_parameters = distribution_parameters, span = span)
			List_of_Simulation_RUNs.append(Simulation_RUN_object)
			print "Done!"

			# Create bed(s) and summary file related to the simulation just run
			print "\n\t\t\t\t* Creating related *.bed file and summary ... ",
			chromosome = str(i) ### REMEMBER TO CHANGE CHROMOSOME EACH LOOP!!
			bedFile_name_and_path, summary_file_path_and_name = Simulation_RUN_object.generate_bedFile(chromosome, IS_distance = IS_distance, group = group, group_internal_distance = group_internal_distance)
			print "Done!"
			print "\t\t\t\t ", bedFile_name_and_path
			print "\n\t\t\t\t* Creating a reference *.bed file ... ",
			reference_bedFile_name_and_path = Simulation_RUN_object.generate_reference_bedFile(chromosome, IS_distance = IS_distance, group = group, group_internal_distance = group_internal_distance)
			print "Done!"
			print "\t\t\t\t ", reference_bedFile_name_and_path
			j += 1
# Now results are in List_of_Simulation_RUNs


### Create Association Files (one for simulated data and another for reference data)
print "\n\n\n{0}\t[ASSOCIATION FILES GENERATION] Processing ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),
unique_bedfile_name = dbtable + ".bed"
unique_reference_bedfile_name = reference_dbtable + ".bed"
associationFile_path_and_name, REFassociationFile_path_and_name = NEW_dataset_module.generateAssociationFiles(List_of_Simulation_RUNs, unique_bedfile_name, unique_reference_bedfile_name)
print "Done!"
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), associationFile_path_and_name
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), REFassociationFile_path_and_name


### Unify *.bed files
print "\n\n\n{0}\t[UNIFY *.BED FILES] Processing ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),

bedfiles_as_list = []
REFERENCE_bedfiles_as_list = []
fileList = os.listdir(".")
for filename in fileList:
	if (filename.split('.')[-1] == 'bed'):
		filePointer = open (filename, 'r')
		if ('REFERENCE' in filename):
			REFERENCE_bedfiles_as_list = REFERENCE_bedfiles_as_list + filePointer.readlines()
		else:
			bedfiles_as_list = bedfiles_as_list + filePointer.readlines()
		filePointer.close()

with open(unique_bedfile_name, 'w') as unique_bedfile:
	for line in bedfiles_as_list:
		if (('\n' not in line) and (line != bedfiles_as_list[-1])):
			unique_bedfile.write(line+'\n')
		else:
			unique_bedfile.write(line)
	unique_bedfile.close()

with open(unique_reference_bedfile_name, 'w') as unique_reference_bedfile:
	for line in REFERENCE_bedfiles_as_list:
		if (('\n' not in line) and (line != REFERENCE_bedfiles_as_list[-1])):
			unique_reference_bedfile.write(line+'\n')
		else:
			unique_reference_bedfile.write(line)
	unique_reference_bedfile.close()

print "done!"
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), unique_bedfile_name
print "{0}\t* ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), unique_reference_bedfile_name


### Launch external app to write result on DB
print "\n\n\n{0}\t[EXPORT DATA TO DB] Processing ... ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),
export_data_log_file = open("export_data_log_file.log","w")

unique_bedfile_name_and_path = os.path.normpath(os.path.join(os.getcwd(), unique_bedfile_name))
std_output, errors = NEW_dataset_module.exportDataToDB (PYTHON_TOOL, unique_bedfile_name_and_path, associationFile_path_and_name, patient, pool, tag, DB, dbschema, dbtable)
export_data_log_file.write("{0}\t### Exporting ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())) + unique_bedfile_name_and_path + " to " + dbschema + "-" + dbtable + "\n\n")
export_data_log_file.write("# std_output:\n")
export_data_log_file.write(std_output)
export_data_log_file.write("\n# errors:\n")
export_data_log_file.write(errors)
export_data_log_file.write("\n\nDone!\n\n\n\n")

unique_reference_bedfile_name_and_path = os.path.normpath(os.path.join(os.getcwd(), unique_reference_bedfile_name))
std_output_REF, errors_REF = NEW_dataset_module.exportDataToDB (PYTHON_TOOL, unique_reference_bedfile_name_and_path, REFassociationFile_path_and_name, patient, pool, tag, DB, dbschema, reference_dbtable)
export_data_log_file.write("{0}\t### Exporting ".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())) + unique_reference_bedfile_name_and_path + " to " + dbschema + "-" + reference_dbtable + "\n\n")
export_data_log_file.write("# std_output:\n")
export_data_log_file.write(std_output_REF)
export_data_log_file.write("\n# errors:\n")
export_data_log_file.write(errors_REF)
export_data_log_file.write("\n\nDone!")
print "Done!"

export_data_log_file.close()

### Final
print "\n\n[QUIT]\tAll tasks performed. Bye!\n\n"