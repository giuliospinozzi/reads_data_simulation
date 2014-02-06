###Import Module(s)#####
import simulation_module
########################

### Set-up variables ##############################
source_distribution = 'gauss'
distribution_parameters = {'gauss':{'st_dev':1.0}}
span = 7
n_of_events = 15
n_of_trial = 100
p_value_threshold = 0.05

show_plot = True
save_plot = False
###################################################



### Starting print
print "\n\n[START]\tTest Environment!\n"

### Loop ##########################################################################################################################
ugly = 0
all_hist = set()
ugly_hist = set()

for i in range(0,n_of_trial):

	# Generate Histogram object
	Histogram_object = simulation_module.GENERATE_IS_HISTOGRAM (source_distribution, distribution_parameters, span, n_of_events)

	# Add tuple(Histogram_object.frequencies) to a set in order to quantify different kinds of shapes
	all_hist.add(tuple(Histogram_object.frequencies))

	# Extract some Histograms from the whole, take advantage of GOF p-value criterion
	if (Histogram_object.p_value <= p_value_threshold):

		# Counter
		ugly=ugly+1
		# Add tuple(Histogram_object.frequencies) to a set in order to quantify different kinds of shapes
		ugly_hist.add(tuple(Histogram_object.frequencies))

		# Print Info about current Histogram
		print "\n\n\n##### Histogram object {0} #####".format(str(i+1))
		print "Attributes:"
		print "source_distribution = ", Histogram_object.source_distribution
		print "expected_value = ", Histogram_object.expected_value
		print "distribution_parameters = ", Histogram_object.distribution_parameters
		print "n_of_events = ", Histogram_object.n_of_events
		print "bins = ", Histogram_object.bins
		print "discrete_realizations = ", Histogram_object.discrete_realizations
		print "occurrencies = ", Histogram_object.occurrencies
		print "frequencies = ", Histogram_object.frequencies
		print "discrete_realizations_beyond_edges = ", Histogram_object.discrete_realizations_beyond_edges
		print "p_value = ", Histogram_object.p_value
		print "another GOF index = ", Histogram_object.another_GOF_indicator

		#Plot
		Histogram_object.bar_plot(bar_width = 0.8, color = 'blue', title = 'auto', show = show_plot, save = save_plot, name = 'auto', path = 'current_location', id_num = ugly)

		#Amplify
		print "\nNow trying to amplify ... ",
		Amplified_Histogram_object = Histogram_object.amplify()
		print "Done!"

		# Print Info about current Histogram
		print "\nNew Attributes updated:"
		print "amplified = ", Amplified_Histogram_object.amplified
		print "n_of_mapped_sequencies = ", Amplified_Histogram_object.n_of_mapped_sequencies
		print "sequence_count = ", Amplified_Histogram_object.sequence_count
		print "abundances = ", Amplified_Histogram_object.abundances
		print "p_value_post_amplification = ", Amplified_Histogram_object.p_value_post_amplification
		print "another_GOF_indicator_post_amplification = ", Amplified_Histogram_object.another_GOF_indicator_post_amplification

		#Plot again
		Amplified_Histogram_object.bar_plot(bar_width = 0.8, color = 'blue', title = 'auto', show = show_plot, save = save_plot, name = 'auto', path = 'current_location', id_num = ugly)

	# Change distribution_parameters
	#distribution_parameters[source_distribution]['st_dev'] = distribution_parameters[source_distribution]['st_dev'] + 0.001

	# Change n_of_events
	#n_of_events = n_of_events + 100

###################################################################################################################################

### Final print
print "\n\n*** SUMMARY ***"
print "- #N Different Histograms generated = {1}/{0} ({2}%)".format(str(n_of_trial),str(len(all_hist)), str(100*float(len(all_hist))/float(n_of_trial)))
print "- #N Shapes under p={3} threshold level = {0}/{1} ({2}%)".format(str(ugly),str(i+1),str(100*float(ugly)/float(i+1)), str(p_value_threshold))
print "  [among which {0} are different ({1}%)]\n\n".format(str(len(ugly_hist)), str(100*float(len(ugly_hist))/float(ugly)))
print "\n\n[QUIT]\tBye!\n\n"






