###Import Module(s)#####
import simulation_module
########################

### Set-up variables ##############################
source_distribution = 'gauss'
distribution_parameters = {'gauss':{'st_dev':1.0}}
span = 7
n_of_events = 20
n_of_trial = 100
###################################################

print "\n\n[START]\tTest Environment!\n\n"

### Loop ##########################################################################################################################
ugly = 0
all_hist = set()
ugly_hist = set()
for i in range(0,n_of_trial):

	# Generate Histogram object
	Histogram_object = simulation_module.GENERATE_IS_HISTOGRAM (source_distribution, distribution_parameters, span, n_of_events)
	all_hist.add(tuple(Histogram_object.frequencies))

	if (Histogram_object.p_value <= 0.05):

		ugly=ugly+1
		ugly_hist.add(tuple(Histogram_object.frequencies))
		print "\nHistogram object {0} ... ".format(str(i+1)),
		print "Attributes:"
		print "source_distribution = ", Histogram_object.source_distribution
		print "expected_value = ", Histogram_object.expected_value
		print "distribution_parameters = ", Histogram_object.distribution_parameters
		print "n_of_events = ", Histogram_object.n_of_events
		print "bins = ", Histogram_object.bins
		print "frequencies = ", Histogram_object.frequencies
		print "realization_beyond_edges = ", Histogram_object.realization_beyond_edges
		print "p_value = ", Histogram_object.p_value
		print "another GOF index = ", Histogram_object.another_GOF_indicator
		print "\n***\t***\t***\t***\n"

		#Plot
		Histogram_object.bar_plot(bar_width = 0.8, color = 'blue', title = 'auto', show = False, save = True, name = 'auto', path = 'current_location', id_num = ugly)

		#Amplify
		Amplified_Histogram_object = Histogram_object.amplify()

		#Plot again
		Amplified_Histogram_object.bar_plot(bar_width = 0.8, color = 'blue', title = 'auto', show = False, save = True, name = 'auto', path = 'current_location', id_num = ugly)

	# Change distribution_parameters
	#distribution_parameters[source_distribution]['st_dev'] = distribution_parameters[source_distribution]['st_dev'] + 0.001

	# Change n_of_events
	#n_of_events = n_of_events + 100
###################################################################################################################################
print "\n\n*** Different Histograms over {0} generated = {1} ***".format(str(n_of_trial),str(len(all_hist)))
print "*** #N Shapes under 0.05 threshold level = {0}/{1} = {2} ***".format(str(ugly),str(i+1),str(float(ugly)/float(i+1)))
print "*** [among which {0} are different] ***\n\n".format(str(len(ugly_hist)))
print "\n\n[QUIT]\tBye!\n\n"






