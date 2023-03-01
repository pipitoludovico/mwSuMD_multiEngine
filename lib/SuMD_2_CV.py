from lib.FolderOps import *
from lib.metrics import *


####################
def metricCompute(par, selection_list, n, mothfolder):
	
	os.chdir('tmp')
	walkers_metrics_1 = []
	walkers_metrics_2 = []
	allMetric_1 = []
	allMetric_2 = []
	for r in range(1, int(par['Walkers'])+1):
		os.chdir('walker_%s' %str(r))
		if (glob.glob('*.xtc') and glob.glob('output.*')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
# First metric is a distance 
			if par['Metric_1'] == 'Distance':
# get all the distances and the last one
				distances, last_distance = distance(par, selection_list[0], selection_list[1], n, mothfolder)
# add the last distance of the walket to the list
				walkers_metrics_1.append(last_distance)
# keep all the distances from all the walkers together
				allMetric_1 = allMetric_1 + distances

			
			if par['Metric_1'] == 'Contacts':
				timeseries, last_timeseries = contacts(par, selection_list[0], selection_list[1], n, mothfolder)
				walkers_metrics_1.append(last_timeseries)
				allMetric_1 = allMetric_1 + timeseries	

			
			if par['Metric_1'] == 'RMSD':
				data,last_rmsd = RMSD(par, selection_list[0], selection_list[1], n, mothfolder)
				walkers_metrics_1.append(last_rmsd)	
				allMetric_1 = allMetric_1 + data		
			

			if par['Metric_2'] == 'Distance':
				distances, last_distance = distance(par, selection_list[2], selection_list[3], n, mothfolder)
				walkers_metrics_2.append(last_distance)
				allMetric_2 = allMetric_2 + distances
			
			if par['Metric_2'] == 'Contacts':
				timeseries, last_timeseries = contacts(par, selection_list[2], selection_list[3], n, mothfolder)
				walkers_metrics_2.append(last_timeseries)
				allMetric_2 = allMetric_2 + timeseries

			
			if par['Metric_2'] == 'RMSD':
				data,last_rmsd = RMSD(par, selection_list[2], selection_list[3], n, mothfolder)
				walkers_metrics_2.append(last_rmsd)
				allMetric_2 = allMetric_2 + data
			
			if par['Metric_2'] == 'HB_score':
				scores,last_score = HB_score(par, selection_list, n, mothfolder)
				print("HBScore", scores, last_score)
				walkers_metrics_2.append(last_score)
				allMetric_2 = allMetric_2 + scores	
		else:
			os.chdir('..')
			walkers_metrics_1.append(None)
			walkers_metrics_2.append(None)

		os.chdir(mothfolder+'/tmp')		

	os.chdir(mothfolder)
	return walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2


############################################# Get list of matric values and computr final score using average value from all walkers
def bestWalker_2metrics(par, walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2):
	
	if all(v is None for v in walkers_metrics_1) or all(v is None for v in walkers_metrics_2):
		return (0,0)
	
	print('INTO THE FUNCTION: orignal 1 and 2 + all metrics 1 and 2', walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2)
	scores_walkers_metric_1 = []
	scores_walkers_metric_2 = []
	
	metric_1_avg = sum(allMetric_1)/len(allMetric_1)
	metric_2_avg = sum(allMetric_2)/len(allMetric_2)
	print('AVERAGE M1', metric_1_avg)
	print('AVERAGE M2', metric_2_avg)
	
# final score computation
	for i in walkers_metrics_1:
		if i != None:
			score_1 = (i - metric_1_avg)*(100/metric_1_avg)

			if par['Transition_1'] == 'negative':
				score_1 = score_1 * -1
			print("original value - score metric 1", i,score_1)
			scores_walkers_metric_1.append(score_1)
		elif i == None:
			score_1 = None
	for i in walkers_metrics_2:
		if i != None:
			score_2 = (i - metric_2_avg)*(100/metric_2_avg)
		
			if par['Transition_2'] == 'negative':
				score_2 = score_2 * -1
			print("original value - score metric 2",i, score_2)	
			scores_walkers_metric_2.append(score_2)	
		elif i == None:
			score_2 = None
		#if par['Transition_2'] == None:
			#score_2 = score_2 * -1
		#print("original value - score metric 2",i, score_2)	
		#scores_walkers_metric_2.append(score_2)			

# put the final scores for the two metrics together for each walker and sum them up
	print("TWO scores LISTS",scores_walkers_metric_1, scores_walkers_metric_2)
	scores_list = []
	for (score1, score2) in zip(scores_walkers_metric_1, scores_walkers_metric_2):
		if score1 != None and score2 != None:
			scores_list.append(score1+score2)
		elif score1 == None or score2 == None:
			scores_list.append(None)
	# get best 	walker accoriding to the final scores sum
	max_value = max([i for i in scores_list if i is not None])
	#max_value = max(scores_list)
	max_index = scores_list.index(max_value)+1
	print("SCORES FINAL", scores_list, max_index)	
	return scores_list, max_index

###############################################################
###############################################################
def SuMD_2_CV(par, selection_list, PARPATH):
	
	
	mothfolder = os.getcwd()

	if par['Restart'] == 'YES':
		n = restart(mothfolder)
		c = 0
	

	if par['Restart'] == 'NO':
		n = 0
		c = 0
		make_folders()

	max_cycles = 1/(int(par['Timewindow'])/10**5) # run for 1 microsecond and then stop

	while c < max_cycles:
		# run acemd for each walker
		runMD(par, n, mothfolder, PARPATH)
		# compute two metrics and get also all the values from all the walkers 
		walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2 = metricCompute(par, selection_list, n, mothfolder)
		# compute scores and decide best walker
		scores_list, best_score = bestWalker_2metrics(par, walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2)
		if scores_list == 0:
			continue
		# write log
		logStep(par, n, scores_list, best_score, mothfolder, walkers_metrics_1, walkers_metrics_2)
		# save trj and restart files
		saveStep(par, best_score, n, mothfolder)
		n += 1
		c +=0
