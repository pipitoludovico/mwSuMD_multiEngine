import numpy as np

scores = {1: {'DISTANCE': {'scoreMetric': 16.47467657646962, 'lastValue': 16.55726280241794,
                           'allMetricValues': [16.436928805717045, 16.18331524275113, 16.55726280241794]}}, 2: {
    'DISTANCE': {'scoreMetric': 16.93916461812932, 'lastValue': 17.2303962097033,
                 'allMetricValues': [16.51962191040289, 16.208548322811097, 17.2303962097033]}}}

initialParameters = {'Parameters': ['/scratch/ludovico3/5UHS/equilibration/5UHS_mwSuMD/run3/system/BGC_combined.par',
                                    '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/par_all36_lipid.prm',
                                    '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/par_all36_prot.prm',
                                    '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/par_all36_carb.prm',
                                    '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/par_all36_na.prm',
                                    '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/par_all36_cgenff.prm'],
                     'Forcefield': 'CHARMM', 'Metric_1': 'DISTANCE', 'Metric_2': 'DISTANCE', 'Restart': 'NO',
                     'CUSTOMFILE': None, 'Timestep': 4, 'Savefreq': 20, 'Wrap': 'segid P0', 'Fails': 5,
                     'Tolerance': 0.3, 'NumberCV': 2, 'RelaxTime': 5.0, 'Relax': False, 'CheckEvery': None,
                     'Output': '5UHS_step_1', 'Cutoff_1': 3.0, 'Transition_1': 'negative', 'Cutoff_2': 5.0,
                     'Transition_2': 'negative', 'Walkers': 2, 'Timewindow': 20, 'coor': 'previous.coor',
                     'xsc': 'previous.xsc', 'vel': 'previous.vel'}

par = {'Transition_1': 'negative', 'Transition_2': 'negative'}
walkers_metrics_1 = [16.55726280241794]
walkers_metrics_2 = [17.2303962097033]

allMetric_1 = [16.436928805717045, 16.18331524275113, 16.55726280241794]
allMetric_2 = [16.51962191040289, 16.208548322811097, 17.2303962097033]


def bestWalker_2metrics(par, walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2):
    print('INTO THE FUNCTION: orignal 1 and 2 + all metrics 1 and 2', walkers_metrics_1, walkers_metrics_2, allMetric_1,
          allMetric_2)
    scores_walkers_metric_1 = []
    scores_walkers_metric_2 = []

    metric_1_avg = sum(allMetric_1) / len(allMetric_1)
    metric_2_avg = sum(allMetric_2) / len(allMetric_2)
    print('AVERAGE M1', metric_1_avg)
    print('AVERAGE M2', metric_2_avg)

    # final score computation
    for i in walkers_metrics_1:
        score_1 = (i - metric_1_avg) * (100 / metric_1_avg)
        print("SCORE 1, ", score_1, " AVERAGE", metric_1_avg)

        if par['Transition_1'] == 'negative':
            score_1 = score_1 * -1
        scores_walkers_metric_1.append(score_1)
    for i in walkers_metrics_2:
        score_2 = (i - metric_2_avg) * (100 / metric_2_avg)

        if par['Transition_2'] == 'negative':
            score_2 = score_2 * -1
        print("original value - score metric 2", i, score_2)
        scores_walkers_metric_2.append(score_2)

    # put the final scores for the two metrics together for each walker and sum them up
    print("TWO scores LISTS", scores_walkers_metric_1, scores_walkers_metric_2)
    scores_list = []
    for (score1, score2) in zip(scores_walkers_metric_1, scores_walkers_metric_2):
        if score1 != None and score2 != None:
            scores_list.append(score1 + score2)
        elif score1 == None or score2 == None:
            scores_list.append(None)
    # get best 	walker accoriding to the final scores sum
    max_value = max([i for i in scores_list if i is not None])
    # max_value = max(scores_list)
    max_index = scores_list.index(max_value) + 1
    print("SCORES FINAL", scores_list, max_index)


def getBestWalker(scores: dict):
    """Returns the walker with the best metric"""
    try:
        if initialParameters['NumberCV'] == 1:

            if initialParameters['Transition_1'] == 'positive':
                best_metric_key = max(scores, key=lambda k: scores[k][next(iter(scores[k]))]['allMetricValues'])
                best_metric_value = scores[best_metric_key][next(iter(scores[best_metric_key]))]['scoreMetric']
                metric_used = next(iter(scores[best_metric_key]))
            else:
                best_metric_key = min(scores, key=lambda k: scores[k][next(iter(scores[k]))]['allMetricValues'])
                best_metric_value = scores[best_metric_key][next(iter(scores[best_metric_key]))]['scoreMetric']
                metric_used = next(iter(scores[best_metric_key]))
            return best_metric_key, best_metric_value, scores[best_metric_key][metric_used]['allMetricValues'][-1]
        else:
            averages = {}
            last_values = {}
            score_ = {}
            keyMax = 0

            for walker, result in scores.items():
                averages[walker] = {}
                last_values[walker] = {}
                for metric, data in result.items():
                    all_metric_values = data['allMetricValues']
                    last_values[walker] = data['lastValue']
                    average = np.average(all_metric_values)
                    averages[walker][metric] = average
                    if initialParameters[f'Transition_{walker}'] == 'negative':
                        k = -1
                    else:
                        k = 1

                    score_[walker] = k * (last_values[walker] - average) * (100 / average)

            return max(score_, key=score_.get), sum(score_.values()), last_values
    except FileExistsError:
        print("no")


bestWalker_2metrics(par, walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2)
print("\n\n")
getBestWalker(scores)
