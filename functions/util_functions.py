import csv, datetime, os
import pickle


def write_lines(csv_file, listoflists, header=False):
    """Write rows to a csv file
    :param csv_file: output csv filename
    :type csv_file: str:
    :param listoflists: data to be written to csv
    :type listoflists: list of lists
    :param header: boolean to enable writing of a single header line to a new csv
    :type header: bool"""

    with open(csv_file, 'a') as f:
        writer = csv.writer(f)
        if header:
            writer.writerow(listoflists)
        else:
            print('Writing results to %s\n' % (csv_file))
            for row in listoflists:
                writer.writerow(row)


def save_pickle(file, destination):
    """Write file to pickle
    :param file: input file
    :param destination: path to savefile
    :type destination: str"""

    print('Saving pickle to ', destination)
    with open(destination, 'wb') as f:
        pickle.dump(file, f, protocol=pickle.HIGHEST_PROTOCOL)
    print('Complete')


def load_pickle(loadfile):
    """Load file from pickle
    :param loadfile: path to file
    :type destination: str
    :return obj: pickled file"""

    print('Loading pickle from', loadfile)
    with open(loadfile, 'rb') as f:
        obj = pickle.load(f)
    print('Complete')
    return obj

def get_time():
    """Generate a timestamp
    :rtype: str"""
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')


def make_resultsfile(path, parameters, pr=False):
    """Generate results file
    :param path: path to csv
    :type path: str
    :param parameters: parameters for readout
    :type parameters: dictionary
    :param pr: flag for performance metrics
    :type pr: bool"""

    if not os.path.exists(path):
        if not pr:
            header=['Experiment',
                    'Model',
                    'N_chains',
                    'retention: actual',
                    'retention: baseline',
                    'retention: dif',
                    'purity: actual',
                    'purity: baseline',
                    'purity: dif',
                    'purity_90: actual',
                    'purity_90: baseline',
                    'purity_90: dif',
                    'consistency: actual',
                    'consistency: baseline',
                    'consistency: dif',
                    'runtime',
                    'Accuracy_train',
                    'Precision_train',
                    'Recall_train',
                    'F1_train',
                    'Precision_weighted_train',
                    'Recall_weighted_train',
                    'F1_weighted_train',
                    'TPR_norm_train',
                    'FPR_norm_train',
                    'Support_train',
                    'Accuracy_test',
                    'Precision_test',
                    'Recall_test',
                    'F1_test',
                    'Precision_weighted_test',
                    'Recall_weighted_test',
                    'F1_weighted_test',
                    'TPR_norm_test',
                    'FPR_norm_test',
                    'Support_test']
                            
            header = header+[x for x in parameters.keys() if x not in ['input_file',
                                                                       'results_file',
                                                                       'chain_selection',
                                                                       'model_selection',
                                                                       'muscle_path',
                                                                       'graphs',
                                                                       'save',
                                                                       'root',
                                                                       'wdir',
                                                                       'name',
                                                                       'experiment']]

        else:
            header = ['Experiment', 'Model',
                      'Train/Test', 'Epitope',
                      'Accuracy', 'Precision',
                      'Recall', 'F1',
                      'Precision_weighted', 'Recall_weighted',
                      'F1_weighted', 'TPR_norm',
                      'FPR_norm', 'Support']

        os.system('touch {}'.format(path))
        write_lines(path, header, header=True)
    