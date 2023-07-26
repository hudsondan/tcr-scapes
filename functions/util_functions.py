import csv, datetime
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
    