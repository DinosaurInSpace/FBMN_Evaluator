#!/usr/bin/env python

"""Feature-Based Molecular Networking (FBMN) is a powerful unsupervised approach for organizing
LCMS data bases on retention time, MS1, and MS2 over multiple samples.

https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

Louis Felix Nothias, et. al.
bioRxiv 812404; doi: https://doi.org/10.1101/812404

A key challenge in generating FBMN's is determining appropriate processing settings when MS1
features are generated, using MZmine or other tools.  If merge settings are too permissive,
different molecular species will be inappropriately merged, whereas if merge settings are
too restrictive, the same molecular species will be inappropriately split.  The first case
artificially hides molecular diversity in the sample, whereas the second over estimates it.
The second case is especially a challenge for downstream visualization (e.g. 'ili).

Here a simple command line Python script is presented for comparing the relative quality of
seveal FBMN's is presented.  FBMN's can be evaluated as a visual plot of retention timeand
m/z features assigned via MS/MS.  As well, figures such as the total number of FBMN features,
and FBMN features within a certain retention time and m/z tolerance are reported.

Usage:

1. Run jobs via FMBN workflow.
2. After job is complete, click on "done" then "Download Cytoscape Data"
3. Unzip the network, and in the "quantification_table" directory there is your network!
4. Save one or multiple networks to a known directory
5. Go to appropriate command line interface: terminal (Mac/Linux) or Windows Power Shell
6. Run command as below, substituting uppercase text:

python FBMN_Evaluator.py --path PATH_TO_NETWORK_DIR --mztol MZ_TOLERANCE_IN_PPM --rttol
RETENTION_TIME_TOL_IN_SEC

Output:

Path to plots, number of FBMN features, number of FBMN features within mz/rt tolerance.

"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import time
from datetime import datetime
import argparse
import pathlib
import re
import textwrap
from os import path, listdir
from scipy.stats import gaussian_kde

__author__ = "Christopher M Baxter Rath"
__copyright__ = "Copyright 2020"
__credits__ = ["Christopher M Baxter Rath"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Christopher M Baxter Rath"
__email__ = "chrisrath@gmail.com"
__status__ = "Development"


def input_parser(input_file, input_dir):
    df = pd.read_csv(input_dir + '/' + input_file)
    col_list = ['row ID', 'row m/z', 'row retention time']
    col_dict = {'row ID': 'id', 'row m/z': 'mz', 'row retention time': 'rt'}
    df = df[col_list].copy(deep=True)
    df = df.rename(columns=col_dict)
    return df


def ppm_to_mz(mz, ppm):
    d_mz = (ppm / 1000000 ) * mz
    return d_mz


def tol_features(fbmn_df, mz_tol, rt_tol):
    # Cytoscape output is in minutes!
    rt_tol = rt_tol / 60
    features_at_tol = 0

    # Range is 0 to number of features
    for row in fbmn_df.itertuples():
        r_mz = ppm_to_mz(row.mz, mz_tol)
        t_df = fbmn_df.copy(deep=True)
        result = t_df[((t_df.mz <= row.mz + r_mz) & (t_df.mz >= row.mz - r_mz) &
                      ((t_df.rt <= row.rt + rt_tol) & (t_df.rt >= row.rt - rt_tol)))]

        # Every point will find themselves!
        if result.shape[0] == 1:
            pass
        else:
            features_at_tol +=1
    return features_at_tol


def mz_rt_plotter(fbmn_df, mz_tol, rt_tol, tot_feat, tol_feat, input_file):
    #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib

    x = fbmn_df.rt
    y = fbmn_df.mz

    title_dict = {'Input file': input_file,
                  'Tested mz tolerance (ppm)': mz_tol,
                  'Tested rt tol (sec.)': rt_tol,
                  'Total features': tot_feat,
                  'Tolerance features': tol_feat}

    filename = str(datetime.now())
    filename = re.sub('[^0-9a-zA-Z]+', '_', filename)
    filename = 'reports/' + filename
    print(filename)

    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=10, edgecolor='')
    ax.set_xlabel('Retention time in minutes')
    ax.set_ylabel('m/z values')
    ax.set_title("\n".join(textwrap.wrap(str(title_dict), 60)))
    plt.tight_layout()

    plt.savefig(filename)
    plt.show()
    plt.close()

    return filename


def filename_path_timestamp():
    path = pathlib.Path().absolute()
    filename = str(datetime.now())
    filename = re.sub('[^0-9a-zA-Z]+', '_', filename)
    filename = str(path) + '/reports/' +str(filename)

    return filename


def fbmn_main(input_dir, mz_tol, rt_tol):
    start_time = time.time()
    output_dict = {}

    if path.isdir(input_dir) and listdir(path=input_dir) != 0:
        for input_file in listdir(input_dir):
            fbmn_df = input_parser(input_file, input_dir)
            tot_feat = fbmn_df.shape[0]
            tol_feat = tol_features(fbmn_df, mz_tol, rt_tol)
            plot_path = mz_rt_plotter(fbmn_df, mz_tol, rt_tol, tot_feat, tol_feat, input_file)

            output_dict[input_file] = {'mz_tol_ppm': mz_tol,
                                       'rt_tol_sec': rt_tol,
                                       'total_feat': tot_feat,
                                       'tolerance_feat': tol_feat,
                                       'plot_path': plot_path}

        output_df = pd.DataFrame(output_dict)
        path_to_df = filename_path_timestamp()
        output_df.to_pickle(path_to_df + '.pickle')

        elapsed_time = time.time() - start_time
        print('Elapsed time:\n', elapsed_time, '\nExecuted without error\n')
        print(output_df)
        print(path_to_df)

        return (output_df, path_to_df)

    else:
        print('Incorrect input directory or empty directory!')
        exit(1)


### Main ###
parser = argparse.ArgumentParser(description='')
parser.add_argument("--path", default='networks/', type=str, help="Directory with GNPS output")
parser.add_argument("--mztol", default=20, type=int, help="Mass error in ppm")
parser.add_argument("--rttol", default=20, type=int, help="RT tolerance in seconds")
args = parser.parse_args()

# Main loop
fbmn_main(args.path, args.mztol, args.rttol)