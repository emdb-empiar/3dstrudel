"""
mapMotifValidation.py

This module runs map versus amino acid residues map motif library validation.

Copyright [2013] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""

__author__ = 'Andrei Istrate'
__email__ = 'andrei@ebi.ac.uk'
__date__ = '2018-05-29'

import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='Map Validation')
    parser.add_argument("-i", "--input", dest="inp", required=True, help="Strudel map motif validation output directory")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Output directory")
    parser.add_argument("-process_log", "--process_log", dest="process_log", default=None, required=False, help="Log file")
    parser.add_argument("-v", "--voxel", dest="voxel", required=False, default=0.25, type=float, help="Segments voxel size")
    parser.add_argument("-r", "--recompute", dest="recompute_scores", action='store_true',
                        help="Recalculate correlations")
    parser.add_argument("-rs", "--recompute_segments", dest="recompute_segments", action='store_true',
                        help="Repeat segmenting of the input model and map if found in the output")
    parser.add_argument("-vc", "--verbose_chimerax", dest="v_c", action='store_true',
                        help="Run ChimeraX in verbose mode")
    parser.add_argument("-wl", "--warning_level", dest="warning_level", default='info',
                        help="Log file warning level [info[debug]")

    args = parser.parse_args()
    if args.log:
        log_file = args.log
    else:
        log_file = os.path.join(args.out, 'map_motif_validation.process_log')

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(levelname)s:  %(message)s')
    date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
    logging.info('\n{:_^100}'.format('MapMotifValidation') + '\n\nStarted: {}\n'.format(date_time))

    args.np = int(args.np)
    start = time.time()

    compute = ComputeScores()
    compute.set_paths(work_dir=args.out, lib=args.lib, in_map=args.in_map, in_model=args.in_model)
    chopped_pairs = compute.segment_structure(n_cores=args.np, voxel=args.voxel, replace=args.recompute_segments)
    json_file_path = compute.compute_correlations(chopped_pairs, args.np, args.recompute_scores, verbose=args.v_c)

    prefix = os.path.splitext(json_file_path)[0]
    if os.path.exists(prefix + '.csv'):
        csv_to_top_csv(prefix+'.csv', prefix+'_top.csv')
        # csv_to_top_csv_scores_only(prefix + '.csv', prefix + '_top_for_web.csv')

    logging.info('\nElapsed: {}\n{:_^100}'.format(func.report_elapsed(start), ''))


if __name__ == '__main__':
    main()