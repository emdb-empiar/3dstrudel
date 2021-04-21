import pytest
import numpy as np
from threed_strudel.utils import bio_utils
from threed_strudel.parse.map_parser import MapParser
from threed_strudel.chop.chop_map import ChopMap
from threed_strudel.validation.map_motif_validation import ComputeScores


def test_check_scores_completeness():

    name_list = ['glu-262-A', 'phe-263-A', 'leu-264-A', 'leu-265-A', 'asn-266-A']
    name_list_plus1 = ['glu-262-A', 'phe-263-A', 'leu-264-A', 'leu-265-A', 'asn-266-A', 'asn-267-A']
    scores_path_incomplete = 'data/in/validate/scores_top_incomplete.csv'
    scores_path = 'data/in/validate/scores_top.csv'

    incomplete1 = ComputeScores.check_scores_completeness(name_list_plus1, scores_path, 'some_log')
    incomplete2 = ComputeScores.check_scores_completeness(name_list, scores_path_incomplete, 'some_log')

    complete = ComputeScores.check_scores_completeness(name_list, scores_path, 'some_log')

    assert not incomplete1
    assert not incomplete2
    assert complete
