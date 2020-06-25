AA_RESIDUES_3LET_TO_FULL = {'arg': 'arginine',
                            'lys': 'lysine',
                            'met': 'methionine',
                            'glu': 'glutamate',
                            'gln': 'glutamine',
                            'asp': 'aspartate',
                            'asn': 'asparagine',
                            'ile': 'isoleucine',
                            'leu': 'leucine',
                            'his': 'histidine',
                            'trp': 'tryptophan',
                            'tyr': 'tyrosine',
                            'phe': 'phenylalanine',
                            'pro': 'proline',
                            'thr': 'threonine',
                            'val': 'valine',
                            'ser': 'serine',
                            'cys': 'cysteine'
                            }

SYMMETRIC_ATOM_PAIRS_ = {'arg': [['nh1', 'nh2']],
                        'phe': [['ce1', 'ce2'], ['cd1', 'cd2']],
                        'tyr': [['ce1', 'ce2'], ['cd1', 'cd2']],
                        'asp': [['od1', 'od2']],
                        'glu': [['oe1', 'oe2']],
                        }
SYMMETRIC_ATOM_PAIRS = {'arg': ['nh1', 'nh2'],
                        'phe': ['ce1', 'ce2', 'cd1', 'cd2'],
                        'tyr': ['ce1', 'ce2', 'cd1', 'cd2'],
                        'asp': ['od1', 'od2'],
                        'glu': ['oe1', 'oe2'],
                        'his': ['nd1', 'cd2', 'ce1', 'ne2']
                        }

SYMMETRIC_CHI = {'arg': [['cd', 'ne', 'cz', 'nh1'], ['cd', 'ne', 'cz', 'nh2']],
                 'phe': [['ca', 'cb', 'cg', 'cd1'], ['ca', 'cb', 'cg', 'cd2']],
                 'tyr': [['ca', 'cb', 'cg', 'cd1'], ['ca', 'cb', 'cg', 'cd2']],
                 'asp': [['ca', 'cb', 'cg', 'od1'], ['ca', 'cb', 'cg', 'od2']],
                 'glu': [['cb', 'cg', 'cd', 'oe1'], ['cb', 'cg', 'cd', 'oe2']],
                 }


ROTAMER_DATA = {'arg': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd'], 'angle_val': None, 'angle_width': None},
                        'chi_3': {'atoms': ['cb', 'cg', 'cd', 'ne'], 'angle_val': None, 'angle_width': None},
                        'chi_4': {'atoms': ['cg', 'cd', 'ne', 'cz'], 'angle_val': None, 'angle_width': None}
                        },
                'lys': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd'], 'angle_val': None, 'angle_width': None},
                        'chi_3': {'atoms': ['cb', 'cg', 'cd', 'ce'], 'angle_val': None, 'angle_width': None},
                        'chi_4': {'atoms': ['cg', 'cd', 'ce', 'nz'], 'angle_val': None, 'angle_width': None}
                        },
                'met': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'sd'], 'angle_val': None, 'angle_width': None},
                        'chi_3': {'atoms': ['cb', 'cg', 'sd', 'ce'], 'angle_val': None, 'angle_width': None}
                        },
                'glu': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd'], 'angle_val': None, 'angle_width': None},
                        'chi_3': {'atoms': ['cb', 'cg', 'cd', 'oe1', 'oe2'], 'angle_val': None, 'angle_width': None}
                        },
                'gln': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd'], 'angle_val': None, 'angle_width': None},
                        'chi_3': {'atoms': ['cb', 'cg', 'cd', 'oe1'], 'angle_val': None, 'angle_width': None}
                        },
                'asp': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'od1', 'od2'], 'angle_val': None, 'angle_width': None}
                        },
                'asn': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'od1'], 'angle_val': None, 'angle_width': None}
                        },
                'ile': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg1'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg1', 'cd1'], 'angle_val': None, 'angle_width': None}
                        },
                'leu': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd1'], 'angle_val': None, 'angle_width': None}
                        },
                'his': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'nd1', 'cd2'], 'angle_val': None, 'angle_width': None}
                        },
                'trp': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd1'], 'angle_val': None, 'angle_width': None}
                        },
                'tyr': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd1', 'cd2'], 'angle_val': None, 'angle_width': None}
                        },
                'phe': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None},
                        'chi_2': {'atoms': ['ca', 'cb', 'cg', 'cd1', 'cd2'], 'angle_val': None, 'angle_width': None}
                        },
                'pro': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg'], 'angle_val': None, 'angle_width': None}
                        },
                'thr': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'og1'], 'angle_val': None, 'angle_width': None}
                        },
                'val': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'cg1'], 'angle_val': None, 'angle_width': None}
                        },
                'ser': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'og'], 'angle_val': None, 'angle_width': None}
                        },
                'cys': {'chi_1': {'atoms': ['n', 'ca', 'cb', 'sg'], 'angle_val': None, 'angle_width': None}
                        },
                }
