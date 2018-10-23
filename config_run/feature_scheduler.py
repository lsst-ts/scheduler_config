"""
This is a configuration for a baseline simulation the Feature Based Scheduler. The if statement is used to bypass the
configuration expect when instantiated by the Feature Scheduler Driver. Note that we import and fill in a
SurveyTopology class that informs the (S)OCS about the projects defined on the configuration.

The only things that cannot be changed here are the names of the variables survey_topoly and scheduler. It is possible
as those are expected by the Driver. The way those objects are configured are entirely up to the user though.

09/07/2018 - Ribeiro, T.
"""
import numpy as np
import healpy as hp
import lsst.sims.featureScheduler as fs
from lsst.ts.scheduler.kernel import SurveyTopology

if __name__ == 'config':
    survey_topology = SurveyTopology()
    survey_topology.num_general_props = 1
    survey_topology.general_propos = ["TMA_test"]
    survey_topology.num_seq_props = 0
    survey_topology.sequence_propos = []

    nside = fs.set_default_nside(nside=32)  # Required

    target_maps = {}
    target_maps['r'] = fs.generate_goal_map(nside=nside, NES_fraction=0.,
                                            WFD_fraction=1.0, SCP_fraction=0.,
                                            GP_fraction=0., WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)

    target_2normfactor = {}
    for filtername in target_maps:
        target_2normfactor[filtername] = target_maps[filtername][0]

    norm_factor = fs.calc_norm_factor(target_2normfactor)

    width = (20.,)
    z_pad = (28.,)
    weight = (1.0,)
    height = (80.,)

    # filters = ['u', 'g', 'r', 'i', 'z', 'y']
    filters = ['r']
    surveys = []

    for filtername in filters:
        bfs = list()
        # bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.HourAngle_bonus_basis_function(max_hourangle=4.))
        bfs.append(fs.Target_map_basis_function(filtername=filtername,
                                                target_map=target_maps[filtername][0],
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=norm_factor))
        bfs.append(fs.MeridianStripeBasisFunction(nside=nside, width=width,
                                                  weight=weight,
                                                  height=height,
                                                  zenith_pad=z_pad))
        bfs.append(fs.Aggressive_Slewtime_basis_function(filtername=filtername, nside=nside, order=6., hard_max=120.))
        bfs.append(fs.Avoid_Fast_Revists(filtername=None, gap_min=480., nside=nside))  # Hide region for 0.5 hours

        weights = np.array([1., 0.1, 1., 1., 1.])
        surveys.append(fs.Greedy_survey_fields(bfs, weights, block_size=1,
                                               filtername=filtername, dither=True,
                                               nside=nside,
                                               tag_fields=True,
                                               tag_map=target_maps[filtername][1],
                                               tag_names=target_maps[filtername][2],
                                               ignore_obs='DD'))

    scheduler = fs.Core_scheduler([surveys], nside=nside)  # Required
