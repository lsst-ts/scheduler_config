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
    survey_topology.num_general_props = 4
    survey_topology.general_propos = ["NorthEclipticSpur", "SouthCelestialPole", "WideFastDeep", "GalacticPlane"]
    survey_topology.num_seq_props = 1
    survey_topology.sequence_propos = ["DeepDrillingCosmology1"]

    nside = fs.set_default_nside(nside=32)

    target_maps = {}
    target_maps['u'] = fs.generate_goal_map(nside=nside, NES_fraction=0.,
                                            WFD_fraction=0.31, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)
    target_maps['g'] = fs.generate_goal_map(nside=nside, NES_fraction=0.2,
                                            WFD_fraction=0.44, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)
    target_maps['r'] = fs.generate_goal_map(nside=nside, NES_fraction=0.46,
                                            WFD_fraction=1.0, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)
    target_maps['i'] = fs.generate_goal_map(nside=nside, NES_fraction=0.46,
                                            WFD_fraction=1.0, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)
    target_maps['z'] = fs.generate_goal_map(nside=nside, NES_fraction=0.4,
                                            WFD_fraction=0.9, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)
    target_maps['y'] = fs.generate_goal_map(nside=nside, NES_fraction=0.,
                                            WFD_fraction=0.9, SCP_fraction=0.15,
                                            GP_fraction=0.15, WFD_upper_edge_fraction=0.,
                                            generate_id_map=True)

    even_year_target = {}
    odd_year_target = {}
    for fname in target_maps:
        even_year_target[fname] = target_maps[fname][0]
        odd_year_target[fname] = target_maps[fname][0]

    up = 1.75
    down = 0.25

    # Let's find the healpix that divides the WFD area in half
    wfd = even_year_target['r'] * 0
    wfd[np.where(even_year_target['r'] == 1)] = 1
    wfd_accum = np.cumsum(wfd)

    split_indx = np.max(np.where(wfd_accum < wfd_accum.max() / 2.))

    indx = np.arange(even_year_target['r'].size)
    top_half_wfd = np.where((even_year_target['r'] == 1) & (indx <= split_indx))
    bottom_half_wfd = np.where((even_year_target['r'] == 1) & (indx > split_indx))

    for filtername in even_year_target:
        even_year_target[filtername][top_half_wfd] *= up
        even_year_target[filtername][bottom_half_wfd] *= down

        odd_year_target[filtername][top_half_wfd] *= down
        odd_year_target[filtername][bottom_half_wfd] *= up

    even_norm = fs.calc_norm_factor(even_year_target)
    odd_norm = fs.calc_norm_factor(odd_year_target)

    surveys = []
    mod_year = 2
    offset = 1
    # Set up observations to be taken in blocks
    filter1s = ['u', 'g', 'r', 'i', 'z', 'y']
    filter2s = [None, 'g', 'r', 'i', None, None]
    for filtername, filtername2 in zip(filter1s, filter2s):
        bfs = []
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        if filtername2 is not None:
            bfs.append(fs.M5_diff_basis_function(filtername=filtername2, nside=nside))
        bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername,
                                                    target_map=even_year_target[filtername],
                                                    mod_year=mod_year, offset=0,
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=even_norm))
        if filtername2 is not None:
            bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername2,
                                                        target_map=even_year_target[filtername2],
                                                        mod_year=mod_year, offset=0,
                                                        out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                        norm_factor=even_norm))

        bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername,
                                                    target_map=odd_year_target[filtername],
                                                    mod_year=mod_year, offset=offset,
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=odd_norm))
        if filtername2 is not None:
            bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername2,
                                                        target_map=odd_year_target[filtername2],
                                                        mod_year=mod_year, offset=offset,
                                                        out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                        norm_factor=odd_norm))
        bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
        bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
        weights = np.array([3.0, 3.0, .3, .3, 0.3, 0.3, 3., 3., 0.])
        if filtername2 is None:
            # Need to scale weights up so filter balancing still works properly.
            weights = np.array([6.0, 0.6, 0.6, 3., 3., 0.])
        # XXX-
        # This is where we could add a look-ahead basis function to include m5_diff in the future.
        # Actually, having a near-future m5 would also help prevent switching to u or g right at twilight?
        # Maybe just need a "filter future" basis function?
        if filtername2 is None:
            survey_name = 'blob, %s' % filtername
        else:
            survey_name = 'blob, %s%s' % (filtername, filtername2)
        surveys.append(fs.Blob_survey(bfs, weights, filtername=filtername, filter2=filtername2,
                                      survey_note=survey_name,
                                     tag_fields=True,
                                     tag_map=target_maps[filtername][1],
                                     tag_names=target_maps[filtername][2]))

    # Set up the greedy surveys for filling time when can't take pairs.
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    greedy_surveys = []
    for filtername in filters:
        bfs = []
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername,
                                                    mod_year=mod_year, offset=0,
                                                    target_map=even_year_target[filtername],
                                                    norm_factor=even_norm))

        bfs.append(fs.Target_map_modulo_basis_function(filtername=filtername,
                                                    mod_year=mod_year, offset=offset,
                                                    target_map=even_year_target[filtername],
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=odd_norm))

        bfs.append(fs.North_south_patch_basis_function(zenith_min_alt=50., nside=nside))
        bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
        bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
        weights = np.array([3.0, 0.3, 0.3, 1., 3., 3., 0.])
        # Might want to try ignoring DD observations here, so the DD area gets covered normally--DONE
        sv = fs.Greedy_survey_fields(bfs, weights, block_size=1, filtername=filtername,
                                     dither=True, nside=nside, ignore_obs='DD',
                                     tag_fields=True,
                                     tag_map=target_maps[filtername][1],
                                     tag_names=target_maps[filtername][2])
        greedy_surveys.append(sv)

    # Set up the DD surveys
    dd_surveys = fs.generate_dd_surveys()

    survey_list_o_lists = [dd_surveys, surveys, greedy_surveys]

    scheduler = fs.Core_scheduler(survey_list_o_lists, nside=nside)

