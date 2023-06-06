from sampling.seq_gen_trimodal import *
import matplotlib.pyplot as plt
import argparse


def analyse_and_sample(seq_gen_temp, tolA, tolB):
    # sample sequence
    seq = seq_gen_temp.sample(seqlen)[:, 1]
    seq += 1

    dev_mat = np.zeros((seqlen, 6))
    cong_mat = np.zeros((seqlen, 4))
    incong_mat = np.zeros((seqlen, 4))

    taccount = 1
    audcount = 1
    viscount = 1

    for e in range(1, seqlen):
        # ----------------------------------------------------------------------
        ev = seq[e]
        ev1 = seq[e-1]

        if ev == 1:
            # ------------------------------------------------------------------
            if ev1 == 2:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                cong_mat[e, 0] = 1
            elif ev1 == 3:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                cong_mat[e, 1] = 1
            elif ev1 == 4:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                cong_mat[e, 2] = 1
            elif ev1 == 1:                                 # Global Rep
                cong_mat[e, 3] = 1

        elif ev == 2:
            # ------------------------------------------------------------------
            if ev1 == 1:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                cong_mat[e, 0] = 1
            elif ev1 == 5:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                incong_mat[e, 1] = 1
            elif ev1 == 6:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                incong_mat[e, 2] = 1
            elif ev1 == 2:                                 # Global Rep
                cong_mat[e, 3] = 1

        elif ev == 3:
            # ------------------------------------------------------------------
            if ev1 == 5:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                incong_mat[e, 0] = 1
            elif ev1 == 1:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                cong_mat[e, 1] = 1
            elif ev1 == 7:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                incong_mat[e, 2] = 1
            elif ev1 == 3:                                 # Global Rep
                cong_mat[e, 3] = 1

        elif ev == 4:
            # ------------------------------------------------------------------
            if ev1 == 6:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                incong_mat[e, 0] = 1
            elif ev1 == 7:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                incong_mat[e, 1] = 1
            elif ev1 == 1:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                cong_mat[e, 2] = 1
            elif ev1 == 4:                                 # Global Rep
                cong_mat[e, 3] = 1

        elif ev == 5:
            # ------------------------------------------------------------------
            if ev1 == 3:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                incong_mat[e, 0] = 1
            elif ev1 == 2:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                incong_mat[e, 1] = 1
            elif ev1 == 8:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                cong_mat[e, 2] = 1
            elif ev1 == 5:                                 # Global Rep
                incong_mat[e, 3] = 1

        elif ev == 6:
            # ------------------------------------------------------------------
            if ev1 == 4:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                incong_mat[e, 0] = 1
            elif ev1 == 8:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                cong_mat[e, 1] = 1
            elif ev1 == 2:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                incong_mat[e, 2] = 1
            elif ev1 == 6:                                 # Global Rep
                incong_mat[e, 3] = 1

        elif ev == 7:
            # ------------------------------------------------------------------
            if ev1 == 8:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                cong_mat[e, 0] = 1
            elif ev1 == 4:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                incong_mat[e, 1] = 1
            elif ev1 == 3:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                incong_mat[e, 2] = 1
            elif ev1 == 7:                                 # Global Rep
                incong_mat[e, 3] = 1

        elif ev == 8:
            # ------------------------------------------------------------------
            if ev1 == 7:                                     # Tactile Deviant
                dev_mat[e, 0] = 1
                dev_mat[e, 1] = taccount
                taccount = 0
                audcount += 1
                viscount += 1
                cong_mat[e, 0] = 1
            elif ev1 == 6:                                   # Auditory Deviant
                dev_mat[e, 2] = 1
                dev_mat[e, 3] = audcount
                audcount = 0
                taccount += 1
                viscount += 1
                cong_mat[e, 1] = 1
            elif ev1 == 5:                                  # Visual Deviant
                dev_mat[e, 4] = 1
                dev_mat[e, 5] = viscount
                viscount = 0
                taccount += 1
                audcount += 1
                cong_mat[e, 2] = 1
            elif ev1 == 8:                                 # Global Rep
                incong_mat[e, 3] = 1


    # Change Probability: Congruent / Incongruent
    # -------------------------------------------------------------------------------

    # tbu probabilities
    p1, p2, p3, p4 = prob_obs_change[0], prob_obs_change[1], prob_obs_change[2], prob_obs_change[3]
    # tbu: global repetition
    grA = 2*p1/8
    grB = 6*p4/8
    # tbu: modality change
    mcA = 4*p2/8
    mcB = 4*p3/8
    mc = (4*p2+4*p3)/8
    # empirical: congruent
    t_cong = np.round(np.sum(cong_mat[:, 0])/seqlen, 3)
    a_cong = np.round(np.sum(cong_mat[:, 1])/seqlen, 3)
    v_cong = np.round(np.sum(cong_mat[:, 2])/seqlen, 3)
    g_cong = np.round(np.sum(cong_mat[:, 3])/seqlen, 3)
    # empirical: incongruent
    t_incong = np.round(np.sum(incong_mat[:, 0])/seqlen, 3)
    a_incong = np.round(np.sum(incong_mat[:, 1])/seqlen, 3)
    v_incong = np.round(np.sum(incong_mat[:, 2])/seqlen, 3)
    g_incong = np.round(np.sum(incong_mat[:, 3])/seqlen, 3)

    # Stimulus Probability: High / Low
    # -------------------------------------------------------------------------------

    taccode = [1, 2, 1, 1, 2, 2, 1, 2]
    audcode = [1, 1, 2, 1, 2, 1, 2, 2]
    viscode = [1, 1, 1, 2, 1, 2, 2, 2]
    tac = np.empty((seq.shape))
    tac[:] = np.NaN
    aud = np.empty((seq.shape))
    aud[:] = np.NaN
    vis = np.empty((seq.shape))
    aud[:] = np.NaN
    for i in range(len(taccode)):
        tac[seq == i+1] = taccode[i]
        aud[seq == i+1] = audcode[i]
        vis[seq == i+1] = viscode[i]

    p_tac_low = np.round(np.count_nonzero(tac == 1)/seqlen, 3)
    p_tac_hi = np.round(np.count_nonzero(tac == 2)/seqlen, 3)

    p_aud_low = np.round(np.count_nonzero(aud == 1)/seqlen, 3)
    p_aud_hi = np.round(np.count_nonzero(aud == 2)/seqlen, 3)

    p_vis_low = np.round(np.count_nonzero(vis == 1)/seqlen, 3)
    p_vis_hi = np.round(np.count_nonzero(vis == 2)/seqlen, 3)

    # Test for conditions
    # -------------------------------------------------------------------------------

    redo_flag = 0

    # 1) Check if cong/incong are roughly equal
    if p2 == p3:
        tol = 0.005
        maxTchange = max(t_cong, t_incong)
        minTchange = min(t_cong, t_incong)
        maxAchange = max(a_cong, a_incong)
        minAchange = min(a_cong, a_incong)
        maxVchange = max(v_cong, v_incong)
        minVchange = min(v_cong, v_incong)
        changediff = np.append([maxTchange-minTchange], [maxAchange-minAchange, maxVchange-minVchange])
        if any(changediff > tol):
            redo_flag = 1

    # 2) Check that general change prob is not too low
    tol = 0.025
    changep = np.append([t_cong+t_incong], [a_cong+a_incong, v_cong+v_incong])
    if any(changep < [mcA+mcB-tol]):
        redo_flag = 1

    # 3) Check that Congruent and Incongruent are not different across modalities
    tol = 0.01
    maxCong = max(t_cong, a_cong, v_cong)
    minCong = min(t_cong, a_cong, v_cong)
    if maxCong-minCong > tol:
        redo_flag = 1

    tol = 0.01
    maxIncong = max(t_incong, a_incong, v_incong)
    minIncong = min(t_incong, a_incong, v_incong)
    if maxIncong-minIncong > tol:
        redo_flag = 1

    #bp()

    # 2) Check that the improbable change prob (congruent or incongruent) is not too different from tbu
    tol = 0.005
    if mcA < mcB:
        # congruent is improbable
        if not mcA-tol <= t_cong <= mcA+tol or not mcA-tol <= a_cong <= mcA+tol or not mcA-tol <= v_cong <= mcA+tol:
            redo_flag = 1
    elif mcA > mcB:
        # incongruent is improbable
        if not mcB-tol <= t_incong <= mcB+tol or not mcB-tol <= a_incong <= mcB+tol or not mcB-tol <= v_incong <= mcB+tol:
            redo_flag = 1

    # 3) Check that low / high probabilities are not too different from tbu
    tol = 0.025
    if not 0.5-tol <= p_tac_low <= 0.5+tol or not 0.5-tol <= p_tac_hi <= 0.5+tol:
        redo_flag = 1
    if not 0.5-tol <= p_aud_low <= 0.5+tol or not 0.5-tol <= p_aud_hi <= 0.5+tol:
        redo_flag = 1
    if not 0.5-tol <= p_vis_low <= 0.5+tol or not 0.5-tol <= p_vis_hi <= 0.5+tol:
        redo_flag = 1

    cong = [t_cong, a_cong, v_cong]
    incong = [t_incong, a_incong, v_incong]
    low = [p_tac_low, p_aud_low, p_vis_low]
    high = [p_tac_hi, p_aud_hi, p_vis_hi]
    analysis = [cong, incong, low, high]
    tbu = [mcA, mcB, 0.5, 0.5]

    return redo_flag, seq, analysis, tbu


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-sj', '--subject', action="store",
                        default=999, type=int,
                        help="subject number")
    parser.add_argument('-p', '--probability', action="store",
                        default=[.3, .1, .05, .6], nargs='+', type=float,
                        help="probability of deviant dependence [dev|congruent, dev|incongruent]")
    parser.add_argument('-v', '--verbose', action="store",
                        default=1, type=int,
                        help="verbosity")
    parser.add_argument('-id', '--seqid', action="store",
                        default='A', type=str,
                        help="Sequence ID: A or B")
    args = parser.parse_args()

    # prepare_sequence:
    seqlen = 600
    order = 1
    prob_catch = 0
    prob_regime_init = 0
    prob_regime_change = 0
    prob_obs_init = 0
    prob_obs_change = args.probability
    verbose = args.verbose
    results_dir = 'E:/TESTING/stimuli/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    title = 'sub-{num:02d}_seq{id}_{a}_{b}_{c}_{d}'.format(num=args.subject, id=args.seqid,
                                          a=str(args.probability[0])[2:], b=str(args.probability[1])[2:],
                                          c=str(args.probability[2])[2:], d=str(args.probability[3])[2:])
    seq_gen_temp = seq_gen(prob_obs_change, verbose)

    tolA = 0.005
    tolB = 0.05

    redo_flag = 1
    whilecount = 0
    while redo_flag:
        redo_flag, seq, analysis, tbu = analyse_and_sample(seq_gen_temp, tolA, tolB)
        whilecount += 1
        # print(whilecount)
        if whilecount >= 50000:
            if verbose:
                print('Reached max. iterations. No appropriate sequence generated.')
            break

    if verbose:
        print(whilecount)

    sequence_meta = {"sample_output": seq-1,
                     "prob_obs_change": seq_gen_temp.prob_obs_change}
    savemat(results_dir + title, sequence_meta)

    # plots
    if verbose:
        fig, ax = plt.subplots(nrows=2, figsize=(4, 3), ncols=1, dpi=300)
        fig.tight_layout()

        labels = ['Tactile', 'Auditory', 'Visual']
        x = np.arange(len(labels))
        width = .2
        ax[0].bar(x - width/2, analysis[0], width, label='Congruent')
        ax[0].bar(x + width/2, analysis[1], width, label='Incongruent')
        ax[0].plot(x, len(labels) * [tbu[0]], 'r--', x, len(labels) * [tbu[1]], 'r--')

        ax[1].bar(x - width/2, analysis[2], width, label='Low')
        ax[1].bar(x + width/2, analysis[3], width, label='High')
        ax[1].plot(x, len(labels) * [tbu[2]], 'r--', x, len(labels) * [tbu[3]], 'r--')

        ax[0].set_title('Sequence probabilities')
        ax[0].set_ylabel('Probability')
        ax[1].set_ylabel('Probability')
        ax[0].set_xticks(x)
        ax[1].set_xticks(x)
        ax[1].set_xticklabels(labels)
        ax[0].legend()
        ax[1].legend()

        plt.show()


    """
    Example execution: python create_and_analyse_sequence.py -sj 11 -p .1 .3 .05 .6 -id A -v 1
    
    probabilities for experiment:
    .1 .3 .05 .6
    .85 .05 .3 .35
    .475 .175 .175 .475
    
    2 sequences for each probability, IDs A & B
    
    """