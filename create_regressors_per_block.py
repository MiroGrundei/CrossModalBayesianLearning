import argparse
from sbl_agents.sbl_cat_dir_LHtoLH_CItoAR import SBL_Cat_Dir
from utils.helpers import *
import scipy.io as sio
from pdb import set_trace as bp


def create_reg_cd(t, model, sub):
    """
    :param subs:
    :param t:
    :return:
    """

    file_dir = 'E:/TESTING/'
    res_dir = file_dir + 'regressors/testing/phseq_2by2TP_' + model + '/' + sub + '/'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    # load empirical sequence of subject
    seq_fn = file_dir + '/post_hoc_sequences/' + sub + '_seq.mat'
    seq = sio.matlab.loadmat(seq_fn)['seq']
    seqlen = seq.shape[0]

    sequence = np.zeros(seqlen)
    for s in range(0, seqlen - 1):
        sequence[s] = seq[s] - 1

    # re-code sequence states according to separate modalities:
    # intensity (0=L, 1=H) or congruency (0=LC, 1=LI, 2=HC, 3=HI)
    tacseq = np.zeros(len(sequence))
    audseq = np.zeros(len(sequence))
    visseq = np.zeros(len(sequence))

    if model == 'LHtoLH':
        # code for intensity:
        taccode = [0, 1, 0, 0, 1, 1, 0, 1]
        audcode = [0, 0, 1, 0, 1, 0, 1, 1]
        viscode = [0, 0, 0, 1, 0, 1, 1, 1]
    elif model == 'CItoAR':
        # code for congruency
        taccode = [0, 2, 1, 1, 3, 3, 0, 2]
        audcode = [0, 1, 2, 1, 3, 0, 3, 2]
        viscode = [0, 1, 1, 2, 0, 3, 3, 2]

    for i in range(0, len(taccode)):
        tacseq[np.where(sequence == i)] = taccode[i]
        audseq[np.where(sequence == i)] = audcode[i]
        visseq[np.where(sequence == i)] = viscode[i]

    # specifics for subjects with fewer blocks
    if int(sub[-2:]) == 23: # or int(sub[-2:]) == 33:
        bmax = 5
    else:
        bmax = 7

    blocklen = 400
    for block in range(1, bmax):
        if block==1:
            indices = np.arange(0, blocklen)
        elif block==2:
            indices = np.arange(blocklen, 2*blocklen)
        elif block==3:
            indices = np.arange(2*blocklen, 3*blocklen)
        elif block==4:
            indices = np.arange(3*blocklen, 4*blocklen)
        elif block==5:
            indices = np.arange(4*blocklen, 5*blocklen)
        elif block==6:
            indices = np.arange(5*blocklen, 6*blocklen)

        tacs = tacseq[indices]
        auds = audseq[indices]
        viss = visseq[indices]
        seq = sequence[indices]

        # generate regressor for different sequences
        # tac
        CD_SBL_temp = SBL_Cat_Dir(tacs, t, model, verbose=True)
        results, alphas = CD_SBL_temp.compute_surprisal(max_T=CD_SBL_temp.T, verbose_surprisal=False)
        # format and save regressors
        results_formatted = {"sequence": results[:, 1],
                             "predictive_surprise": results[:, 2],
                             "bayesian_surprise": results[:, 3],
                             "confidence_corrected_surprise": results[:, 4],
                             "alphas": alphas}
        res_fn = sub + '_tac_reg_tau_{:.4f}_block{:d}'.format(t, block)
        sio.savemat(res_dir + res_fn, results_formatted)

        # aud
        CD_SBL_temp = SBL_Cat_Dir(auds, t, model, verbose=True)
        results, alphas = CD_SBL_temp.compute_surprisal(max_T=CD_SBL_temp.T, verbose_surprisal=False)
        # format and save regressors
        results_formatted = {"sequence": results[:, 1],
                             "predictive_surprise": results[:, 2],
                             "bayesian_surprise": results[:, 3],
                             "confidence_corrected_surprise": results[:, 4],
                             "alphas": alphas}
        res_fn = sub + '_aud_reg_tau_{:.4f}_block{:d}'.format(t, block)
        sio.savemat(res_dir + res_fn, results_formatted)

        # vis
        CD_SBL_temp = SBL_Cat_Dir(viss, t, model, verbose=True)
        results, alphas = CD_SBL_temp.compute_surprisal(max_T=CD_SBL_temp.T, verbose_surprisal=False)
        # format and save regressors
        results_formatted = {"sequence": results[:, 1],
                             "predictive_surprise": results[:, 2],
                             "bayesian_surprise": results[:, 3],
                             "confidence_corrected_surprise": results[:, 4],
                             "alphas": alphas}
        res_fn = sub + '_vis_reg_tau_{:.4f}_block{:d}'.format(t, block)
        sio.savemat(res_dir + res_fn, results_formatted)

        print("Saved in dir: {}, block{:d}".format(res_dir, block))


if __name__ == "__main__":

    # argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--forget_param', action="store", default=str(0.14),
                        type=str,
                        help="Exponentially weighting parameter for memory/posterior updating; Either number or log for logscale taus.")
    parser.add_argument('-model', '--model', action="store", default="SP",
                        type=str,
                        help='Categorical Dirichlet Transition Probability Model (LHtoLH or CItoAR)')
    parser.add_argument('-sj', '--sj', default="1", nargs='+',
                        help='Subject')

    args = parser.parse_args()

    t = args.forget_param
    model = args.model
    subjects = args.sj

    for s in subjects:
        sub = 'sub-{:02d}'.format(int(s))

        # generate regressors for different taus
        if t == 'log':
            n = 4
            base = 10 ** -n
            # taus = np.logspace(1, 0, num=100, base=base)
            taus = [0, 0.0008, 0.0014, 0.0022, 0.0034, 0.0055, 0.0087, 0.0138, 0.0221, 0.0559, 0.0890, 1.0000]
            for tau in taus:
                tau = round(tau, n)
                create_reg_cd(tau, model, sub)
        else:
            tau = float(t)
            create_reg_cd(tau, model, sub)

    """
    How to run:
        python create_dc_regressors.py -t 0.1 -model TP
        python create_dc_regressors.py -t log -model TP
    """