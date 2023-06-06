from scipy.io import savemat

from sampling.seq_analysis import *
from utils.helpers import *

mpl.rcParams['keymap.save'] = ''


class seq_gen():
    """
    DESCRIPTION:

    """
    def __init__(self, prob_obs_change, verbose):
        # Initialize parameters of sequence generation instance
        self.order = 1
        self.prob_catch = 0
        self.obs_space = 8
        self.regime_space = 0
        self.prob_regime_init = [0.5]
        self.prob_obs_init = np.tile(1/(self.obs_space), (1, self.obs_space))[0]
        self.prob_obs_change = prob_obs_change
        self.prob_regime_change = 0
        self.verbose = verbose

        # Construct the transition matrices based on input parameters
        self.transition_matrices = self.construct_transitions()

    def construct_transitions(self):
        # function constructs transition matrix and checks row-stoch

        # initialize transition matrix: p(j|i) = P_{i,j}
        B_0 = np.zeros((self.obs_space**self.order, self.obs_space))

        p1, p2, p3, p4 = self.prob_obs_change[0], self.prob_obs_change[1], self.prob_obs_change[2], self.prob_obs_change[3]

        # TAV states: S0=000, S1=100, S2=010, S3=001, S4=110, S5=101, S6=011, S7=111
        # Transition Matrix (states by states):
        # p1 p2 p2 p2 0  0  0  0
        # p2 p4 0  0  p3 p3 0  0
        # p2 0  p4 0  p3 0  p3 0
        # p2 0  0  p4 0  p3 p3 0
        # 0  p3 p3 0  p4 0  0  p2
        # 0  p3 0  p3 0  p4 0  p2
        # 0  0  p3 p3 0  0  p4 p2
        # 0  0  0  0  p2 p2 p2 p1

        # first row:
        B_0[0, 0] = p1       # global repetition
        B_0[0, 1:4] = p2     # unisensory change
        B_0[0, 4:7] = 0.     # bimodal change
        B_0[0, 7] = 0.       # global change

        # fill transition matrix with corresponding transitions:
        B_0[1, 0] = B_0[2, 0] = B_0[3, 0] = B_0[4, 7] = B_0[5, 7] = B_0[6, 7] = B_0[7, 4] = B_0[7, 5] = B_0[7, 6] = B_0[0, 1] # p2
        B_0[1, 4] = B_0[1, 5] = B_0[2, 4] = B_0[2, 6] = B_0[3, 5] = B_0[3, 6] = B_0[4, 1] = B_0[4, 2] = B_0[5, 1] = B_0[5, 3] = B_0[6, 2] = B_0[6, 3] = p3
        # global repetitions:
        B_0[7, 7] = B_0[0, 0]
        for i in range(1, 7):
            B_0[i, i] = p4
        # global changes:
        for i in range(1, 8):
            B_0[i, 7 - i] = B_0[0, 7]
        # bimodal changes:
        B_0[1, 2] = B_0[1, 3] = B_0[1, 7] = B_0[2, 1] = B_0[2, 3] = B_0[2, 7] = B_0[3, 1] = B_0[3, 2] = B_0[3, 7] = B_0[4, 0] \
            = B_0[4, 5] = B_0[4, 6] = B_0[5, 0] = B_0[5, 4] = B_0[5, 6] = B_0[6, 0] = B_0[6, 4] = B_0[6, 5] = B_0[7, 1] = B_0[7, 2] = B_0[7, 3] = B_0[0, 4]

        # check row stochastic
        # to check row stochasticity (sum to 1) decimals of transition matrix are rounded to 10th decimal place
        if (np.round(np.sum(B_0, axis=1), 10) != np.ones(B_0.shape[1])).any():
            print("TP matrix: \n {}".format(B_0))
            raise ValueError("Matrices are not row stochastic")

        if self.verbose:
            print("HHMM correctly initialized. Ready to Sample.")
            print("--------------------------------------------")
            if self.order == 1:
                print("TP matrix: \n {}".format(B_0))
        return [B_0]


    def sample(self, seq_length):
        """
        INPUT:
            * seq_length: Length of desired observed sequence
        OUTPUT:
            * sample: (t x 4) array: index, hidden, observed, alternation indicator
        DESCRIPTION:
            1. Sample inital regime and first trial from initial vectors
            2. Loop through desired time steps
        """
        Q = np.zeros((seq_length, 1)).astype(int)

        # Sample first states and observations uniformly
        Q[0:self.order, 0] = np.random.multinomial(self.order, self.prob_obs_init).argmax()

        # Run sampling over the whole sequence
        for t in range(self.order, seq_length):
            idx = Q[t-1, 0]
            obs_sample = np.random.multinomial(1, self.transition_matrices[0][idx, :]).argmax()

            # prevent change following change in sequence: (only if resamp = 1)
            T0 = [0, 2, 3, 6]
            T1 = [1, 4, 5, 7]
            A0 = [0, 1, 3, 5]
            A1 = [2, 4, 6, 7]
            V0 = [0, 1, 2, 4]
            V1 = [3, 5, 6, 7]

            # resample if double change at t
            resamp = 1
            while resamp == 1 and t > 1:
                obs_sample = np.random.multinomial(1, self.transition_matrices[0][idx, :]).argmax()
                resamp = 0

                if obs_sample in T0:
                    if Q[t-1, 0] in T1:
                        if Q[t-2, 0] in T0:
                            resamp = 1
                elif obs_sample in T1:
                    if Q[t-1, 0] in T0:
                        if Q[t-2, 0] in T1:
                            resamp = 1

                if obs_sample in A0:
                    if Q[t-1, 0] in A1:
                        if Q[t-2, 0] in A0:
                            resamp = 1
                elif obs_sample in A1:
                    if Q[t-1, 0] in A0:
                        if Q[t-2, 0] in A1:
                            resamp = 1

                if obs_sample in V0:
                    if Q[t-1, 0] in V1:
                        if Q[t-2, 0] in V0:
                            resamp = 1
                elif obs_sample in V1:
                    if Q[t-1, 0] in V0:
                        if Q[t-2, 0] in V1:
                            resamp = 1

            # enter final sample
            Q[t, 0] = obs_sample

        # Add column with trial/obs/time
        self.sample_seq = np.column_stack((np.arange(seq_length), Q))

        if self.verbose:
            calc_stats(self.sample_seq, self.verbose)
        return self.sample_seq


def save(sequence, seq_gen_temp, matlab_out, stats):

    sequence_meta = {"sample_output": sequence,
                     "prob_regime_init": seq_gen_temp.prob_regime_init,
                     "prob_obs_init": seq_gen_temp.prob_obs_init,
                     "prob_obs_change": seq_gen_temp.prob_obs_change,
                     "prob_regime_change": seq_gen_temp.prob_regime_change,
                     "deviants_analysis": stats["deviants"]}

    if matlab_out:
        savemat(results_dir + title, sequence_meta)
    else:
        save_obj(sequence_meta, results_dir + title)
    print('Saved data and outfiled file')


def sample_and_save(seq_gen_temp, seq_length, title, matlab_out, plot_seq):
    sequence = seq_gen_temp.sample(seq_length)
    stats = calc_stats(sequence, False)

    np.set_printoptions(precision=2)
    np.set_printoptions(suppress=True)

    save(sequence, seq_gen_temp, matlab_out, stats)

    if plot_seq:
        pltlen = 200
        fig, ax = plt.subplots(nrows=5, figsize=(8, 8), ncols=1, dpi=300)
        fig.tight_layout()

        tac = np.zeros(seq_length)
        aud = np.zeros(seq_length)
        vis = np.zeros(seq_length)
        tac[np.where(np.logical_or.reduce((sequence[:, 1] == 1, sequence[:, 1] == 4, sequence[:, 1] == 5, sequence[:, 1] == 7)))] = 1
        aud[np.where(np.logical_or.reduce((sequence[:, 1] == 2, sequence[:, 1] == 4, sequence[:, 1] == 6, sequence[:, 1] == 7)))] = 1
        vis[np.where(np.logical_or.reduce((sequence[:, 1] == 3, sequence[:, 1] == 5, sequence[:, 1] == 6, sequence[:, 1] == 7)))] = 1

        ax[0].scatter(np.arange(pltlen), tac[:pltlen], s=1, c='g')
        ax[0].set_title("Tactile", fontsize="xx-small")

        ax[1].scatter(np.arange(pltlen), aud[:pltlen], s=1, c='y')
        ax[1].set_title("Auditory", fontsize="xx-small")

        ax[2].scatter(np.arange(pltlen), vis[:pltlen], s=1, c='b')
        ax[2].set_title("Visual", fontsize="xx-small")

        # histogram of train lengths
        ax[3].hist(stats["n_tac_trains"], density=False, alpha=0.3, bins=range(1, int(stats["n_tac_trains"].max())), color="g",
                   label=r"Tactile trains ({} > 1)".format(np.sum(stats["n_tac_trains"]>1)))
        ax[3].hist(stats["n_aud_trains"], density=False, alpha=0.3, bins=range(1, int(stats["n_tac_trains"].max())), color="y",
                   label=r"Auditory trains ({} > 1)".format(np.sum(stats["n_aud_trains"]>1)))
        ax[3].hist(stats["n_vis_trains"], density=False, alpha=0.3, bins=range(1, int(stats["n_tac_trains"].max())), color="b",
                   label=r"Visual trains ({} > 1)".format(np.sum(stats["n_vis_trains"]>1)))
        ax[3].legend(ncol=1, fontsize="xx-small")

        # Add extra info as additional subplot with label in legend
        ax[4].plot([], [], ' ', label=r"Avg global train length: {}".format(round(stats["avg_train"], 3)))
        ax[4].plot([], [], ' ', label=r"Global Devs: {} (p={})".format(int(stats["n_dev_glob"]), round((stats["p_dev_glob"]), 3)))
        ax[4].plot([], [], ' ', label=r"Global trains > 1: {}".format(np.sum(stats["n_trains"]>1)))

        ax[4].plot([], [], ' ', label=r"Tactile Devs: {} (p={})".format(int(stats["n_dev_tac"]), round((stats["p_dev_tac"]), 3)))
        ax[4].plot([], [], ' ', label=r"Auditory Devs: {} (p={})".format(int(stats["n_dev_aud"]), round((stats["p_dev_aud"]), 3)))
        ax[4].plot([], [], ' ', label=r"Visual Devs: {} (p={})".format(int(stats["n_dev_vis"]), round((stats["p_dev_vis"]), 3)))

        ax[4].plot([], [], ' ', label=r"Tactile-Auditory Devs: {} (p={})".format(int(stats["n_dev_tac-aud"]), round((stats["p_dev_tac-aud"]), 3)))
        ax[4].plot([], [], ' ', label=r"Tactile-Visual Devs: {} (p={})".format(int(stats["n_dev_tac-vis"]), round((stats["p_dev_tac-vis"]), 3)))
        ax[4].plot([], [], ' ', label=r"Auditory-Visual Devs: {} (p={})".format(int(stats["n_dev_aud-vis"]), round((stats["p_dev_aud-vis"]), 3)))

        ax[4].legend(ncol=3, fontsize="xx-small")
        ax[4].axis('off')

        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-reg_init', '--prob_regime_init', action="store",
                        default=0.5, type=float,
                        help="Initial regime probability")
    parser.add_argument('-reg_change', '--prob_regime_change', action="store",
                        default=0.01, type=float,
                        help="Probability of changing regime")
    parser.add_argument('-obs_init', '--prob_obs_init', action="store",
                        default=0.5, type=float,
                        help="Initial regime probability")
    parser.add_argument('-obs_change', '--prob_obs_change', nargs='+',
                        help="Probability of sampling observations",
                        action="store", type=float)
    parser.add_argument('-catch', '--prob_catch', action="store",
                        default=0.05, type=float,
                        help="Probability of changing regime")
    parser.add_argument('-t', '--title', action="store",
                        default="temporary_sample_title", type=str,
                        help='Title of file which stores sequence')
    parser.add_argument('-seq', '--sequence_length', action="store",
                        default=200, type=int,
                        help='Length of binary sequence being processed')
    parser.add_argument('-matlab', '--mat_file_out',
                        action="store_true",
                        default=True,
                        help='Save output as a .mat file')
    parser.add_argument('-order', '--markov_order', action="store",
                        default=1, type=int,
                        help='Markov dependency on observation level')
    parser.add_argument('-p', '--plot_seq',
                        action="store_true",
                        default=False,
                        help='View/Plot the sampled sequence')
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=True,
                        help='Get status printed out')

    args = parser.parse_args()

    prob_regime_init = np.array([args.prob_regime_init, 1-args.prob_regime_init])
    prob_regime_change = args.prob_regime_change
    prob_obs_init = np.array([args.prob_obs_init, 1-args.prob_obs_init, 0])
    prob_obs_change = args.prob_obs_change
    prob_catch = args.prob_catch

    order = args.markov_order
    seq_length = args.sequence_length
    title = args.title
    matlab_out = args.mat_file_out
    plot_seq = args.plot_seq
    verbose = args.verbose

    gen_temp = seq_gen(prob_obs_change, verbose)

    # sequence = gen_temp.sample(seq_length)
    sample_and_save(gen_temp, seq_length, title, matlab_out, plot_seq)

    """
    
    NOTES: 6 sequences with 600 stim | 3 different probability settings
    
    python seq_gen_trimodal.py -t sub-01_seq_A_1_3_05_6 -seq 600 -obs_change .1 .3 .05 .6
    
         'A_475_175_175_475', ...
         'B_475_175_175_475', ...
         'A_1_3_05_6', ...
         'B_1_3_05_6', ...
         'A_85_05_3_35', ...
         'B_85_05_3_35'
        
    """
