import itertools
import numpy as np
import modification_rules as mr
import enzymeanalysis as da
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import enzymeplots as dp
import regex as re

matplotlib.rcParams['pdf.fonttype'] = 42
big_font = {'family' : 'sans-serif', 'size'   : 6}
small_font = {'family' : 'sans-serif', 'size' : 4}
matplotlib.rc('font', **big_font)

def overlap_score(*positions):
    #positions is a list enzyme recognition site positions
    positions = [p[::-1] for p in positions]
    overlaps = [0 for i in range(len(positions))]
    lengths = [len(i.replace("x", "")) for i in positions]
    #go through each position, where "position" represents each of the enzymes at position ind
    for ind, position in enumerate(itertools.zip_longest(*positions, fillvalue="x")):
        #'rs_num' will be the index of the enzyme/rs we're looking at, and enzyme_position is the char at
        # that position for that enzyme/rs
        for rs_num, enzyme_position in enumerate(position):
            #if the enzyme is not particular at that position, we don't care
            if enzyme_position == "x":
                continue
            #otherwise...
            else:
                #for all the other enzymes, track the mismatches between this RS and the others
                others = [ep for i, ep in enumerate(position) if (i != rs_num) &
                                                                 (ep != 'x') &
                                                                 (ep != enzyme_position)]
                #if there are mismatches, reocrd this position as a mismatch for this enzyme
                if len(others):
                    overlaps[rs_num] += 1
    scores = [(l-o) / l for l, o in zip(lengths, overlaps)]
    return np.prod(scores)

def spring(x, k1, k2):
    delta = np.array([1 if xi > 0 else 0 for xi in x])
    
    return np.exp(((-k1*delta*(x**2) - k2*(1-delta)*(x**2))/2000)/(8.314462*310))

def position_score(*positions):
    #positions is a list of lists (one for each enzyme, each sublist is
    # [spacing variant, optimal spacing, core_adjustment, and spacing penalty]
    score = 1
    for enzyme_position in positions:
        variant, optimal, core_adj, penalty = enzyme_position
        spacing = variant.count("x") + core_adj
        delta = spacing - optimal
        score = score * spring(np.array([delta]), penalty[0], penalty[1])[0]
    score = max([0, score])
    return score

def position_constraint(enzyme, core_modification_position, query_position):
    enzyme_motif = mr.core_rules[enzyme]
    motif_modification_position = mr.spacing_rules[enzyme][1]
    motif_query_position = query_position - (core_modification_position - motif_modification_position)
    if motif_query_position < 0 or motif_query_position >= len(enzyme_motif):
        return set(mr.aas)
    else:
        return enzyme_motif[motif_query_position]

def hybrid_core_motif(*modifications, **other):
    """takes a list of modifications and positions and returns a schematic of the core motif
    each modification is a tuple of (enzyme, (position1, position2, position3, .... ))
    where each position is where modification is desired from that enzyme
    """
    core_length = 15 if 'core_length' not in other else other['core_length']
    core = []
    for i in range(core_length):
        position_sets = [set(mr.aas)]
        for enzyme, core_positions in modifications:
            for core_position in core_positions:
                position_sets.append(position_constraint(enzyme, core_position, i))
        core.append(set.intersection(*position_sets))
    return core


def enumerate_cores(core_schematic, offtarget_list=[], bad_list=[], good_list=[]):
    """Function takes a core_schematic, which is a list of lists (or list of sets/tuples).
    Each element of the list is a different position, and each element of the nested
    list/set/tuple is the amino acids allowed at that position. Function enumerates all
    possible discrete sequences from that schematic.
    
    BE CAREFUL with this function. It will enumerate ALL possible cores in RAM. That may
    crash the kernel (or computer) if the list if possible cores is extremely large.
    
    offtarget_list is the modifications list used to generate the core_schematic (passed
                    to hybrid_core_motif), enumerated cores with modification sites at
                    positions other than those specified will be removed
    bad_list is a list of strings or regular expressions that are not allowed. Enumerated
                    cores that match any of these strings/regexs will be removed
    good_list is a list of strings or regular expressions that are REQUIRED. Enumerated
                    cores that do not contain at least one of these strings/regex will be
                    removed. Position of the match does not matter here. If position is
                    important, it should be encoded as a core_rule and incorporated
                    into the core_schematic using hybrid_core_motif.
    """
    if np.product([len(i) for i in core_schematic]) > 1000000:
        print('Caution, there are {} cores to enumerate'.format(np.product([len(i) for i in core_schematic])))
    cores = []
    for c in itertools.product(*core_schematic):
        c = "".join(c)
        to_add=True
        for good in good_list:
            if not re.findall(good, c):
                to_add=False
                break
        if not to_add:
            continue
        for mod, positions in offtarget_list:
            if mod == 'tevp':
                continue
            positions = set([p - mr.spacing_rules[mod][1] for p in positions])
            for hit in re.finditer(mr.core_rules_r[mod], c, overlapped=True):
                if hit.start() not in positions:
                    to_add=False
                    break
        if not to_add:
            continue
        for bad in bad_list:
            if re.findall(bad, c):
                to_add=False
                break
        if not to_add:
            continue
        else:
            cores.append(c)
    return cores

def hybrid_leader(*modifications, **other):
    """takes a list of modifications and positions and returns possible leader schematics
    each modification is a tuple of (enzyme, (position1, position2, position3, .... ))
    where each position is where modification is desired from that enzyme
    """
    max_length = 40 if 'max_length' not in other else other['max_length']
    #This is a list of constraints for each enzyme
    # (enzyme name, optimal distance from RS to start of core, list of core residues to start of motif)
    core_distances_to_motif = []
    #This is a list of all modification positions in long form
    # (enzyme name, optimal distance from RS to start of core, modification position in core)
    all_leader_positions = []
    for enzyme, positions in modifications:
        
        if mr.spacing_rules[enzyme][0] == None:
            continue
        optimal_spacing, core_adjustment, _, _ = mr.spacing_rules[enzyme]
        core_distances_to_motif.append((enzyme, optimal_spacing-core_adjustment,
                                        [pos - core_adjustment for pos in positions]))
        for pos in positions:
            all_leader_positions.append((enzyme, optimal_spacing-core_adjustment, pos - core_adjustment))
        
    rs_position_options = []
    for enzyme, _, _ in core_distances_to_motif:
        enzyme_options = []
        for i in range(max_length+1-len(mr.recognition_sites[enzyme])):
            enzyme_options.append((enzyme, mr.recognition_sites[enzyme] + i*'x'))
        rs_position_options.append(enzyme_options)
    
    enzyme_combos = []
    for enz_set in itertools.product(*rs_position_options):
        enz_set_dictionary = dict()
        for enz, rs_position in enz_set:
            enz_set_dictionary[enz] = rs_position
        enz_set_dictionary['overlap_score'] = overlap_score(*[i[1] for i in enz_set])
        
        spacing_list = []
        for enz, opt, core_adj in all_leader_positions:
            optimal_spacing, core_adjustment, k1, k2 = mr.spacing_rules[enz]
            spacing_list.append((enz_set_dictionary[enz], opt, core_adj, (k1, k2)))
        enz_set_dictionary['position_score'] = position_score(*spacing_list)
        
        enz_set_dictionary['score'] = enz_set_dictionary['overlap_score'] * enz_set_dictionary['position_score']
        enz_set_dictionary['length'] = max([len(i[1]) for i in enz_set])
        
        if enz_set_dictionary['score'] > 0.1:
            enzyme_combos.append(enz_set_dictionary)
    return pd.DataFrame(enzyme_combos)
        
def plot_choose_leader(leader_df, ind=None):
    """function takes a pandas dataframe of different enzyme recognition site spacings and plots them so
    the desired leader can be chosen. In the dataframe, each row is a different combination of enzyme
    RS spacings, where a column is labeled with the enzyme and contains the RS followed by "xxx..." where
    the number of x's is the space to the core
    The dataframe should also have columns corresponding to 'score' and 'length'
   
    if ind is supplied, a plot is not generated, rather the leader dictionary of spacing corresponding
    to the leader at that index value is returned (this can be resolved to a single leader sequence using
    build_leader()
    """
    if ind:
        choice = leader_df.loc[ind]
        return dict(choice.drop(['overlap_score', 'position_score', 'score', 'length']))
    
    fig = plt.figure(figsize=(2,2), dpi=300)
    ax = plt.axes([0,0,1,1])
    
    sns.scatterplot(data=leader_df, x='length', y='score', ax=ax, color='black', s=6)
    max_length = leader_df['length'].max()
    max_length = max_length if max_length % 10 == 0 else (int(max_length/10) + 1) * 10
    ax.set_xlim(10, max_length)
    ax.set_ylim(0, 1)
    dp.format_axis(ax)
    
    inds = leader_df.round(decimals=2).sort_values('score').drop_duplicates(['length'],keep='last').\
                    sort_values('length').drop_duplicates(['score'],keep='first').index
    boundry = leader_df.loc[inds]
    for i, combo in boundry.iterrows():
        ax.annotate('{}'.format(i) , xy=(combo['length'], combo['score']-0.031), fontsize=4, ha='center')
    ax.tick_params(axis='both', which='both', labelsize=4)
    ax.set_xlabel("Leader Length", fontdict=dp.small_font)
    ax.set_ylabel("Score", fontdict=dp.small_font)
    plt.show()
    
def fill_gaps(rs_placing, enzyme, fill=None):
    if fill:
        wt_pre_rs = fill.lower() * len(rs_placing)
        wt_post_rs = fill.lower() * len(rs_placing)
    else:
        wt_pre_rs, wt_post_rs = mr.leaders[enzyme]
    to_fill_pre_rs, to_fill_post_rs = rs_placing.split(mr.recognition_sites[enzyme])

    wt_pre_rs = ('x'*100) + wt_pre_rs
    wt_post_rs = wt_post_rs + ('x' * 100)

    pre_rs = wt_pre_rs[-len(to_fill_pre_rs):] if len(to_fill_pre_rs) else ""
    post_rs = wt_post_rs[:len(to_fill_post_rs)] if len(to_fill_post_rs) else ""

    return pre_rs.lower() + mr.recognition_sites[enzyme] + post_rs.lower()
    
def build_leader(leader_dict, fill='random', prioritize='random'):
    """Function takes a leader dictionary and resolves multiple potentially overlapping recognition
    sites, while also filling in spacer regions with amino acids. The dictionary is formatted as
    key-value pairs, where each key is the name of the enzyme, and the value is the recognition site
    sequence in uppercase, with lower-case 'x' between the recognition sites and the core
        for example "lynd" - LAELSEEALxxxxx, would specify the lynd recognition site with a 5aa spacer
    each of the enzymes, and the recognition sites, need to be present in the recognition_sites
    dictionary in modification_rules.py
    
    other options:
    fill - specifies how to fill spacer regions:
                random - default - choose a random aa at each spacer position
                enzyme - name of an enzyme to fill around the recognition site using wild-type peptide
                            sequence. Any remaining unspecified positions will be filled randomly
                list of enzymes - multiple names of enzymes to fill around their respective
                            recognition sites. Priority of which peptide amino acids to use can be
                            specified using the prioritize option, otherwise it is random
                amino acids - string of amino acids, which will be tiled and inserted in the spacer
                            regions (for example 'ggs' would insert ggsggsggs... linkers)
    prioritize - specifies which enzyme constraints take priority:
                enzyme - name of an enzyme, recognition site amino acids and fill amino acids
                            (if specified) will take priority from this enzyme
                list of enzymes - name of enzymes, in order of importance (most important to least
                            important), to be prioritized
    """
    length = max([len(v) for v in leader_dict.values()])
    for k, v in leader_dict.items():
        leader_dict[k] = ("x" * (length-len(v))) + v

    if fill == 'random':
        pass
    elif (type(fill) == str):
        if (type(fill) in mr.leaders):
            leader_dict[fill] = fill_gaps(leader_dict[fill], fill)
        else:
            leader_dict[k] = fill_gaps(leader_dict[k], k, fill=fill)
    else:
        for f in fill:
            leader_dict[f] = fill_gaps(leader_dict[f], f)
    
    if prioritize == 'random':
        prioritize = []
    elif type(prioritize) == str:
        prioritize = [prioritize]
    elif type(prioritize) == list:
        pass
    else:
        print("prioritize setting is not correct")
        return
    
    positions = []
    for enz, pep in leader_dict.items():
        for ind, aa in enumerate(pep):
            positions.append(dict([('enzyme', enz), ('aa', aa), ('position', ind)]))
    positions = pd.DataFrame(positions)
    positions['priority'] = positions['enzyme'].apply(lambda x: 1000 if x not in prioritize else prioritize.index(x))
    positions['uppercase'] = positions['aa'].apply(lambda x: int(x.islower()))
    positions['wild'] = positions['aa'].apply(lambda x: int((x == 'x') or (x == 'X')))
    random_ordering = np.random.shuffle(list(range(len(positions))))
    positions['random'] = random_ordering
    leader = list(positions.sort_values(['position',
                       'uppercase',
                       'wild',
                       'priority',
                       'random']).drop_duplicates('position', keep='first')['aa'])
    for leader_ind, leader_aa in enumerate(leader):
        if leader_aa == 'x' or leader_aa == 'X':
            leader[leader_ind] = np.random.choice(list(mr.aas.lower()))
    return "".join(leader)

def build_precursor_peptides(leader, cores):
    return [leader + core for core in cores]
    