#!usr/bin/env python3

import os, random, sys
import numpy
import scipy.stats
import plotter
import matplotlib.pyplot as plt
from nltk import ngrams

import prosody as pros
import nclasses as pnc
import segments as sgs


'''
A bit of garbage in here that I'm too lazy to clean up for a proceedings :) but the key simulation function is finc_syllcount_monte, which runs a monte carlo simulation for syllable count and keeps track of which final consonants are drawn. To reproduce the results reported in Gouskova 2025 (AMP proceedings), run the following:

$ python lex_comparison.py --lexicon russian/freq_noun_stems --sublexicon russian/freq_astyj --nsamples 100000 --last True

Or replace freq_astyj with freq_izm and freq_ist to do same for the suffixes in the second set of simulations.
'''

def ransample(**kwargs):
    '''
    samples N times from a dictionary of words and tracks how often word length (in syll) exceeds some max size. 
    '''
    random.seed(4)
    lex = kwargs.get('lex')
    samsize = kwargs.get('samsize')
    nsamples = kwargs.get('nsamples')
    maxlenth = kwargs.get('maxsize')
    verbosity = kwargs.get('verbosity', 1)
    sdic = {'max_length':[], 'min_length':[], 'number_over_maxsize':[]}
    for n in range(nsamples):
        wds = random.choices(lex, k=samsize)
        sdic['max_length'].append(max(wds))
        sdic['min_length'].append(min(wds))
        if max(wds)>maxlenth:
            sdic['number_over_maxsize'].append(1)
        else:
            sdic['number_over_maxsize'].append(0)
    return sdic

def finc_syllcount_monte(**kwargs):
    '''
    samples N times from a dictionary of words and tracks how often segments with certain final segments *and* a certain syllable count appear in the sample at the same time.
    needs a lexicon and a sublexicon, a natural class dictionary, and a number of simulations (samples)
    this also tracks how often segments from various natural classes fail to occur in stem-final position, comparing the sublexicon and the reference lexicon.
    '''
    random.seed(55)
    samsize = kwargs.get('samsize')
    lex = kwargs.get('lex') # this is a set of transcribed words
    sublex = kwargs.get('sublex')
    nsamples = kwargs.get('nsamples')
    kwargs['vowels']=pnc.get_vowels(**kwargs)
    fclassdic = pnc.featclassdic(**kwargs['nclassdic'])['featclassdic']
    nclasses = {}
    kwargs['segset'] = [x.strip('# ').split(" ")[-1] for x in sublex] #list of final segs in sublex
    sublexnatclass = list(list(pnc.tightest_class(**kwargs).values())[0])[0]
    if sublexnatclass !=None:
        print(f"Final segs in sublexicon form the natural class {sublexnatclass}")
    for cl in fclassdic:# here we check which classes are absent from the sublexicon:
        if len(fclassdic[cl]&set(kwargs['segset']))==0: #keep track of nat classes absent both from sublex and then later in MC sims:
            nclasses[cl]={'sim':0, 'segs':fclassdic[cl]}
    sylls = set()
    for wd in sublex:
        sylls.add(len([x for x in wd if x in kwargs['vowels']]))
    sublexmaxsyll = max(sylls)
    print(f"The maximum syllable count in the sublexicon is {sublexmaxsyll}")
    outdic = {'lastnclass': {}, 'maxlenth': {}, 'joint':0}
    #drawing random samples now and checking for feat co-occurrence:
    for n in range(nsamples):
        wds = random.choices(lex, k=samsize)
        #collect actual lenths and put them in a dict:
        lenths = [len([x for x in wd if x in kwargs['vowels']]) for wd in wds]
        #collect nat classes and put them in a dict:
        kwargs['segset']=[x.strip("# ").split(" ")[-1] for x in wds]
        for cl in nclasses:
            if len(fclassdic[cl]&set(kwargs['segset']))==0:
                nclasses[cl]['sim']+=1 #means it was absent from that simulation!!
        #check against sublex and enter result in 'joint':
        simmaxsyll = max(lenths)
        simnatclass = list(list(pnc.tightest_class(**kwargs).values())[0])[0]
        if simnatclass in outdic['lastnclass']:
            outdic['lastnclass'][simnatclass]+=1
        else:
            outdic['lastnclass'][simnatclass]=1
        if simmaxsyll in outdic['maxlenth']:
            outdic['maxlenth'][simmaxsyll]+=1
        else:
            outdic['maxlenth'][simmaxsyll]=1
        if simmaxsyll <= sublexmaxsyll and kwargs['featdic'][simnatclass].issubset(kwargs['featdic'][sublexnatclass]):#check simnatclass is a subset of sublexnatclass
            outdic['joint']+=1
    cl_to_remove = []
    for cl in nclasses:
        for xcl in nclasses:
            if cl != xcl and nclasses[cl]['segs'].issuperset(nclasses[xcl]['segs']):
                 cl_to_remove.append(xcl)
    kwargs['finc']=outdic
    kwargs['nclinc']={k:nclasses[k] for k in nclasses if not k in cl_to_remove}
    return kwargs


def analyze_word(**kwargs):
    '''
    return some specified property of the word (e.g., the last segment or a series of ngrams it contains)
    current options include:
    "last": last segment of the word. (string)
    "ngrams": supplied with integer value. returns a list of segments of length N, padded with pound signs for word edges. 
    '''
    vowels = kwargs.get('vowels') # should be a set
    wdsegs = kwargs.get('word').strip("# ").split(' ')
    return {'last':wdsegs[-1], 'syllcount':len([x for x in wdsegs if x in vowels])}    

def run_segsyll_monte(**kwargs):
    '''
    DOES NOT WORK
    runs segment-based simulations on a lexicon and a sublexicon and returns segments that occur in some defined position in the lexicon only; this tells us what is banned from the sublexicon.
    if it's possible to identify a natural class based generalization, this will do it.
    '''
    print(f"Analyzing {kwargs['lexicon']}")
    kwargs['lex_to_analyze']=kwargs['lex']
    lex = natclass_monte(**kwargs)
    print(f"Analyzing {kwargs['sublexicon']}")
    kwargs['lex_to_analyze']=kwargs['sublex']
    sublex = natclass_monte(**kwargs)
    del kwargs['lex_to_analyze']
    if set(lex.keys()).issuperset(set(sublex.keys())):
        not_chosen = set(lex.keys())-set(sublex.keys())
        print("Not chosen:")
        print(not_chosen)
    elif set(sublex.keys()).issuperset(set(lex.keys())):
        print("Is everything ok? The sublexicon has more variety than the lexicon!!")
    else: #they overlap partially, presumably? or don't overlap at all
        not_chosen = set(lex.keys())-set(sublex.keys())
        subrestricted = set(sublex.keys())-set(lex.keys())
        print(f"Your lexicon and sublexicon overlap only partially.\n {not_chosen} occur only in the reference lexicon.\n {subrestricted} occur only in the sublexicon.")

    

def runsim(**kwargs):
    '''
    runs the monte carlo simulation that focuses on syllable count/size and plots the results
    '''
    sim = ransample(**kwargs)
    ci = ci_long(sim['max_length'])
    for k in sim:
        print(f'{k}\t{numpy.mean(sim[k])}')
    print(f'confidence intervals: {ci}')
    plotter.plot_sim_with_ci(sim['max_length'], ci, abline=kwargs.get('maxsize'), fname=kwargs.get('fname'), color=kwargs.get('color'))

def ci_long(inlist):
    '''
    manual method for confidence interval calculation
    '''
    mean = sum(inlist)/len(inlist)
    std = numpy.std(inlist)
    lowerb = mean-1.96*(std/numpy.sqrt(len(inlist)))
    upperb = mean+1.96*(std/numpy.sqrt(len(inlist)))
    return (lowerb, upperb)

def plot_hist(values, ci, bins=30):
    '''
    plots the frequency of your chosen feature's occurrence in a histogram, with confidence intervals marked in red.
    '''
    fig = plt.hist(values, bins)
    plt.axvline(ci[0], color='r')
    plt.axvline(ci[1], color='r')
    plt.show()
    plt.clf()
    plt.close()

def fisher_test(nums):
    x = scipy.stats.fisher_exact(nums)
    print(f'fisher exact test: {x[0]}, p = {x[1]:5f}')
    return x 

def chisq(nums):
    x = scipy.stats.chi2_contingency(nums)
    #x[0] is odds ratio; x[1] is p val; x[2] is degs of freedom
    print(f'X2{x[2]} = {x[0]}, p={x[1]}')
    return x

def compare_dists(**kwargs): 
    '''
    suppose we have a sublexicon where only words that look like x and xx occur.
    and a lexicon where a bigger range of words occurs, but most of them are still x and xx.
    sublexicon dic: {'x': 30, 'xx': 20}
    lexicon dic: {'x': 1000, 'xx': 2000, 'xxx': 300, 'xxxx': 200}
    this function will sample 50 words (length of sublexicon dic: 30 + 20) from the types of "words" that occur in the lexicon dic, 10,000 times. we'll see how often we get a distribution like that in the sublexicon
    '''
    random.seed(5)
    sublex = kwargs.get('sublex')
    lex = kwargs.get('lex')
    lensublex = sum(sublex.values())
    lenlex = sum(lex.values())
    verbosity = kwargs.get('verbosity', 1)
    customnumber=kwargs.get('customnumber')
    if verbosity>0:
        #some basic lexical statistics
        print(f"\nSublexicon: {lensublex}")
        for x in sorted(sublex):
            print(f"{x}\t{sublex[x]}\t{round(100*(sublex[x]/lensublex),1)}%")
        print(f"\nLexicon: {lenlex}")
        for x in sorted(lex):
            print(f"{x}\t{lex[x]}\t{round(100*(lex[x]/lenlex),1)}%")
    findic = {}.fromkeys(lex.keys(), 0)
    nsamples = kwargs.get('nsamples', 100)
    dumblex = []
    poshits = 0
    if customnumber:
        custdix = {}
    #now sample from the summary lexicon: if there are 20 words of length x, the likelihood of hitting that length is 20/length of lex (note, this samples only abstract descriptions of the words, not the words themselves. thus, 'x' not 'cat', or 'xx' and not 'doggy')
    #this re-inflates the dictionary into a list with the same abstract structure.
    for wd in lex:
        for i in range(lex[wd]):
            dumblex.append(wd)
    if verbosity>2:
        print("\nHere are the individual sims\n\n")
    for n in range(nsamples):
        wds = random.choices(dumblex, k=lensublex) #this is the randomly drawn list for this iteration in the simulation
        owds = {}.fromkeys(wds, 0)
        if verbosity>2:
            print(max([len(x.replace(" ","")) for x in owds]))
        #now count how often the list includes the same types of x-grids as in the sublexicon.
        for wd in wds:
            owds[wd]+=1
            findic[wd]+=1
        if verbosity>2:
            for wd in sorted(owds):
                print(f'{wd}\t{owds[wd]}')
            print('\n\n')
        if set(wds)==set(sublex.keys()):
            poshits+=1
        if customnumber:
            if max([len(wd.replace(" ","")) for wd in owds])==customnumber:
                k = ','.join({wd.replace(" ","") for wd in owds})
                if k in custdix:
                    custdix[k]+=1
                else:
                    custdix[k]=1
        # and here, we do a manual/ad-hoc assessment of whether the monte carlo draw results in the same distrib as the extended sublexicon. this requires an extra switch:
    print('\n\n\n')
    ratio = round(poshits/nsamples, 2)
    print(f"Number of draws with the same inventory as sublexicon: {poshits}")
    print(f'Ratio: {ratio}')
    if customnumber:
        print(f'\nNumber of draws with length cap you specified, {customnumber}, is {sum(custdix.values())}')
        if custdix !={}:
            for wd in sorted(custdix):
                print(f'{" ".join(list(wd))}\t\t{custdix[wd]}')
        else:
            print("\nNone")
    print(f"\nThis is how often each word type was chosen:")
    for wd in sorted(findic, key=findic.get, reverse=True):
        print(f'{wd}\t\t{findic[wd]}\t\t{100*round(findic[wd]/sum(findic.values()),2)}%')
    return findic


def ld_process(**kwargs):
    with open(kwargs['ld'], 'r', encoding='utf-8') as f:
        ld = set()
        for line in f:
            ld.add("# " + line.strip() + " #")
    kwargs['ld']=ld
    return kwargs


if __name__=="__main__":
    import argparse
    basepath=os.path.dirname(os.getcwd())
    datapath = os.path.join(basepath, 'data')
    parser = argparse.ArgumentParser(description="Some Monte Carlo simulations and simple stats for assessing what's missing from a sublexicon")
    parser.add_argument("--lexicon", help="the baseline lexicon: a path pointing to a directory in 'data', with LearningData.txt and Features.txt files")
    parser.add_argument("--sublexicon", help="the sublexicon: a path pointing to a directory in 'data' containing LearingData.txt and Features.txt files")
    parser.add_argument("--verbosity", help="set to 1 or higher if you want to print stuff to screen (default=1)", type=int, default=1)
    #parser.add_argument("--", help="", type=bool, default=True)
    parser.add_argument('--nsamples', help="number of simulations (default is 1,000)", nargs="?", type=int, default=1000)
    parser.add_argument("--compare", help="compares sublexicon and lexicon distributions using monte carlo simulations", type=bool, default=False)
    parser.add_argument("--customnumber", help="enter a cap for max number of syllables to compare monte carlo distributions to (e.g., 3)", type=int, default=None)
    parser.add_argument("--plotsims", help="run a monte carlo simulation with a given lexicon and sublexicon, and plot the results", type=bool, default=False)
    parser.add_argument('--last', help="run a monte carlo simulation with a lexicon and sublexicon and count how often segments occur in stem-final position.", type=bool, default=None)
    args = parser.parse_args()
    kwargs=vars(args)
    kwargs['featpath']=os.path.join(datapath, args.lexicon, 'Features.txt')
    kwargs['nclassdic']=pnc.compactdic(**kwargs)
    kwargs['ignore_stress']=False
    if args.lexicon:
        lexpath = os.path.join(datapath, args.lexicon, 'LearningData.txt')
    if args.sublexicon:
        sublexpath = os.path.join(datapath, args.sublexicon, 'LearningData.txt')
    if args.compare:
        kwargs['ld']=sublexpath
        kwargs=ld_process(**kwargs)
        kwargs['sublex']=pros.count_x_grids(**kwargs)
        kwargs['ld']=lexpath
        kwargs=ld_process(**kwargs)
        kwargs['lex']=pros.count_x_grids(**kwargs)
        outdic = compare_dists(**kwargs)
        print(f"running {kwargs['nsamples']} sims")
    if args.plotsims:
        kwargs['ld']=sublexpath
        kwargs=ld_process(**kwargs)
        kwargs['sublex'] = list(pros.make_syllcount_dic(**kwargs).values())
        kwargs['ld']=lexpath
        kwargs = ld_process(**kwargs)
        kwargs['lex']=list(pros.make_syllcount_dic(**kwargs).values())
        kwargs['ld']=None
        kwargs['maxsize']=max(kwargs['sublex'])
        kwargs['samsize']=len(kwargs['sublex'])
        kwargs['fname']='_'.join([args.lexicon.split('/')[1], 'vs', args.sublexicon.split('/')[1]])
        kwargs['color']=False
        runsim(**kwargs)
    if args.last:
        fstuff=pnc.make_featdic(**{'featpath':os.path.join(os.path.dirname(lexpath), 'Features.txt')})
        kwargs['featdic']=fstuff['featdic']
        kwargs['segdic']=fstuff['segdic']
        with open(lexpath, 'r', encoding='utf-8') as f:
            kwargs['lex']=[x.strip() for x in f.readlines()]
        with open(sublexpath, 'r', encoding='utf-8') as f:
            kwargs['sublex']=[x.strip() for x in f.readlines()]
        kwargs['print']=True
        kwargs['samsize']=len(kwargs['sublex'])
        k = finc_syllcount_monte(**kwargs)
        print(f"How often the nat class of last seg and the max syll count were the same in simulation as in the sublexicon:\n{k['finc']['joint']}/{kwargs['nsamples']}")
        print(f"The natural classes of stem-final segments in the simulations:\n")
        for i in sorted(k['finc']['lastnclass'], key=k['finc']['lastnclass'].get, reverse=True):
            print(f"{i}\t{k['finc']['lastnclass'][i]}")
        print(f"Max stem lengths")
        for i in sorted(k['finc']['maxlenth'], key = k['finc']['maxlenth'].get, reverse=True):
            print(f"{i}\t{k['finc']['maxlenth'][i]}")
        print(f"There were {len(k['nclinc'])} natural classes out of {len(kwargs['nclassdic']['nclassdic'])} that did not occur in stem-final position in the sublexicon\n")
        print(f"Here are the nat classes of final segments absent from the sublexicon and number of times they were drawn in MC simulation\n")
        print(f"class\tsegs\tsize\tn_sims_not_drawn\tratio\n")
        for i in k['nclinc']:
            try:
                print(f"{i}\t{','.join(list(k['nclinc'][i]['segs']))}\t{len(k['nclinc'][i]['segs'])}\t{k['nclinc'][i]['sim']}\t{len(k['nclinc'][i]['segs'])*kwargs['nsamples']/k['nclinc'][i]['sim']}")
            except ZeroDivisionError:
                print(f"{i}\t{','.join(list(k['nclinc'][i]['segs']))}\t{len(k['nclinc'][i]['segs'])}\t{k['nclinc'][i]['sim']}\t{kwargs['nsamples']}")

