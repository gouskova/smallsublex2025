#!/usr/bin/env python3

import os, argparse, itertools

'''
a module for phonological feature wrangling and natural class calculations.
fixed the earlier problem of pynatclasses where certain classes were accidentally left out
'''


def check_feats(**kwargs):
    '''
    returns False if the feature specifications of one seg are a proper subset of the other. when this holds, the first seg cannot be uniquely identified using its features.
    prints a message if the file is well-formed and correctly structured
    '''
    segdic = kwargs.get('segdic', make_segdic(**kwargs)['segdic'])
    problemsegs = []
    for seg, otherseg in itertools.combinations(segdic.keys(), 2):
        if segdic[seg].issubset(segdic[otherseg]):
            problemsegs.append((seg, otherseg))
    if not problemsegs:
         kwargs['message'] = "\nThe feature file is well-formed. all the segments can be uniquely identified.\n"
         print(kwargs['message']) 
         return True
    else:
         kwargs['message'] = "\nThe feature file does not allow certain segments to be distinguished from each other:\n\n"
         print(kwargs['message'])
         for x in problemsegs:
             kwargs['message']=f'\n{x[0]}  has a subset of the features of {x[1]}'
             print(kwargs['message']) 
         return False 
 

def make_segdic(**kwargs):
    '''
    requires 'featpath' argument (which is a full path to a feature file)
    returns a dictionary of segments with feature values:
    {'f': {'+labial', '+continuant', ...} 
    and also a list of feature names from the file.
    '''
    with open(kwargs['featpath'], 'r', encoding='utf-8') as f:
        x = f.readlines()
    featnames = x[0].strip().split('\t')
    seglines = [y.strip().split('\t') for y in x[1:]]
    segdic={}
    for line in seglines:
        segdic[line[0]] = set() 
        for feat in featnames:
            if not line[1:][featnames.index(feat)]=='0':
                segdic[line[0]].add(f'{line[1:][featnames.index(feat)]}{feat}')
    kwargs['segdic']=segdic
    kwargs['featnames']=featnames
    return kwargs


def make_featdic(**kwargs):
    '''
    returns a dictionary of all the + and - valued features (not zero-valued) and the segments they extend to
    {+son : {m, n, l, ...}}
    '''
    segdic = kwargs.get('segdic', make_segdic(**kwargs)['segdic'])
    verbosity = kwargs.get('verbosity', 1)
    featdic = {}
    for k in segdic:
        for f in segdic[k]:
            if f in featdic:
                featdic[f].add(k)
            else:
                featdic[f]=set()
                featdic[f].add(k)
    kwargs['featdic']=featdic
    kwargs['segdic']=segdic
    return kwargs

def make_feat_vectors(**kwargs):
    '''
    like make_featdic, only makes its own featnames and segdic on the fly. meant to be called externally
    '''
    kwargs = make_segdic(**kwargs)
    kwargs = make_featdic(**kwargs)
    return kwargs

def segs_to_feats(segdic, segs):
    '''
    given two or more segs, looks up whether they have features in common
    '''
    allset = segdic[segs[0]] 
    for seg in segs:
        allset=allset.intersection(segdic[seg])
    return allset


def feats_to_segs(featdic, feats):
    '''
    takes a feature-to-segment dictionary and a list of features;
    finds segments that have those feature specifications. for example, given "+nasal" an "-syllabic", might return {m, n}
    '''
    allsegs=featdic[feats[0]]
    for feat in feats[1:]:
        allsegs=allsegs.intersection(featdic[feat])
    return allsegs

def powerset(x):
    '''
    powerset recipe from itertools doc.
    powerset (a, b, c):
    (a)
    (b)
    (c)
    (a b)
    (b c)
    (a c)
    (a b c)
    '''
    return(itertools.chain.from_iterable(itertools.combinations(x, r) for r in range(len(x)+1)))


def make_verbose_dic(**kwargs):
    '''
    identifies exhaustively all the natural classes in a feature file, without attempting to find shortest descriptions for the classes
    keys are comma-joined strings of segments; values are lists of feature values (+syll, -son, etc)
    '''
    verbosedic = {}
    segdic = kwargs.get('segdic') 
    featdic = kwargs.get('featdic')
    verbosity = kwargs.get('verbosity', 3)
    for seg in segdic: #every seg is a natural class
        verbosedic[seg]=segdic[seg]
        for tup in powerset(segdic[seg]): #take all the feature specs for the seg and analyze subsets they describe
            if tup!=tuple(): #for any group of features in a seg's description
                segs = sorted(feats_to_segs(featdic, [x for x in tup]))
                if len(segs)>1:
                    verbosedic[','.join(segs)]=segs_to_feats(segdic, segs)
    if verbosity>2:
        print(f'{len(verbosedic)} natural classes')
    kwargs['verbosedic']=verbosedic
    return kwargs 


def avg_cl_size(featdict, featlist):
	'''
	calculates average number of segs that a class refers to
	'''
	lenths = [len(featdict[feat]) for feat in featlist]
	try:
		return sum(lenths)/len(featlist)
	except ZeroDivisionError:	
		return 1 

def compactdic(**kwargs):
    '''
    finds the shortest featural descriptions for each natural class, prioritizing features that refer to the largest natural classes themselves.
    for example, in Russian, all vowels are [-nasal], so an exhaustive feature listing for [a e i o u á é í ó ú] includes [+syllabic, -consonantal, +voice, -nasal, +cont, +son]. since [+syllabic] alone also refers to all and only those segments, this would be a better descriptor for that class. 
    this function is designed to produce the same natural class descriptions each time, however imperfect. Ir prioritizes the shortest of equivalent descriptions, and in the case of a tie, it favors the description that uses the features that refer to the largest classes.
    '''
    if not 'verbosedic' in kwargs:
        kwargs = make_verbose_dic(**make_featdic(**make_segdic(**kwargs)))
        verbosedic = kwargs['verbosedic']
    else:
        verbosedic = kwargs['verbosedic']
    featnames = kwargs.get('featnames')
    featdic = kwargs.get('featdic')
    nclassdic = {}.fromkeys(verbosedic)
    for nclass in verbosedic:
        powerfeats = powerset(verbosedic[nclass])
        good_desc = []
        nset = set(nclass.split(','))
        for fp in powerfeats:
            if fp!=tuple() and fp!=tuple(verbosedic[nclass]) and fp!=("/",):
                kwargs['feats']=fp
                fts = feats_to_segs(featdic, fp)
                if fts==nset and len(fp)<len(verbosedic[nclass]):
                    good_desc.append(fp)
        if len(good_desc)==0:
            nclassdic[nclass]=verbosedic[nclass]
        elif len(good_desc)==1:
            nclassdic[nclass]=set(good_desc[0])
        else:
            lenths = [len(x) for x in good_desc]
            shortest = [x for x in good_desc if len(x)==min(lenths)]
            if len(shortest)==1:
                nclassdic[nclass]=set(shortest[0])
            else:
                cl_sizes = [avg_cl_size(featdic, cl) for cl in good_desc]
                generalest = [x for x in good_desc if avg_cl_size(featdic, x) == min(cl_sizes)]
                nclassdic[nclass]=set(generalest[0])
    del fts
    del kwargs['verbosedic']
    kwargs['nclassdic']=nclassdic
    return kwargs

def featclassdic(**kwargs):
    '''
    takes in a dictionary of segment groups and feature values, and returns the inverse (feature values map to segment lists)
    '''
    if 'nclassdic' in kwargs:
        nclassdic = kwargs.get('nclassdic')
    else:
        nclassdic = compactdic(**kwargs)['nclassdic']
    featclassdic = {}
    for k in nclassdic:
        featclassdic[','.join(sorted(list(nclassdic[k])))] = set(k.split(','))
    kwargs['featclassdic']=featclassdic
    return kwargs

def sclassdic(**kwargs):
    '''
    returns a dictionary of all the natural classes that each seg belongs to. at a minimum needs a feature path; does the rest.
    '''
    if 'nclassdic' in kwargs:
        fdic = kwargs.get('nclassdic')
    else:
        fdic = compactdic(**kwargs)['nclassdic'].values()
    if 'segdic' in kwargs:
        segdic = kwargs.get('segdic')
    else:
        segdic = make_segdic(**kwargs)['segdic']
    segclassdic={}
    for seg in segdic:
        segclassdic[seg] = [cl for cl in fdic if cl.issubset(segdic[seg])]
    kwargs['segclassdic']=segclassdic
    return kwargs

def nclasses(**kwargs):
    '''
    prints natural classes to a file or returns the natural class dictionary 
    '''
    featpath = kwargs.get('featpath')
    outpath = kwargs.get('outpath')
    pr = kwargs.get('print', True)
    kwargs = compactdic(**kwargs)
    outdic = kwargs['nclassdic']
    if pr:
        with open(outpath, 'w', encoding='utf-8') as f:
            for cl in sorted(outdic):
                f.write(f"{cl}\t{','.join(list(outdic[cl]))}\n")
        print(f"{len(outdic)} classes written to {outpath}")
    else:
        return outdic

def get_vocoids(**kwargs):
    '''
    gets vocoids, i.e. vowels and glides (-consonantal or -cons)
    '''
    if 'featdic' in kwargs:
        featdic = kwargs['featdic']
    else:
        featdic = make_featdic(**kwargs)['featdic']
    if '-cons' in featdic:
        return featdic['-cons']
    elif '-consonantal' in featdic:
        return featdic['-consonantal']
    else:
        print("Your feature table does not include [consonantal]. Please fix this and try again.")

def get_consonants(**kwargs):
    '''
    gets consonants! returns a list of -syllabic segments
    '''
    if 'featdic' in kwargs:
        featdic = kwargs['featdic']
    else:
        featdic = make_featdic(**kwargs)['featdic']
    if '-syll' in featdic:
        return featdic['-syll']
    elif '-syllabic' in featdic:
        return featdic['-syllabic']
    else:
        print("Your feature table does not include [syllabic]. Please fix this and try again.")

def get_vowels(**kwargs):
    '''
    gets vowels! returns a list of +syllabic segments
    '''
    if 'featdic' in kwargs:
        featdic = kwargs['featdic']
    else:
        featdic = make_featdic(**kwargs)['featdic']
    if '+syll' in featdic:
        return featdic['+syll']
    else:
        return featdic['+syllabic']

def tightest_class(**kwargs):
    '''
    gets a list of segments, returns the smallest natural class that includes those segments and the smallest number of other segments
    '''
    if type(kwargs['segset'])==str:
        segset = [x.strip(" ") for x in kwargs.get('segset').split(",")]
    elif type(kwargs['segset'])==list:
        segset = kwargs.get('segset')
    if 'nclassdic' in kwargs:
        nclassdic = kwargs.get('nclassdic')['nclassdic']
    else:
        nclassdic = compactdic(**kwargs)['nclassdic']
    posscls = {}
    lenths = []
    for cl in nclassdic:
        if set(segset)==set(cl.split(",")):
            return {cl:nclassdic[cl]}
        elif set(segset).issubset(cl.split(",")):
            posscls[cl]=nclassdic[cl]
            lenths.append(len(cl.split(",")))
    if len(posscls)==0:
        print(f"{segset} do not form a natural class")
        return None
    else:
        tightest = {cl:posscls[cl] for cl in posscls if len(cl.split(","))==min(lenths)} 
        return tightest 
       
def missing_classes(**kwargs):
    '''
    gets a list of segments, returns all the classes missing from the list
    '''
    if type(kwargs['segset'])==str:
        segset = [x.strip(" ") for x in kwargs.get('segset').split(",")]
    elif type(kwargs['segset'])==list:
        segset = kwargs.get('segset')
    if 'nclassdic' in kwargs:
        nclassdic = kwargs.get('nclassdic')['nclassdic']
    else:
        nclassdic = compactdic(**kwargs)['nclassdic']
    posscls = {}
    for cl in nclassdic:
        if len(set(segset) & set(cl.split(",")))==0:
            posscls[cl] = nclassdic[cl]
    outdic = {}
    return posscls
 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="natural classes and phonological feature queries")
    parser.add_argument("--featpath", help="full path to the feature file. Must be tab-separated and use only +, -, and 0 values for features, a la Hayes & Wilson's format")
    parser.add_argument("--ld", help="path to the learning data file. Must contain one word per line, segs separated by spaces. e.g. 'tʃ i p s' for <chips>")
    parser.add_argument("--language", help="partial path to a language folder with learning data and features files. Assumes rigid filename convention, LearningData.txt and Features.txt", default=None)
    parser.add_argument("--nc", help="switch for producing a natural classes dictionary; key is comma-joined feat list and value is comma-joined seg list.")
    parser.add_argument('--verbosity', help="on a scale of 0 to 3, specify how much you want to print to screen. Less means faster performance.", type=int, default=1)
    parser.add_argument("--outpath", help="path to the file where you want the natural classes to be written. Any file by that name will be overwritten without a prompt.", default=os.path.expanduser("~/Desktop/natclasses.txt"))
    parser.add_argument("--segclassdic", help="produce a segment-to-nat-class dictionary", type=bool, default=False)
    parser.add_argument("--segset", help="return the smallest natural class that contains all the segments in a given list", type=str, default=None)
    args = parser.parse_args()
    kwargs = vars(args)
    basepath=os.path.dirname(os.getcwd())
    if kwargs['language']!=None: 
        kwargs['featpath']=os.path.join(basepath, 
                                        'data', 
                                        kwargs['language'], 
                                        'Features.txt')
        kwargs['ld']=os.path.join(basepath, 
                                    'data', 
                                    kwargs['language'], 
                                    'LearningData.txt')
        kwargs['outpath'] = os.path.join(basepath, 
                                            'data', 
                                            kwargs['language'], 
                                            'natclasses.txt')
        if kwargs['segclassdic']==True:
            x = sclassdic(**compactdic(**kwargs))
            with open(kwargs['outpath'].replace('natclasses', 'segclasses'), 'w', encoding='utf-8') as f:
                for seg in x['segclassdic']:
                    for cl in x['segclassdic'][seg]:
                        f.write(f"{seg}\t{cl}\n")
        else:
            nclasses(**kwargs)
    if kwargs['segset']!=None:
        kwargs['featpath']='/home/maria/git/smallsublex/data/russian/Features.txt'
        kwargs['nclassdic']=compactdic(**kwargs)
        tightest = tightest_class(**kwargs)
        print(f"The smallest set including {kwargs['segset']} is {tightest}")
        missing = missing_classes(**kwargs)
        for cl in sorted(missing, key=missing.get, reverse=True):
            print(f"{cl}\t{missing[cl]}")
    else:
        nclasses(**kwargs)
