#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os 
from itertools import product
from nltk import ngrams

# should be in the same code directory
import nclasses as pnc

'''
takes in a learning data file and a feature file and counts up segmental ngrams (up to 3), as well as the natural class sequences they correspond to. Thus, given 

p a l i
n e k t u
gg...

it will count up # p a l, p a, p, # p, etc. 
as well as [wb] [lab] [son] [lat], etc

requires the custom nclasses module, and relies on the features file containling the following features:
    if there are morpheme boundaries in the data, the feature for the symbols needs to be called [mb].
    the module will assume that word boundaries are marked with "#", so that should not be used as a segment label in the data.

'''
def count_seg_ngrams(**kwargs):
    '''
    takes as input a dictionary of xgrids and counts from a lexicon
    returns counts of ngrams
    '''
    inddic = kwargs.get('ld')
    outdic = {} 
    rng_bottom = 1
    rng_top = 4 #maximally trigrams
    for strw in inddic:
        for i in range(rng_bottom, rng_top): #can be tweaked later 
            for x in ngrams(strw.split(" "), i):
                strx = " ".join(x)
                if strx in outdic:
                    outdic[strx]+=1
                else:
                    outdic[strx]=1
    if kwargs.get('printout'):
        for ngram in sorted(outdic, key=outdic.get, reverse=True):
            print(f'{ngram}\t{outdic[ngram]}')
    kwargs['seg_ngrams']=outdic
    print(f'segmental {rng_top-1}-grams: {len(outdic)}') 
    return(kwargs)


def make_natclass_ngrams(**kwargs):
    seg_ngrams = kwargs.get('seg_ngrams')
    if 'segclassdic' in kwargs:
        segclassdic=kwargs['segclassdic']
    else:
        segclassdic = pnc.sclassdic(**kwargs)['segclassdic']
    attested_ngrams = {}
    for gram in seg_ngrams:
        subs = gram.split(" ")
        nclist = [segclassdic[i] for i in subs] # get natural class expansions for each seg
        its_grams = [str(x) for x in product(*nclist)] # combine them into ngrams
        for gr in its_grams:
            if gr in attested_ngrams:
                attested_ngrams[gr]+=seg_ngrams[gram]
            else:
                attested_ngrams[gr]=seg_ngrams[gram]
    kwargs['natclass_ngrams']=attested_ngrams
    print(f'natural class ngrams: {len(attested_ngrams)}')
    return kwargs 


def find_seg_diff(**kwargs):
    '''
    gets two lexicons to compare, and returns segmental ngrams found in the lexicon but not in the sublexicon (assuming they're in a subset relationship). if there is no subset relationship, it will tell you.
    '''
    dummy=kwargs
    dummy['ld']=dummy['lex']
    lex_seg_ngrams = count_seg_ngrams(**dummy)['seg_ngrams']
    dummy['ld']=dummy['sublex']
    sublex_seg_ngrams = count_seg_ngrams(**dummy)['seg_ngrams']
    segdiff = {}.fromkeys(lex_seg_ngrams) #assume a subset relationship by default
    for k in sublex_seg_ngrams:
        if k in lex_seg_ngrams:
            segdiff.pop(k) # remove all the sublexicon ngrams, leaving diff
        else:
            print(f"{k} is in the sublexicon but not in the lexicon")
            segdiff[k]=None
    kwargs['seg_diff']=segdiff
    print(f'There are {len(segdiff)} ngrams in the lexicon that do not occur in the sublexicon') 
    return kwargs


def lexsublex_seg_ngrams(**kwargs):
    '''
    gets two lexicons, and returns a dictionary of segmental ngrams with counts in each lexicon 
    '''
    dummy=kwargs
    dummy['ld']=dummy['lex']
    lex_seg_ngrams = count_seg_ngrams(**dummy)['seg_ngrams']
    dummy['ld']=dummy['sublex']
    sublex_seg_ngrams = count_seg_ngrams(**dummy)['seg_ngrams']
    outdic = {}
    for k in lex_seg_ngrams:
        outdic[k] = {'sublex':0, "lex":lex_seg_ngrams[k]}
    for k in sublex_seg_ngrams:
        if k in outdic:
            outdic[k]['sublex']=sublex_seg_ngrams[k]
        else:
            outdic[k]={'lex':0, 'sublex':sublex_seg_ngrams[k]}
    kwargs['lsub_seg_ngrams']=outdic
    return kwargs




if __name__=='__main__':
        import argparse
        basepath = os.path.dirname(os.path.dirname(os.getcwd()))
        parser = argparse.ArgumentParser(description="A command line utility for extracting prosodic profiles from a one-word-per-line data list. Requires a feature file, to be passed as a command line argument. Try 'python3 prosody.py --help' for more.")
        parser.add_argument('--ld', help="full path to learning data file")
        parser.add_argument('--feats', help="full path to feature file")        
        parser.add_argument('--outdir', help='full path to the location of output files. Warning: any folder called "simulation" in that location will be overwritten without a prompt.')
        parser.add_argument('--language', help='if only this argument is specified, the learner will look for an appropriately named folder within "data" (located at the same level as "code") and will run the simulation on the learning data and features files inside that folder. For example, python compseg.py --lang=english/celex/broad runs the learner on the learning data and features inside ../data/english/celex/broad', default=None)
        parser.add_argument("--printout", help="prints to screen", type=bool, default=False)
        parser.add_argument("--do_ngrams", help="count up segmental and natural class ngrams and save them to file", type=bool, default=False)
        parser.add_argument("--find_lex_segdiff", help="find differences between segmental ngrams in the lexicon and the sublexicon", type=bool, default=False)
        parser.add_argument("--lex", help="partial path to a reference lexicon", type=str, default=None)
        parser.add_argument("--sublex", help="partial path to the sublexicon", type=str, default=None)
        parser.add_argument('--countall', help="get all ngram counts for the lexicon and the sublexicon", type=bool, default=False)
        args=parser.parse_args()
        kwargs = vars(args)
        if args.language!=None:
            lgpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.language)
            ld = os.path.join(lgpath, 'LearningData.txt')
            kwargs['ld']=set()
            with open(ld, 'r', encoding='utf-8') as f:
                for line in f:
                    word = line.strip().split('\t')[0]
                    if not word in kwargs['ld']:
                        kwargs['ld'].add("# " + word+ " #")
                    else:
                        print(f"your word list has doublets! {word} appears at least twice")
            kwargs['featpath'] = os.path.join(lgpath, 'Features.txt')
            if 'do_ngrams' != False:
                try:
                    print("\ncounting segmental ngrams...")
                    kwargs = make_natclass_ngrams(**count_seg_ngrams(**kwargs))
                    ngramdic = kwargs['seg_ngrams']
                    outpath=os.path.join(lgpath, 'segngrams.txt')
                    with open(outpath, 'w', encoding='utf-8') as f:
                        for i in sorted(ngramdic, key=ngramdic.get, reverse=True):
                            f.write(f'{i}\t{ngramdic[i]}\n')
                    ngramdic = kwargs['natclass_ngrams']
                    outpath = os.path.join(lgpath, 'nclassngrams.txt')
                    with open(outpath, 'w', encoding='utf-8') as f:
                        for i in sorted(ngramdic, key=ngramdic.get, reverse=True):
                            f.write(f'{i}\t{ngramdic[i]}\n')
                except FileNotFoundError:
                    print(f"There is no file in the specified location {kwargs['feats']} or {kwargs['ld']}")
        if kwargs['find_lex_segdiff']==True or kwargs['countall']==True:
            lexpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.lex, 'LearningData.txt')
            sublexpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.sublex, 'LearningData.txt')
            kwargs['featpath']=os.path.join(os.path.dirname(os.getcwd()), 'data', args.lex, 'Features.txt')
            temp = {'lex':lexpath, 'sublex':sublexpath}
            for s in temp:
                kwargs[s] = set()
                with open(temp[s], 'r', encoding='utf-8') as f:
                    for line in f:
                        word = line.strip().split('\t')[0].replace(" |", "")
                        if not word in kwargs[s]:
                            kwargs[s].add("# " + word + " #")
                        else:
                            print(f"Your word list has doublets! {word} appears at least twice")
            if kwargs['find_lex_segdiff']==True:
                find_seg_diff(**kwargs)
            if kwargs['countall']==True:
                k = lexsublex_seg_ngrams(**kwargs)['lsub_seg_ngrams']
                outpath = os.path.join(os.path.dirname(sublexpath), 'lex_sublex_segmental_ngrams.txt')
                if kwargs['printout']:
                    print(f"NGRAM\tLEX\tSUBLEX")
                    for i in k:
                        outstring = f"{i}\t{k[i]['lex']}\t{k[i]['sublex']}\n"
                        print(outstring)
                else:
                    with open(outpath, 'w', encoding='utf-8') as f:
                        f.write(f"NGRAM\tLEX\tSUBLEX\n")
                        for i in k:
                            outstring = f"{i}\t{k[i]['lex']}\t{k[i]['sublex']}\n"
                            f.write(outstring)
