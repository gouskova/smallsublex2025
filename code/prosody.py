#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys
import re
from nltk import ngrams

# should be in the same code directory
import nclasses as pnc

'''


examines a word list and counts up words of different syllable lengths and 
prosodic shapes, by using grids and CV skeleta.

for example:
    pataKA (with final stress) is ssS, or cvcvCV
    PAta (init stress) is Ss, or CVcv
    panDA (final stress) is sS, or cvcCV

requires the custom nclasses module, and relies on the features file containling the following features:
    [syllabic] : for consonant-vowel distinctions (that's how it counts syllables and generates CV skeleta)
    [stress]: for two-way stress oppositions within vocoids
    if there are morpheme boundaries in the data, the feature for the symbols needs to be called [mb].
'''

def count_syllables(**kwargs):
    '''
    returns a dictionary with syllable counts, assuming [syllabic] is in the feature file passed in kwargs['feats']. prints them to screen optionally
    '''
    ld = kwargs.get('ld')
    kwargs=pnc.make_feat_vectors(**kwargs)
    feats = kwargs.get('featdic') 
    vowels = feats['+syllabic']
    lenthdic = {}
    for wd in ld:
        word = wd.split(" ")
        sylls =len([x for x in word if x in vowels])
        if sylls in lenthdic:
            lenthdic[sylls]+=1
        else:
            lenthdic[sylls]=1
    if kwargs.get('printout'):
        for lenth in sorted(lenthdic):
            print(f"{lenth}\t{lenthdic[lenth]}")
    return lenthdic


def make_syllcount_dic(**kwargs):
    '''
    like count_syllables, only returns a dictionary of words augmented with syllable count information
    '''
    ld = kwargs.get('ld')
    kwargs=pnc.make_feat_vectors(**kwargs)
    feats= kwargs.get('featdic')
    vowels = feats['+syllabic']
    outdic = {}
    for word in ld:
        word = word.split(" ")
        sylls = len([x for x in word if x in vowels])
        outdic[' '.join(word)] = sylls
    return outdic

def count_cv_skeleta(**kwargs):
    '''
    returns a dictionary with CV skeleta. needs a [syllabic] feature; if there is a morpheme boundary in the feature file, it must be named 'mb'.
    if you want morpheme boundaries to be ignored in CV skeleton reckoning, pass 'ignore_mb' as a kwargs[] key
    '''
    ld = kwargs.get('ld')
    kwargs = pnc.make_feat_vectors(**kwargs)
    feats = kwargs.get('featdic')
    vowels = feats['+syllabic']
    #mb = list(feats['+mb'])[0]
    consonants = sorted(feats['-syllabic'], reverse=True)
    #to handle palatalization correctly (sorting in reverse ensures palatalized segs get replaced first):
    cvdic = {}
    for word in sorted(ld, reverse=True):
        for cons in consonants:
            word = word.replace(cons, "C")
        for vow in vowels:
            if len(vow)>1:
                word = word.replace(vow, "VV")
            else:
                word = word.replace(vow, "V")
        if word in cvdic:
            cvdic[word]+=1
        else:
            cvdic[word]=1
    if kwargs.get('printout'):
        for word in sorted(cvdic, key=cvdic.get, reverse=True):
            print(f'{word}\t{cvdic[word]}')
    kwargs['inddic']=cvdic
    return kwargs['inddic']


def count_x_grids(**kwargs):
    '''
    returns a dictionary with x grids. Only distinguishes a two-way stress option, and ignores consonants (therefore cannot be weight-sensitive)
    '''
    ld = kwargs.get('ld')
    lower = kwargs.get('ignore_stress') #it's called 'lower' because things will be lowercase (xx, not xXx, etc.). fewer types to count
    kwargs = pnc.make_feat_vectors(**kwargs)
    feats = kwargs.get('featdic')
    vowels = feats['+syllabic']
    stress = feats['+stress']
    xgriddic = {}
    for word in ld:
        word = word.split(" ")
        word_vowels = [x for x in word if x in vowels]
        xgrid = "# "+" ".join(word_vowels)+ " #"
        for v in stress:
            xgrid = xgrid.replace(v, 'X')
        for v in vowels:
            xgrid = xgrid.replace(v, 'x')
        if lower==True:
            xgrid = xgrid.lower()
        if xgrid in xgriddic:
            xgriddic[xgrid]+=1
        else:
            xgriddic[xgrid]=1
    if kwargs.get('slice'):
        syllcount=int(kwargs['slice'])
        for word in sorted(xgriddic, key=xgriddic.get, reverse=True):
            if len(word)>=syllcount:
                print(f'{word}\t{xgriddic[word]}')
    if kwargs.get('printout'):
        for word in sorted(xgriddic, key=xgriddic.get, reverse=True):
            print(f'{word}\t{xgriddic[word]}')
    kwargs['xgriddic']=xgriddic
    kwargs['inddic']=xgriddic
    return kwargs['xgriddic']

def count_cv_grid_ngrams(**kwargs):
    '''
    takes as input a dictionary of xgrids and counts from a lexicon
    returns counts of ngrams (default up to 5)
    '''
    inddic = kwargs['inddic']
    outdic = dict()
    for strw in inddic:
        for i in range(1,5): 
            for x in ngrams(strw.split(" "), i):
                strx = " ".join(x)
                if strx in outdic:
                    outdic[strx]+=inddic[strw]
                else:
                    outdic[strx]=inddic[strw]
    if kwargs.get('printout'):
        if not 'lexics' in kwargs:
            pass
        elif kwargs['lexics']==True:
            print("PROSODIC NGRAMS IN YOUR LEXICON")
        else:
            print("PROSODIC NGRAMS IN YOUR SUBLEXICON")
        try:
            for ngram in sorted(outdic, key=outdic.get, reverse=True):
                print(f'{ngram}\t{outdic[ngram]}')
        except TypeError:
            for ngram in outdic:
                print(f"{ngram}\t{outdic[ngram]}")
    kwargs['ngramdic']=outdic
    return kwargs


def ngram_diff(**kwargs):
    '''
    gets two lexicons to compare, returns xgrid ngrams found in the lexicon but not in the sublexicon (assuming there is a subset relationship). if there isn't a subset relationship, it will print the differences.
    '''
    if kwargs['CV']==True:
        fnc = count_cv_skeleta
    elif kwargs['xgrids']==True:
        fnc = count_x_grids
    dummy = kwargs
    dummy['ld']=dummy['lex']
    dummy['lexics']=True
    lex_ngrams = count_cv_grid_ngrams(**fnc(**dummy))['ngramdic']
    dummy['ld']=dummy['sublex']
    dummy['lexics']=False
    sublex_ngrams = count_cv_grid_ngrams(**fnc(**dummy))['ngramdic']
    ngramdiff = {}.fromkeys(lex_ngrams)
    for k in sublex_ngrams:
        if k in lex_ngrams:
            ngramdiff.pop(k)
        else:
            print(f'{k} is in the xublexicon but not in the lexicon')
            ngramdiff[k]=None
    kwargs['ngram_diff']=ngramdiff
    print(f"There are {len(ngramdiff)} ngrams in the lexicon that do not occur in the sublexicon")
    if kwargs.get('printout'):
        print(f"here are the ngrams that are not in the sublexicon:")
        for k in ngramdiff:
            print(k)
    return kwargs

def lexsublex_pros_ngrams(**kwargs):
    if kwargs['CV']==True:
        fnc = count_cv_skeleta
    elif kwargs['xgrids']==True:
        fnc = count_x_grids
    dummy = kwargs
    dummy['ld']=dummy['lex']
    dummy['lexics']=True
    lex_ngrams=count_cv_grid_ngrams(**fnc(**dummy))['ngramdic']
    dummy['ld']=dummy['sublex']
    dummy['lexics']=False
    sublex_ngrams=count_cv_grid_ngrams(**fnc(**dummy))['ngramdic']
    outdic = {}
    for k in lex_ngrams:
        outdic[k]={'sublex':0, "lex":lex_ngrams[k]}
    for k in sublex_ngrams:
        if k in outdic:
            outdic[k]['sublex']=sublex_ngrams[k]
        else:
            outdic[k]={'lex':0, 'sublex':sublex_ngrams[k]}
    if kwargs['CV']:
        kwargs['lsub_cv_ngrams']=outdic
    if kwargs['xgrids']:
        kwargs['lsub_xgrid_ngrams']=outdic
    return kwargs

if __name__=='__main__':
        import argparse
        basepath = os.path.dirname(os.path.dirname(os.getcwd()))
        parser = argparse.ArgumentParser(description="A command line utility for extracting prosodic profiles from a one-word-per-line data list. Requires a feature file, to be passed as a command line argument. Try 'python3 prosody.py --help' for more.")
        parser.add_argument('--ld', help="full path to learning data file")
        parser.add_argument('--feats', help="full path to feature file")
        parser.add_argument('--outdir', help='full path to the location of output files. Warning: any folder called "simulation" in that location will be overwritten without a prompt.')
        parser.add_argument('--CV', help='produces CV skeleta for each word', type=bool, default=False)
        parser.add_argument('--xgrids', help='counts metrical x grids', type=bool, default=False)
        parser.add_argument('--language', help='if only this argument is specified, the learner will look for an appropriately named folder within "data" (located at the same level as "code") and will run the simulation on the learning data and features files inside that folder. For example, python compseg.py --lang=english/celex/broad runs the learner on the learning data and features inside ../data/english/celex/broad', default=None)
        parser.add_argument("--printout", help="prints to screen", type=bool, default=True)
        parser.add_argument("--slice", help="prints to screen only a subset of dictionary whose syllable count is equal or greater than slice", type=int, default=False)
        parser.add_argument("--lex", help="partial path to the reference lexicon (inside 'data' directory)", default=None)
        parser.add_argument("--sublex", help="partial path to the reference sublexicon", default=None)
        parser.add_argument("--ignore_stress", help="ignore stress when counting x-grids (basically becomes syllable count", default=False)
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
            try:
                kwargs['featdic']=pnc.make_feat_vectors(**kwargs)
                #kwargs['featdict']=pnc.make_featdic(**kwargs)
                print("\ncounting syllables...")
                count_syllables(**kwargs) 
                print("\ncounting CV skeleta...")
                kwargs['cvngrams']=count_cv_skeleta(**kwargs)
                kwargs['inddic']=kwargs['cvngrams']
                print('\ncounting CV ngrams...')
                count_cv_grid_ngrams(**kwargs)
                print('\ncounting X grids...')
                kwargs['inddic'] = count_x_grids(**kwargs)
                print('\ncounting X grid ngrams...')
                kwargs['xgriddic']=count_cv_grid_ngrams(**kwargs)
            except FileNotFoundError:
                print(f"There is no file in the specified location {kwargs['feats']} or {kwargs['ld']}")
        #alternatively counting ngrams
        if args.lex !=None:
            lexpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.lex, 'LearningData.txt')
            sublexpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.sublex, "LearningData.txt")
            kwargs['featpath']=os.path.join(os.path.dirname(lexpath), "Features.txt")
            temp = {'lex':lexpath, "sublex":sublexpath}
            for s in temp:
                kwargs[s]=set()
                with open(temp[s], 'r', encoding='utf-8') as f:
                    for line in f:
                        word = line.strip().split('\t')[0].replace(" | ", " ")
                        if not word in kwargs[s]:
                            kwargs[s].add("# "+word + " #")
                        else:
                            print(f"Your word list has doublets! {word} appears at least twice")
            #ngram_diff(**kwargs)
            kwargs['CV'] = True
            out = lexsublex_pros_ngrams(**kwargs)['lsub_cv_ngrams']
            outpath = os.path.join(os.path.dirname(sublexpath), "cv_ngrams.txt")
            with open(outpath, 'w', encoding='utf-8') as f:
                f.write("NGRAM\tLEX\tSUBLEX\n")
                for d in out:
                    f.write(f"{d}\t{out[d]['lex']}\t{out[d]['sublex']}\n")
            kwargs['CV'] = False
            kwargs['xgrids']=True
            out = lexsublex_pros_ngrams(**kwargs)['lsub_xgrid_ngrams']
            outpath = os.path.join(os.path.dirname(sublexpath), "xgrid_ngrams.txt")
            with open(outpath, 'w', encoding='utf-8') as f:
                f.write("NGRAM\tLEX\tSUBLEX\n")
                for d in out:
                    f.write(f"{d}\t{out[d]['lex']}\t{out[d]['sublex']}\n")

