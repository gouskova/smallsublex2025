#!/usr/bin/env python3

'''
A bunch of utility functions for extracting counts from word lists
'''


import os, sys, numpy, re

shpath = os.path.expanduser('~/git/morphology/russian/input/sharoff_freq.txt')

basedir = os.path.dirname(os.getcwd())
#print(basedir)
subdir = 'data/corpus_searches'


astpath = os.path.join(basedir, subdir, 'search_results_aranea/lemma_astyj.txt')
atpath = os.path.join(basedir, subdir, 'search_results_aranea/lemma_atyj.txt')
bpath = os.path.join(basedir, subdir, 'search_results_aranea/lemma_bases.txt')
rmdpath = os.path.join(basedir, subdir, 'russian_morpheme_dictionary.txt')
yerpath = os.path.join(basedir, subdir, 'yers.txt')
nounpath = os.path.join(basedir, subdir, 'nouns.txt')

def tolerance(N):
    '''
    calculates yang's formula for any N
    '''
    try:
        return N/numpy.log(N)
    except RuntimeWarning:
        return NA

def dumbcounts(**kwargs):
    '''
    dumb because it's not morphology-aware, just looks for whether a word ends in a particular string.
    '''
    count = 0
    outdic={}
    with open(shpath, 'r', encoding='utf-8') as f:
        for line in f:
            line=line.strip().split('\t')
            if line[2].endswith(kwargs['affix']):
                #print(line[2])
                outdic[line[2]]=float(line[1])
                count+=1
            elif kwargs['affix']=='астый' and line[2].endswith('ястый'):
                outdic[line[2]]=float(line[1])
                count+=1
    if kwargs['affix']=='ист':
        for wd in ['лист', 'свист', 'антихрист', 'аист', 'твист', 'посвист']:
            outdic.pop(wd)
            count-=1
    if kwargs['affix']=='астый':
        outdic.pop('частый')
        count-=1
        outdic.pop('нечастый')
        count-=1
        print(outdic)
    print(f'{count} words found with {kwargs["affix"]}')
    print(f'mean occurrences per million: {numpy.mean(list(outdic.values()))}')
    return outdic
            

def find_frequent_ast(toler=True):
    '''
    filters the aranea results so that only the most frequent Sharoff lemmas are kept (in a list of 32000)
    '''
    freqlist = set()
    astlist = set()
    atlist = set()
    baselist = set()
    with open(shpath, 'r', encoding='utf-8') as f:
        for line in f:
            freqlist.add(line.strip().split('\t')[2].strip())
    with open(astpath, 'r', encoding='utf-8') as f:
        for line in f:
            astlist.add(line.strip().split('\t')[0].strip())
    with open(atpath, 'r', encoding='utf-8') as f:
        for line in f:
            atlist.add(line.strip().split('\t')[0].strip())
    with open(bpath, 'r', encoding='utf-8') as f:
        for line in f:
            baselist.add(line.strip().split('\t')[0].strip())
    freqasts = freqlist.intersection(astlist)
    print("astyj")
    N = len(freqasts)
    for x in sorted(freqasts):
            print(x)
    print('\n\n\n\n')
    if toler:
        print(f'{N} words,  \t{tolerance(N)} exceptions allowed')
    freqats = freqlist.intersection(atlist)
    print("atyj")
    N = len(freqats)
    for x in sorted(freqats):
            print(x)
    print('\n\n\n\n')
    if toler:
        print(f'{N} words,  \t{tolerance(N)} exceptions allowed')
    freqbases = freqlist.intersection(baselist)
    print(f'bases\t {len(freqbases)}')
    for x in sorted(freqbases):
            print(x)
    outdic = {'ast':freqasts, 'at':freqats, 'bases':freqbases}
    return outdic

def outwrite_freqs(**kwargs):
    for l in kwargs:
        fname = f'freqsorted_{l}.txt'
        with open(os.path.join(basedir, 'data', fname), 'w') as f:
                for wd in sorted(kwargs[l]):
                    f.write(wd+'\n')


def filter_ost_by_freq(**kwargs):
    '''
    a utility function that checks which of the sharoff ost forms also occur in the russian morpheme dictionary.
    used mostly to add to the RMD so it includes frequent words!
    '''
    freqlist = set()
    rmdlist = set()
    ostlist = set()
    with open(shpath, 'r', encoding='utf-8') as f:
        for line in f:
            wd = line.strip().split('\t')[2].strip()
            if len(wd)>5 and (wd.endswith('ость') or wd.endswith('есть')):
                freqlist.add(wd)
    with open(rmdpath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.split()
            if line[2].endswith('/ость/') or line[2].endswith('/есть/'):
                rmdlist.add(line[0])
    ostlist = rmdlist.intersection(freqlist)
    residue = freqlist - rmdlist
    if len(residue)>0:
        print(sorted(list(residue)))
    else:
        print("everything is in the output set")
    return ostlist

    
def make_frequent_cat_list(**kwargs):
    '''
    genericized 'make frequent nouns' function that takes any category writes to a path or produces a bolean
    '''
    cat = kwargs['cat']
    writ = kwargs['writ']
    outpath=kwargs['outpath']
    if kwargs['freq']:
        fcat = {}
        with open(shpath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[3]==cat:
                    fcat[line[2]]=float(line[1])
    else:
        fcat = set()
        with open(shpath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[3]==cat:
                    fcat.add(line[2])
    if writ:
        with open(outpath, 'w', encoding='utf-8') as f:
            if kwargs['freq']:
                for x in sorted(fcat, key=fcat.get, reverse=True):
                    f.write(f'{x}\t{fcat[x]}\n')
            else:
                for x in sorted(fcat):
                    f.write(x+'\n')
    else:
        return fcat


def get_yers(inlist, monos=False):
    '''
    very specific function to get yer stems for generating -ast forms. inlist is the output of make_frequent_nouns function above.

    protects monosyllabic yer stems from deletion!
    '''
    yerns = {}
    with open(yerpath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().split('\t')
            wd = line[6]
            pl = line[8].replace("'", "")
            syllcount = len(re.findall('[vV]', line[-2]))
            if monos:
                if wd in inlist and syllcount>1:
                    stem = pl.rstrip('иыа')
                    yerns[wd]=stem
            else:
                if wd in inlist:
                    stem=pl.rstrip('иыa')
                    yerns[wd]=stem
    print(len(yerns))
    print('yer stems')
    return yerns

           
def generate_ast_forms(inlist, outlist, bases=False):
    ast = []
    yerstems = get_yers(inlist, monos=True)
    if bases==True:
        outdic = {}
    for wd in inlist:
        #class 2:
        if wd[-1] in ['а', 'я']:
            astform = wd+'стый'
        elif wd[-1] =='и':
            astform = wd[:-1]+'ястый'
        elif wd[-1] == 'ы':
            astform=wd[:-1]+'астый'
        elif wd[-1] in ['е', 'о', 'ё']:
            #class 1b:
            if wd[-2] in ['ь', 'и']:
                astform=wd[:-1]+'ястый'
            else:
                astform=wd[:-1]+'астый'
        elif wd[-1] in ['й','ь']:
            #this is for c-final masc declension only:
            if wd in yerstems:
                astform=yerstems[wd]+'ястый'
            else:
                astform=wd[:-1]+'ястый'
        else:
            if wd in yerstems:
                astform=yerstems[wd]+'астый'
            else:
                astform=wd+'астый'
        if bases==True:
            outdic[astform]=wd
        ast.append(astform)
        if astform.endswith('ястый') and not astform[-6] in ['а', 'я', 'и', 'ы', 'у', 'ю', 'э', 'е', 'о','ё']:
            astform = astform.split('ястый')[0]+'астый'
        if not 'ьа' in astform:
            ast.append(astform)
            if bases==True:
                outdic[astform]=wd
    ast = sorted(list(set([x for x in ast if not 'ьа' in x])))
    with open(outlist, 'w', encoding='utf-8') as f:
        if bases:
            for astform in sorted(outdic.keys()):
                f.write(f'{astform}\t{outdic[astform]}\n')
        else:
            for wd in ast:
                f.write(wd+'\n')

def generate_onok_forms(inlist, outlist, bases=True):
    '''
    note: this has to have so many lex-specific provisions that it might not ever work
    '''
    onok = []
    yerstems = get_yers(inlist, monos=False)
    if bases==True:
        outdic = {}
    for wd in inlist:
        if wd in yerstems: #угорь-угря-угренок, пес-псенок
            if yerstems[wd].endswith('ц') or yerstems[wd].endswith('к'):
                form = yerstems[wd][:-1]+'чонок'
            else:
                form = yerstems[wd]+'енок'
        #most frequent scenario: final C palatalization
        elif wd[-1] in ['а', 'о', 'ы', 'я', 'и', 'ь', 'й', 'у']: # голубь-голубенок, лиса-лисенок
            form = wd[:-1]+'енок'
        else:
            form = wd+'енок'
        #velar palatalizations: k, g, x --> ch, zh, sh
        form = form.replace('кенок', 'чонок')
        form = form.replace('генок', 'жонок')
        form = form.replace('хенок', 'шонок')
        form = form.replace('ценок', 'чонок')
        form = form.replace('чч', 'ч')
        if bases==True:
            outdic[form]=wd
        onok.append(form)
    manuals = ['жеребенок', 'поросенок', 'цыпленок', 'окоренок', 'бельчонок', 'верблюжонок', 'медвежонок', 'опенок', 'курчонок', 'дошколенок', 'зайчонок', 'крольчонок', 'ягненок', 'мальчонок', 'миленок', 'ребятенок', 'индюшонок', 'щенок', 'утенок', 'навальненок', 'негритенок', 'бесененок', 'теленок'] 
    onok=onok+manuals
    onok = sorted(list(set([x for x in onok if not 'ьо' in x])))
    with open(outlist, 'w', encoding='utf-8') as f:
        if bases:
            for form in sorted(outdic.keys()):
                f.write(f'{form}\t{outdic[form]}\n')
        else:
            for wd in ast:
                f.write(wd+'\n')

def generate_ist_izm(inlist, outlist):
    '''
    Shvedova: nouns as bases,
        take off -ик if bolshevik, 
        take off -ор in konservator/konservatizm,
        take off -um in ultimatum/ultimatizm,
        tsia /tsionizm
        take off й, ь (героизм, царизм)
    For -ist and -izm:
        take off adjectival desinence (stankovij / stankovist/stankovizm)
        take off -ichesk (atlanticheskij/atlantist)
        arxaichnij/arxaist/arxaizm
        latinskij/latinist/latinizm
        take off -n (retsidivnij/retsidivizm, akkuratnij/akkuratist/izm)
    '''
    outist = {}
    outizm = {}
    for wd in inlist:
        




def multidumbs(**kwargs):
    affixes = [#'астый', 
            #'ство', 
            #'ник', 
            #'ость', 
            #'ский', 
            #'онок', 
            #'енок'
            #'атый',
            #"ыня",
            #"иня",
            #"изм",
            #"ист",
            #'ский',
            #'онький',
            'ун',
               ]
    for aff in affixes:
        kwargs['affix']=aff
        fname = aff+'.txt'
        with open(os.path.join('/home/maria', 'Desktop', fname), 'w', encoding='utf-8') as f:
            x = dumbcounts(**kwargs)
            for wd in x:
                f.write(wd+'\n')


if __name__=='__main__':
    '''
    to re-generate astyj_auto and atyj_auto files
    '''
    #outpath = sys.argv[1]
    #kwargs={}
#    z = filter_ost_by_freq(**kwargs)
#    with open(outpath, 'w', encoding='utf-8') as f:
#        for line in sorted(z):
#            f.write(line+'\n')
    #kwargs={'writ':True, 'outpath':'/home/maria/Desktop/freq_nouns.txt','cat':'noun', 'freq':True}
#   make_frequent_cat_list(**kwargs)

#    outpath = sys.argv[1]
#    x = make_frequent_cat_list()
#    if "bases" in sys.argv:
#        generate_ast_forms(x, outpath, bases=True)
#    else:
#        generate_ast_forms(x, outpath)


    #else:
    #    x = find_frequent_ast()
    #    outwrite_freqs(**x)
    kwargs={}
    if sys.argv[1]=='multi':
        multidumbs(**kwargs)

#    kwargs = {'writ':False, "cat":"noun", "outpath":None}
#    inlist = make_frequent_cat_list(**kwargs) 
#    outlist ='/home/maria/Desktop/onok_w_bases.txt'
#    generate_onok_forms(inlist, outlist)
