#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:13:31 2022

@author: sidra
"""

__author__ = Sidra Sohail

# ignore comments

''' completed - first complete this part
UrbTaxa = open('/Users/sidra/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile_R.txt','r').read()
HiekTaxa = open('/Users/sidra/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021_R.txt','r').read()
ChanTaxa = open('/Users/sidra/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChannewtaxafromR12302021_R.txt','r').read()


HiekUrbCombinedTaxa = HiekTaxa + UrbTaxa
HiekUrbChanCombinedTaxa = HiekUrbCombinedTaxa + ChanTaxa

outfile = open('HiekUrbChanCombinedTaxa.txt','w')
outfile.write(HiekUrbChanCombinedTaxa)
outfile.close()
'''

HiekUrbChanCombinedTaxa = open('HiekUrbChanCombinedTaxa.txt','r').read().split('\n')
removethis = 'Taxonomy	"Kingdom"	"Phylum"	"Class"	"Order"	"Family"	"Genus"'
print(HiekUrbChanCombinedTaxa[1:16020].index(removethis)) # 6111+1 = 6112 plus one bc python has 0 indexing
#for i in range(1,len(HiekUrbChanCombinedTaxa)):
    #if HiekUrbChanCombinedTaxa[i] == removethis:
    #HiekUrbChanCombinedTaxa[i].replace(removethis,'')
HiekUrbChanCombinedTaxa[6112] = ''
removethis = 'Taxonomy	Kingdom	Phylum	Class	Order	Family	Genus'
print(HiekUrbChanCombinedTaxa[1:16020].index(removethis)) # 13055+1 = 13056
HiekUrbChanCombinedTaxa[13056] = ''

with open('HiekUrbChan_CombinedTaxa_fin.txt','w') as f:
    for item in HiekUrbChanCombinedTaxa:
        f.write("%s\n" % item)

# outfile = open('HiekUrbChan_CombinedTaxa_fin.txt','w')
# outfile.write(HiekUrbChanCombinedTaxa)
# outfile.close()
