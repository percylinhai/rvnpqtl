#!/usr/bin/python 
#regression-based QTL mapping in extended families (Sham 2002)
#Copied from HE_extend_all_combined_v01.py
#Implement recombination
from __future__ import division
import argparse
import time
import os
import re
import math
from joblib import Parallel, delayed
import Queue
import ccalibd
import random
import numpy as np
import copy
import ibd_m                         #Linhai/others/python/ibd_m.py
import missingparents
from scipy import stats
#################
#parse arguments
parser = argparse.ArgumentParser(description = "Regression-based QTL mapping in extended families")
parser.add_argument('--path',metavar='PATH',required=True, help = "Path for input pedigree information.")
parser.add_argument('--null',action='store_true', help = "Type 1 error analysis without prespecified trait values")
parser.add_argument('--snv', action='store_true',help = 'Calculate on SNV markers')
parser.add_argument('--mname', metavar='Marker Name',help = 'Marker name')
parser.add_argument('--maffile',metavar='FILE',help = 'file containing MAF info for CHP region')
args = parser.parse_args()
maf_by_marker={}
if args.maffile:
    loci_maf = []
    MafFile=file(args.maffile,'r').read()
    MafLines=MafFile.split('\n')
    MafLines.pop()
    for line in MafLines:
        tmp = line.split()
        name=tmp[1]+':'+tmp[2]
        maf_by_marker[name]=float(tmp[3])

def null_generator(founderid,fam,traits):
    #generate permutations by shuffling founder genotypes
    foundergt_list = []
    foundergt = {}
    gpid,gmid = sorted(set(fam.parentsid))[0:2]
    if fam.missing_parent != []:
        gpgt_config = fam.conditional_prob.keys()
        prob=[0]
        accum_prob=0
        for p in gpgt_config:
            accum_prob += fam.conditional_prob[p]
            prob.append(accum_prob)
        randpick = random.uniform(0,1)
        i = 0
        while randpick >= prob[i]:
            i += 1
        gpgt = [int(x) for x in list(gpgt_config[i-1])]
    for f in founderid:
        if fam.missing_parent != [] and f in [gpid,gmid]:
            idx = [gpid,gmid].index(f)
            foundergt_list += gpgt[2*idx:2*idx+2]
        else:
            foundergt_list += fam.fam_dict[f]['gt'][2*fam.mid:2*fam.mid+2]
    alleles = {}
    ibdall = []    #pairwise ibdall
    fnum = len(founderid)
    nfnum = len(fam.nonfounder)
    #Randomly shuffle founder genotypes
    random.shuffle(foundergt_list)
    for idx,f in enumerate(founderid):
        foundergt[f] = foundergt_list[2*idx:2*idx+2]
    #Randomly assign an inheritance vector
    i = random.randint(0,pow(2,2*nfnum)-1)
    tmp = bin(i)[2:]
    inv = '0'*(2*nfnum-len(tmp))+tmp
    alleles = fam.getGT(inv,foundergt)
    conditional_prob = {}
    #Assign missing again and re-impute
    if fam.missing_parent != []:
        mp_pgt=[]
        mp_offgt=[]
        for off in fam.fam_dict[gpid]['offspring']:
            mp_offgt.append(str(alleles[off][0])+str(alleles[off][1]))
        if 'F' in fam.missing_parent:
            mp_pgt.append(None)
            alleles[gpid]=[0,0]
        else:
            mp_pgt.append(str(alleles[gpid][0])+str(alleles[gpid][1]))
        if 'M' in fam.missing_parent:
            mp_pgt.append(None)
            alleles[gmid]=[0,0]
        else:
            mp_pgt.append(str(alleles[gmid][0])+str(alleles[gmid][1]))
        conditional_prob=missingparents.probability(sorted(fam.mp_freq,reverse=True),mp_offgt,mp_pgt)
    ibdall = fam.classify_affect(alleles=alleles,conditional_prob=conditional_prob) #calculate pairwise IBD
    n,d=fam.matrix_cal(alleles,ibdall,traits,conditional_prob)
    return (n,d)


#######################################
class Family:
	def __init__(self,mid):
		self.fam_dict = {}                  #dictionary for stating family relationship
		self.fam=[]                         #original pedigree information
		self.un=[]                          #pairs of uncle-nephew,cousin,sibpair
		self.cousin=[]
		self.sib=[]
		self.hsib=[]
		self.parentsid=[]
		self.nonfounder=[]
                self.parents={}
                self.pairs=[]
		self.ibd_cousin=[]                  #pairwise IBD estimation for cousin pairs
		self.ibd_sib=[]
		self.ibd_hsib=[]
		self.ibd_un=[]
                self.ibdall=[]                  #estimated IBD value based on marker information
                self.null_nd=[]
                self.mid = mid
                #matrix
                self.S = np.matrix([])
                self.D = np.matrix([])
                self.ES = np.matrix([])
                self.ED = np.matrix([])
                self.ET = np.matrix([])
                self.ibdm = np.matrix([])
                self.covSS = np.matrix([])
                self.covDD = np.matrix([])
                self.covSD = np.matrix([])
                self.covTT = np.matrix([])
                self.Y = np.matrix([])
                #inheritance vectors
                self.postinv_prob = []
                self.postibd = {}
                self.prioribd = {}          #list of values of the 1st part in calculate CovTT
                self.prior = []        #sum of prior prob*tau(i,j)*tau(k,l)
                self.traitflag = 0
                self.numerator = 0
                self.denominator = 0
                #family structures
                self.famstruct = {'2nd':0,'3rd':[]}     
                self.info = 1
                #missing parents
                self.conditional_prob={}
                self.missing_parent=[]
                self.mp_freq = []
                self.err = 0  #Mendelian Error
	def setdict(self,family):
		self.fam=family
		fam_num=len(self.fam)
		original_trait=[self.fam[i][4] for i in xrange(fam_num)]
		for i in range(fam_num):        #for each family member
			iid = self.fam[i][0]
			#print 'mid:'+repr(iid)
                        tmp = {'parents':[],'offspring':[],'mate':[],'sex':[],'trait':[],'gt':[]}
			tmp['parents'] = self.fam[i][1:3]
			fid,mid = tmp['parents'][:]       #parents id
                        tmp['sex']=self.fam[i][3]
                        tmp['trait']= self.fam[i][4]
                        if self.traitflag ==0 and set(original_trait) == set([0]):          #Need to specify traits
                            self.traitflag = 1
                        tmp['gt']=self.fam[i][5:]
			#print fid,mid
			if fid != 0:
                            #if nonfounder
                            for pid in [m for m in [fid,mid] if m not in self.fam_dict]:
                                self.fam_dict[pid]={'parents':[],'offspring':[],'mate':[]}
                            #fam_dict={id:{parents:[],offspring:[],mate:[]},id:{},{}}
                            self.fam_dict[fid]['offspring'].append(iid)    
                            self.fam_dict[mid]['offspring'].append(iid)
                            self.fam_dict[fid]['mate']=mid
                            self.fam_dict[mid]['mate']=fid
                            #print self.fam_dict[fid]
                            #print self.fam_dict[mid]
                        if iid in self.fam_dict:
                            self.fam_dict[iid]['parents']=tmp['parents']
                            self.fam_dict[iid]['sex']=tmp['sex']
                            self.fam_dict[iid]['trait']=tmp['trait']
                            self.fam_dict[iid]['gt']=tmp['gt']
                        else:
			    self.fam_dict[iid]=tmp
                        #print self.fam_dict[iid]
                gpid = min(self.fam_dict.keys())
                # get the nonfounders       
		second = [s for s in self.fam_dict[gpid]['offspring']]
                for i in sorted(second, key = lambda x: len(self.fam_dict[x]['offspring'])):
                    # get the 2nd generation individuals in order of their offspring size
                    self.nonfounder.append(i)
		    for o in self.fam_dict[i]['offspring']:
			self.nonfounder.append(o)
                # set value for self.pairs
                nfcount = len(self.nonfounder)
                for i in range(nfcount):
			for j in range(i+1,nfcount):   #for each pair of affecteds
				aid = self.nonfounder[i]
				bid = self.nonfounder[j]
                                self.pairs.append([aid,bid])
                paircount = len(self.pairs)
                self.prior = [[0 for x in range(paircount)] for x in range(paircount)]
		for n,mid in enumerate(self.nonfounder):   
		        self.parents[mid]=self.fam_dict[mid]['parents']	
                # famstruct={'2nd':INTEGER,'3rd':[int,int]}
                tmp = []
                for i in self.nonfounder:
                        if self.fam_dict[i]['parents'][0] == gpid:
                            self.famstruct['2nd'] +=1
                            tmp.append(len(self.fam_dict[i]['offspring']))
                self.famstruct['3rd']=sorted(tmp)
                foundergt = []
                for iid in self.fam_dict.keys():
                    if self.fam_dict[iid]['parents'] == [0,0]:
                        foundergt += self.fam_dict[iid]['gt'][2*self.mid:2*self.mid+2]
                for allele in set(foundergt):
                    self.info *= foundergt.count(allele)/len(foundergt)


	def firstparent(self,parent,gp,gm):                           
            #make the first parent the important one 
		if self.fam_dict[parent[1]]['parents'] == [gp,gm]:        
                #if mother is the offspring of grandparents
			parent.reverse() #put the important parent in the first

	def dcousin(self,ap,bp,gp,gm):                          #check if they are cousin
		if ap != [0,0] and ap != [gp,gm] and bp != [0,0] and bp !=[gp,gm] \
				and ap != bp:                               
                #if the affecteds are 1st cousin
			return True
		else:
			return False

	def dun(self,aid,bid,ap,bp,gp,gm):                      
            #check if it is uncle-nephew pair
		i = False
		if ap !=[0,0] and bp !=[0,0] and bid not in ap and aid not in bp:
			if ap == [gp,gm] and bp != [gp,gm] or ap != [gp,gm] and bp == [gp,gm]:
				i = True
		return i		

	def dsib(self,ap,bp):                       			#check if it is sibpair
		i = False
		if ap == bp and ap != [0,0]:
			i = True
		return i

	def dhsib(self,ap,bp):
		i = False
		if ap != bp and len(set(ap)&set(bp)) > 0:
			i = True
		return i

        def conv(self,dic):
            #combine the inheritance vectors of nonfounders
            d = copy.copy(dic)
            k=d.keys()
            result = []
            if len(k) > 2:
                i = min(sorted([self.nonfounder.index(key) for key in k]))   
                #in the order in nonfounder
                pred=d[self.nonfounder[i]]
                del d[self.nonfounder[i]]
                con = self.conv(d)
            elif len(k) == 2:
                i,j = sorted([self.nonfounder.index(key) for key in k])
                pred = d[self.nonfounder[i]]
                con = d[self.nonfounder[j]]
            for ele1 in pred:
                for ele2 in con:
                    result.append(ele1+ele2)
            return result

	def cal_ibd(self,v):          
            #calculate IBD among nonfounders under a specific inheritance vector
            #assuming fully informative markers
                ibdinv=[]
                ibdinv=ccalibd.apply(self.parents,v,self.nonfounder)
		return ibdinv         #pairwise IBD

        def inherit_vec(self,alleles={}):
            #calculate the posterior inheritance vectors and their corresponding probability
            tmpv = {}
            v = {}
            pos={}
            offspring=0
            postinv = []
            gpid,gmid = sorted(set(self.parentsid))[0:2]
            for nf in self.nonfounder:    #for each nonfounder
                #print 'nf:'+repr(nf)
                fid,mid = self.fam_dict[nf]['parents']
                pgeno = alleles[fid]+alleles[mid]  #parents genotype 
                if fid == gpid:                      #if 2nd generation
                    v[nf] = []
                    tmpv[nf] = []               #inheritance vector for small nuclear fams
                    offspring=0
                tmpos=[]
                pos[nf]=[]
                for allele in alleles[nf]:
                    tmpos.append([i for i,x in enumerate(pgeno) if x == allele])
                for p in tmpos[0]:
                    for q in tmpos[1]:
                        if p<2 and q>1 or p>1 and q<2:
                            if [p,q] not in pos[nf] and [q,p] not in pos[nf]:
                                pos[nf].append([p,q])     #which allele is from which position [0,1,2,3]
                if fid == gpid:                     #if the nonfounder is the 2nd generation
                    for element in pos[nf]:
                        inh = [ele%2 for ele in sorted(element)]   #inheritance vector
                        if inh not in tmpv[nf]:
                            tmpv[nf].append(inh) 
                    v[nf] = copy.copy(tmpv[nf])
                else:                                   #if it is 3rd generation nonfounder
                    offspring +=1
                    if self.fam_dict[fid]['parents'] != [0,0]:  #father important
                        for element in pos[nf]:
                            pa, ma = sorted(element)[:]      #paternal & maternal allele
                            for idx,fpos in enumerate(pos[fid]):
                                if fpos[pa] < 2:             #if paternal allele is from grandpa
                                    tmpinv=[0,ma%2]
                                elif fpos[pa] > 1:
                                    tmpinv=[1,ma%2]
                                for cv in [vi for vi in v[fid] if vi[0:2] == tmpv[fid][idx] and len(vi)==2*offspring]:
                                    v[fid].append(cv+tmpinv)
                        for mv in [vi for vi in v[fid] if len(vi)<2*offspring+2]:    #remove redundant vectors
                            v[fid].remove(mv)
                        #print v
                    elif self.fam_dict[fid]['parents'] == [0,0]:    #if mother important
                        for element in pos[nf]:
                            pa, ma = sorted(element)[:]
                            for idx,mpos in enumerate(pos[mid]):
                                if mpos[ma-2] < 2:             #if paternal allele is from grandpa
                                    tmpinv=[pa%2,0]
                                elif mpos[ma-2] > 1:
                                    tmpinv=[pa%2,1]
                                for cv in [vi for vi in v[mid] if vi[0:2] == tmpv[mid][idx] and len(vi)==2*offspring]:
                                    v[mid].append(cv+tmpinv)
                        for mv in [vi for vi in v[mid] if len(vi)<2*offspring+2]:
                            v[mid].remove(mv)
            #print tmpv
            #print v
            list_postinv = self.conv(v)
            for tmp_list_inv in list_postinv:
                postinv.append(''.join([str(x) for x in tmp_list_inv]))
            return postinv
        
        def combine_inherit_vec(self,m=-1,alleles={},conditional_prob={}):
            #calculate the overall posterior inv probability
            combined_postinv_prob = {}
            combined_postibd = {}
            gpid,gmid = sorted(set(self.parentsid))[0:2]
            if m != -1:
                for iid in self.fam_dict.keys():
                    alleles[iid] = self.fam_dict[iid]['gt'][2*m:2*m+2]
                conditional_prob=self.conditional_prob
            if self.missing_parent != []:
                #if there is missing grandparents
                for gpgt in conditional_prob.keys():
                    # reassign genotypes for grandparents
                    alleles[gpid] = [int(x) for x in gpgt[:2]]
                    alleles[gmid] = [int(x) for x in gpgt[2:]]
                    tmp_postinv=self.inherit_vec(alleles)
                    len_postinv = len(tmp_postinv)
                    for tmpinv in tmp_postinv:
                        #for each posterior inv
                        #calculate the weighted over-all probability
                        if not combined_postinv_prob.has_key(tmpinv):
                            combined_postinv_prob[tmpinv] = 1/len_postinv*conditional_prob[gpgt]
                        else:
                            combined_postinv_prob[tmpinv] += 1/len_postinv*conditional_prob[gpgt]
            else:
                postinv = self.inherit_vec(alleles)
                for tmpinv in postinv:
                    combined_postinv_prob[tmpinv] = 1/len(postinv)
            if self.info != 1:
                #if informative
                for tmpinv in combined_postinv_prob.keys():
                    list_inv = [int(x) for x in tmpinv]
                    combined_postibd[tmpinv] = self.cal_ibd(list_inv)
            if m != -1:
                self.postinv_prob = combined_postinv_prob
                self.postibd = combined_postibd
            else:
                return combined_postinv_prob, combined_postibd



        def null_inv(self):
            ##get the pairwise IBD under each possible inheritance vector
            nfnum = len(self.nonfounder)
            pairnum = len(self.pairs)
            for i in xrange(pow(2,2*nfnum)):
                tmp = bin(i)[2:]
                tmpinv = '0'*(2*nfnum-len(tmp))+tmp
                inv = [int(ele) for ele in tmpinv]         #prior inheritance vectors
                self.prioribd[i]=self.cal_ibd(inv)      #pairwise IBD under each inv
            #print len(self.prioribd) 
            #########################
                for idx1 in range(pairnum):
                    for idx2 in range(idx1,pairnum):
                        self.prior[idx1][idx2] += self.prioribd[i][idx1]*self.prioribd[i][idx2]/(2**(2*nfnum))
                        self.prior[idx2][idx1] = self.prior[idx1][idx2]

        #classify affecteds and calculate the IBD allele number
        def classify_affect(self,m=-1,alleles={},conditional_prob={}):                   
                #m: marker index
                self.parentsid = [p for nf in self.nonfounder for p in self.fam_dict[nf]['parents']]      
                #print parentsid
                ibdall = []
                gpid,gmid = sorted(set(self.parentsid))[0:2]
                #print gpid,gmid
                nfnum = len(self.nonfounder)
                if m != -1:
                    alleles = {}
                    for iid in self.fam_dict.keys():
                        alleles[iid]=self.fam_dict[iid]['gt'][2*m:2*m+2]
                    ##Dealing with missing parents
                    gpgt = alleles[gpid]+alleles[gmid]
                    tmp_offgt=[]
                    mp_offgt=[]
                    for off in self.fam_dict[gpid]['offspring']:
                        tmp_gt=alleles[off]
                        tmp_offgt+=tmp_gt
                        mp_offgt.append(str(tmp_gt[0])+str(tmp_gt[1]))
                    if 0 in gpgt and 0 not in tmp_offgt:
                        #if there is missing grandparent
                        mp_pgt=[]
                        if 0 in gpgt[:2]:
                            mp_pgt.append(None)
                            self.missing_parent.append('F')
                        else:
                            mp_pgt.append(str(gpgt[0])+str(gpgt[1]))
                        if 0 in gpgt[2:]:
                            mp_pgt.append(None)
                            self.missing_parent.append('M')
                        else:
                            mp_pgt.append(str(gpgt[2])+str(gpgt[3]))
                        if self.mp_freq == []:
                            raise Exception('marker frequency info missing')
                        self.conditional_prob=missingparents.probability(sorted(self.mp_freq,reverse=True),mp_offgt,mp_pgt)
                        conditional_prob = self.conditional_prob
                #########
                for i in range(nfnum):
                    if self.err == 0:
                        for j in range(i+1,nfnum):   
                            #for each pair of affecteds
                            aid = self.nonfounder[i]
                            bid = self.nonfounder[j]
                            #print aid,bid
                            enroll = 0
                            ap = copy.copy(self.fam_dict[aid]['parents'])
                            bp = copy.copy(self.fam_dict[bid]['parents'])
                            mean_ibd = 0
                            if self.dcousin(ap,bp,gpid,gmid):       # if a and b are cousin
                                enroll = 1
                                self.firstparent(ap,gpid,gmid)
                                self.firstparent(bp,gpid,gmid)
                                #print bp,bid
                                if self.missing_parent != []:
                                    for gpgt in conditional_prob.keys():
                                        tmp_gpgt=[int(x) for x in list(gpgt)]
                                        genotype = tmp_gpgt+alleles[ap[0]]+alleles[ap[1]]\
                                                +alleles[aid]+alleles[bp[0]]+\
                                                alleles[bp[1]]+alleles[bid]
                                        ibd = ibd_m.cousin_ibd(genotype) 
                                        try:
                                            mean_ibd+=ibd*conditional_prob[gpgt]
                                        except:
                                            self.err = 1
                                            break
                                else:
                                    genotype = alleles[gpid]+alleles[gmid]+\
                                        alleles[ap[0]]+alleles[ap[1]]+alleles[aid]+\
                                        alleles[bp[0]]+alleles[bp[1]]+alleles[bid]
                                        #GT for [grandparents,parents,kid]
                                    #calculate IBD for each cousin pair
                                    mean_ibd = ibd_m.cousin_ibd(genotype)
                                    if mean_ibd is None:
                                        print "Mendelian Inconsistency"
                                        self.err = 1
                                        break
                                if m != -1:
                                    self.cousin.append([aid,bid])
                                    self.ibd_cousin.append(mean_ibd)	
                            elif self.dun(aid,bid,ap,bp,gpid,gmid):                 
                                # uncle-nephew
                                    enroll = 1
                                    uid = 0
                                    nid = 0
                                    if bp == [gpid,gmid]:                                                                           #put uncle at first
                                        uid = bid
                                        nid = aid
                                    else:
                                        uid = aid
                                        nid = bid
                                    unp = copy.copy(self.fam_dict[uid]['parents'])
                                    nep = copy.copy(self.fam_dict[nid]['parents'])
                                    self.firstparent(nep,gpid,gmid)
                                    if self.missing_parent != []:
                                        for gpgt in conditional_prob.keys():
                                            tmp_gpgt=[int(x) for x in list(gpgt)]
                                            genotype = tmp_gpgt+alleles[uid]+alleles[nep[0]]+\
                                                    alleles[nep[1]]+alleles[nid] 
                                            ibd = ibd_m.un_ibd(genotype) 
                                            try:
                                                mean_ibd+=ibd*conditional_prob[gpgt]
                                            except:
                                                self.err = 1
                                                break 
                                    else:
                                        genotype = alleles[gpid]+alleles[gmid]+\
                                            alleles[uid]+alleles[nep[0]]+\
                                            alleles[nep[1]]+alleles[nid] 
                                            #grandparents+uncle+parents+kid
                                        mean_ibd = ibd_m.un_ibd(genotype)
                                        #calculate IBD for each uncle-nephew pair 
                                        if mean_ibd is None:
                                            print "Mendelian Inconsistency"
                                            self.err = 1
                                            break
                                    if m != -1:
                                        self.un.append([uid,nid])
                                        self.ibd_un.append(mean_ibd)
                            elif self.dsib(ap,bp):          # sibpair
                                    enroll = 1
                                    if ap == [gpid,gmid] and self.missing_parent != []:
                                        # if there is missing parent
                                        for gpgt in conditional_prob.keys():
                                            tmp_gpgt=[int(x) for x in list(gpgt)]
                                            genotype = tmp_gpgt+alleles[aid]+alleles[bid]
                                            ibd = ibd_m.sib_ibd(genotype)  #calculate IBD for each cousin pair
                                            try:
                                                mean_ibd+=ibd*conditional_prob[gpgt]
                                            except:
                                                self.err = 1
                                    else:
                                        genotype = alleles[ap[0]]+alleles[ap[1]]+\
                                            alleles[aid]+alleles[bid]
                                        mean_ibd = ibd_m.sib_ibd(genotype)                                                               #calculate IBD for each sibpair
                                        if mean_ibd is None:
                                            self.err = 1
                                            break
                                    if m != -1:
                                        self.sib.append([aid,bid])
                                        self.ibd_sib.append(mean_ibd)
                            elif self.dhsib(ap,bp):                             
                                # half-sib
                                    enroll = 1
                                    sp = set(ap)&set(bp)                                                                    #shared parent between half-sibs
                                    ap.remove(sp)
                                    bp.remove(sp)
                                    genotype = alleles[sp[0]]+alleles[ap[0]]+\
                                            alleles[aid]+alleles[bp[0]]+alleles[bid]                                             #[shared parent, non-share parent, kid]
                                    mean_ibd = ibd_m.hsib_ibd(genotype)
                                    if m != -1:
                                        self.hsib.append([aid,bid])
                                        self.ibd_hsib.append(mean_ibd)
                            else:             #parent-offspring
                                enroll = 2
                                mean_ibd = 1
                            if enroll > 0:
                                ibdall.append(mean_ibd/2)
                if m != -1:
                    self.ibdall=ibdall     #proportion of alleles sharing IBD
                else:
                    return ibdall


        def corre(self,pair):       
            #return the correlation of traits between relative pairs
            aid, bid = pair
            r = 0
            h = 0.5      #estimated heritability
            if [aid,bid] in self.cousin or [bid,aid] in self.cousin:
                r = h*2*1/16
            elif [aid,bid] in self.sib or [bid, aid] in self.sib:
                r = h*2*1/4
            elif [aid,bid] in self.un or [bid, aid] in self.un:
                r = h*2*1/8
            elif [aid,bid] in self.hsib or [bid, aid] in self.hsib:
                r = h*2*1/8
            elif aid == bid:    #self
                r = h*2*1/2
            else:               #parent-offspring
                r = h*2*1/4
            return r

        
        def expect(self,pair):       #return the expected IBD between a relative pair
            aid, bid = pair
            tau = 0
            if [aid,bid] in self.cousin or [bid,aid] in self.cousin:
                tau = 0.125
            elif [aid,bid] in self.sib or [bid, aid] in self.sib:
                tau = 0.5
            elif [aid,bid] in self.un or [bid, aid] in self.un:
                tau = 0.25
            elif [aid,bid] in self.hsib or [bid, aid] in self.hsib:
                tau = 0.25
            else:       #parent-offspring
                tau = 0.5
            return tau

        def setTrait(self,mod=1):
            if mod == 0: #null 
                nulltraits = {}
            ###STEP 1: set traits for founders######
            for iid in self.fam_dict.keys():
                if self.fam_dict[iid]['parents'] == [0,0]:
                    #founders
                    if mod == 1:
                        self.fam_dict[iid]['trait'] = float('%.3f'%random.gauss(0,1))
                    elif mod == 0:
                        nulltraits[iid] = float('%.3f'%random.gauss(0,1))
            ###STEP 2: set traits for nonfounders from multivariate normal distribution####
            cov = np.identity(len(self.nonfounder))
            #Variance-Covariance matrix for traits
            for aid,bid in self.pairs:
                r = self.corre([aid,bid])
                idxa,idxb=self.nonfounder.index(aid),self.nonfounder.index(bid)
                cov[idxa,idxb] = cov[idxb,idxa] = r
            mean = [0]*len(self.nonfounder)
            traits = np.random.multivariate_normal(mean,cov)
            for idx,t in enumerate(traits):
                if mod == 1:
                    self.fam_dict[self.nonfounder[idx]]['trait'] = float('%.3f'%t)
                elif mod == 0:
                    nulltraits[self.nonfounder[idx]] = float('%.3f'%t)
            if mod == 0:
                return nulltraits
    
        def setmatrix(self): 
            ##set the values for the matrix S, D and IBDm
            s = []
            es = []  #expectation of S
            d = []
            ed = []  #expectation of D
            tau = [] 
            for aid,bid in self.pairs:            
                #for each pair of relatives (include parent-offspring)
                trait1 = self.fam_dict[aid]['trait']
                trait2 = self.fam_dict[bid]['trait']
                s.append((trait1+trait2)**2)        #squared sum of traits
                d.append((trait1-trait2)**2)        #squared difference of traits
                r = self.corre([aid,bid])
                es.append(2*(1+r))
                ed.append(2*(1-r))
                tau.append(self.expect([aid,bid]))  #expected IBD between a relative pair
            self.S = np.matrix(s)
            self.D = np.matrix(d)
            self.ES = np.matrix(es)
            self.ED = np.matrix(ed)
            self.ibdm = np.matrix(self.ibdall)
            self.ET = np.matrix(tau)


        def setcov(self):
            ##set the covariance matrix
            pairnum=len(self.pairs)
            nfnum = len(self.nonfounder)
            ss = np.zeros((pairnum,pairnum))
            dd = np.zeros((pairnum,pairnum))
            sd = np.zeros((pairnum,pairnum))
            tt = np.zeros((pairnum,pairnum))
            #for idx1, pair1 in enumerate(self.pairs):
            #    for idx2, pair2 in enumerate(self.pairs):
            for idx1 in range(pairnum):
                for idx2 in range(idx1,pairnum):
                    ##get the covariance matrix for tau(i,j),tau(k,l)
                    #print pair1,pair2
                    pair1 = self.pairs[idx1]
                    pair2 = self.pairs[idx2]
                    prior,post = 0,0
                    prior_expect = 0
                    prior_expect=self.expect(pair1)*self.expect(pair2)
                    prior = self.prior[idx1][idx2]
                    if self.info == 1: #if posterior inv equals to priterior (marker not informative)
                        post = prior
                    else:
                        for inv in self.postinv_prob.keys():
                        #print inv
                        #print self.postibd[ivx][idx1],self.postibd[ivx][idx2]
                            post += self.postibd[inv][idx1]*self.postibd[inv][idx2]*\
                                    self.postinv_prob[inv]
		    #covariance for IBD of two pairs
                    tt[idx2,idx1]=tt[idx1,idx2] = (prior-prior_expect)-\
                            (post-self.ibdall[idx1]*self.ibdall[idx2])
                    #print tt[idx1,idx2]
                    #print [idx1,idx2]
                    ##get the covariance matrix for SS, DD, SD
                    #get the correlations for different pairs
                    rik = self.corre([pair1[0],pair2[0]])       
                    ril = self.corre([pair1[0],pair2[1]])
                    rjk = self.corre([pair1[1],pair2[0]])
                    rjl = self.corre([pair1[1],pair2[1]])
                    ss[idx2,idx1]=ss[idx1,idx2] = 2*(rik+ril+rjk+rjl)**2
                    dd[idx2,idx1]=dd[idx1,idx2] = 2*(rik+rjl-ril-rjk)**2
                    sd[idx1,idx2] = 2*(rik+rjk-ril-rjl)**2
                    sd[idx2,idx1] = 2*(rik+ril-rjk-rjl)**2
            self.covTT = np.matrix(tt)
            self.covSS = np.matrix(ss)
            self.covDD = np.matrix(dd)
            self.covSD = np.matrix(sd)
            #print self.covSD

        def reg(self):
            ##Regression Step
            nfnum = len(self.nonfounder)
            pairnum = len(self.pairs)
            n=pairnum+nfnum
            #covY = np.zeros((n,n))
            self.Y= np.hstack([self.S,self.D[0,0:nfnum]]).T
            covss = self.covSS
            covsd = self.covSD[:,0:nfnum]
            covdd = self.covDD[0:nfnum,0:nfnum]
            ###########
            covY = np.vstack([np.hstack([covss,covsd]), np.hstack([covsd.T,covdd])])
            Yc = self.Y-np.hstack([self.ES,self.ED[0,0:nfnum]]).T
            Tc = (self.ibdall-self.ET).T
            diag1 = np.diag([2]*pairnum)
            diag2 = np.diag([-2]*pairnum)
            H = np.hstack([diag1,diag2[:,0:nfnum]])
            B = H*covY.getI()*Yc
            self.numerator = B.T*Tc
            self.denominator = B.T*self.covTT*B
            
        def getGT(self,v,foundergt):
            #assign nonfounders' genotypes based on inv and founder genotypes
                alleles=copy.copy(foundergt)
                gpid,gmid = sorted(set(self.parentsid))[0:2]
                for idx,mid in enumerate(self.nonfounder):   
                    #pick out affected nonfounders
			alleles[mid]=[None,None]
			#corresponding inheritance vector for this member
                        mv=v[2*idx:2*idx+2]                     
                        mp = self.fam_dict[mid]['parents']
			flag = 0                           #flag=0 => father important
                        if mp[0] in self.nonfounder:
                            pn = self.nonfounder.index(mp[0])
                        elif mp[1] in self.nonfounder:
                            flag = 1                       #flag=1 => mother important
                            mn = self.nonfounder.index(mp[1])
                        else:
                            flag = 2                       #flag=2 => member's parents are grandparents
			if flag == 2:                           #if parents are grandparents
                            alleles[mid][0] = foundergt[gpid][0] if mv[0] == '0' else \
                                    foundergt[gpid][1]
                            alleles[mid][1] = foundergt[gmid][0] if mv[1] == '0' else \
                                    foundergt[gmid][1]
			else:
                            if flag == 0:        #if father is the kid of grandparents
                                if mv[0] == '0':
                                    alleles[mid][0] = foundergt[gpid][0] if v[2*pn] == '0' \
                                            else foundergt[gpid][1]
                                elif mv[0] == '1':
                                    alleles[mid][0] = foundergt[gmid][0] if v[2*pn+1] == '0'\
                                            else foundergt[gmid][1]
                                alleles[mid][1] = foundergt[mp[1]][0] if mv[1]=='0' \
                                        else foundergt[mp[1]][1]
                            elif flag == 1:      #if mother is the kid of grandparents
                                if mv[1] == '0':
                                    alleles[mid][0] = foundergt[gpid][0] if v[2*mn] == '0' \
                                            else foundergt[gpid][1]
                                elif mv[1] == '1':
                                    alleles[mid][0] = foundergt[gmid][0] if v[2*mn+1] == '0'\
                                            else foundergt[gmid][1]
                                alleles[mid][1] = foundergt[mp[0]][0] if mv[0]=='0' \
                                        else foundergt[mp[0]][1]
			alleles[mid] = sorted(alleles[mid])
                        #sort the alleles since order is not important here
                return alleles

        def matrix_cal(self,alleles,ibdall,traits,conditional_prob):
	    #the whole matrix calculation step
            ##set the values for the matrix S, D and IBDm
            s = []
            es = []
            d = []
            ed = []
            tau = []
            for aid,bid in self.pairs:            
                #for each pair of relatives (include parent-offspring)
                trait1 = traits[aid]
                trait2 = traits[bid]
                s.append((trait1+trait2)**2)        #squared sum of traits
                d.append((trait1-trait2)**2)        #squared difference of traits
                r = self.corre([aid,bid])
                es.append(2*(1+r))
                ed.append(2*(1-r))
                tau.append(self.expect([aid,bid]))
            S = np.matrix(s)
            D = np.matrix(d)
            ES = np.matrix(es)
            ED = np.matrix(ed)
            ibdm = np.matrix(ibdall)
            ET = np.matrix(tau)
            ###########################
            postinv_prob,postibd=self.combine_inherit_vec(alleles=alleles,conditional_prob=conditional_prob)
            #print alleles
            #print postibd
            ##set the covariance matrix
            pairnum=len(self.pairs)
            nfnum = len(self.nonfounder)
            ss = np.zeros((pairnum,pairnum))
            dd = np.zeros((pairnum,pairnum))
            sd = np.zeros((pairnum,pairnum))
            tt = np.zeros((pairnum,pairnum))
            for idx1 in range(pairnum):
                for idx2 in range(idx1,pairnum):
                    ##get the covariance matrix for tau(i,j),tau(k,l)
                    #print pair1,pair2
                    pair1 = self.pairs[idx1]
                    pair2 = self.pairs[idx2]
                    prior,post = 0,0
                    prior_expect = 0
                    prior_expect=self.expect(pair1)*self.expect(pair2)
                    prior = self.prior[idx1][idx2]
                    if self.info == 1: #if posterior inv equals to priterior (marker not informative)
                        post = prior
                    else:
                        for inv in postinv_prob.keys():
                        #print inv
                        #print self.postibd[ivx][idx1],self.postibd[ivx][idx2]
                            post += postibd[inv][idx1]*postibd[inv][idx2]*postinv_prob[inv]
                    tt[idx2,idx1]=tt[idx1,idx2] = (prior-prior_expect)-(post-ibdall[idx1]*ibdall[idx2])
                    #print tt[idx1,idx2]
                    #print [idx1,idx2]
                    ##get the covariance matrix for SS, DD, SD
                    #get the correlations for different pairs
                    rik = self.corre([pair1[0],pair2[0]])       
                    ril = self.corre([pair1[0],pair2[1]])
                    rjk = self.corre([pair1[1],pair2[0]])
                    rjl = self.corre([pair1[1],pair2[1]])
                    ss[idx2,idx1]=ss[idx1,idx2] = 2*(rik+ril+rjk+rjl)**2
                    dd[idx2,idx1]=dd[idx1,idx2] = 2*(rik+rjl-ril-rjk)**2
                    sd[idx1,idx2] = 2*(rik+rjk-ril-rjl)**2
                    sd[idx2,idx1] = 2*(rik+ril-rjk-rjl)**2
            covTT = np.matrix(tt)
            covSS = np.matrix(ss)
            covDD = np.matrix(dd)
            covSD = np.matrix(sd)
            ##Regression Step
            n=pairnum+nfnum
            #covY = np.zeros((n,n))
            Y= np.hstack([S,D[0,0:nfnum]]).T
            covss = covSS
            covsd = covSD[:,0:nfnum]
            covdd = covDD[0:nfnum,0:nfnum]
            ###########
            covY = np.vstack([np.hstack([covss,covsd]), np.hstack([covsd.T,covdd])])
            Yc = Y-np.hstack([ES,ED[0,0:nfnum]]).T
            Tc = (ibdall-ET).T
            diag1 = np.diag([2]*pairnum)
            diag2 = np.diag([-2]*pairnum)
            H = np.hstack([diag1,diag2[:,0:nfnum]])
            B = H*covY.getI()*Yc
            numerator = B.T*Tc
            denominator = B.T*covTT*B
            return numerator[0,0], denominator[0,0]


        def clean(self):
            self.ibdall = []
            self.un=[]                                      
            self.cousin=[]
            self.sib=[]
            self.hsib=[]
            self.ibd_cousin=[]                  
            self.ibd_sib=[]
            self.ibd_hsib=[]
            self.ibd_un=[]

        def nulldist(self):
            print "calculating null dist for a single fam.."
            rep = 3000
            gpid,gmid = sorted(set(self.parentsid))[0:2]
            founderid = []
            outnd = []
            nfnum = len(self.nonfounder)
            traits = {}
            for iid in self.fam_dict.keys():
                traits[iid] = self.fam_dict[iid]['trait']
                if iid not in self.nonfounder:
                    founderid.append(iid)
            fnum = len(founderid)
            ###parallel processing
            if self.info != 1:
                #if foundergt informative
                outnd=Parallel(n_jobs=16)(delayed(null_generator)(founderid,self,traits) for i in xrange(rep))
            else:
                #if foundergt not informative
                (n,d) = null_generator(founderid,self,traits)
                outnd=[(n,d) for x in xrange(rep)]
            self.null_nd = outnd
            ###test result from matrix_cal###############
            #realallele = {}
            #realtrait = {}
            #for iid in self.fam_dict.keys():
            #    realallele[iid]=self.fam_dict[iid]['gt'][2*fam.mid:2*fam.mid+2]
            #    realtrait[iid] = self.fam_dict[iid]['trait']
            #n,d = self.matrix_cal(realallele,self.ibdall,realtrait)
        

        def execute(self,m):
	    #set the value for matrixs
            self.setmatrix()
            self.combine_inherit_vec(m) #calculate the overall posterior inv probability and corresponding IBD #time-consuming
            self.setcov()      #set covariance matrix #time-consuming
            self.reg()    #regression
            
def randomsample(null):
    rep = 20000
    cdist = []
    for trial in xrange(rep):
        numerator=0
        denominator=0
        Q = 0
        for dist in null:
            pick=random.randint(0,2999)
            n,d=dist[pick]
            numerator +=n
            denominator +=d
        Q=numerator/denominator
        cdist.append(Q)
    return sorted(cdist)


print args.path
famstruct = []
famstruct_null = {}
fam_null_nd = []             #[{'info':(n,d);'info2':(n,d)},{fam2}...]
traits=[]
for (dirpath,dirname,filename) in os.walk(args.path):
    if re.search(r'rep',dirpath):                           
        #find directories containing replicates
        fped = [f for f in filename if re.search(r'ped',f)]
        fdat = [d for d in filename if re.search(r'dat',d)]
        if len(fped) == 1 and len(fdat) == 1:
            pedfile = fped[0]
            markerfile = fdat[0]
            m = re.search(r'\/(rep\d+)\/',dirpath)
            rep = m.group(1)                   #get the replicate id
            chr_m = re.search(r'(chr\d+)',markerfile)
            chr_id = chr_m.group(1)                   #get the chromosome id
            print rep,chr_id
            cdist=[]    #combined distribution for each replicate
            pdirname = os.path.dirname(os.path.dirname(dirpath))
            #get the marker names
            markername = []
            markerfile = '{}/'.format(dirpath)+markerfile
            pedfile = '{}/'.format(dirpath)+pedfile
            cache_pdirname = '{}/cache/CACHE'.format(pdirname)
            mfile = file(markerfile,'r').read()
            mlines = mfile.split('\n')[1:-1]
            for element in mlines:
                markername.append(element.split()[-1])
            #output original result
            result = file("{}/original_result.txt".format(pdirname),'a')
            pv = file("{}/pvalue_HEextend.txt".format(pdirname),'a')
            rtext = '%s:\n'%rep
            #take in ped files
            f = file(pedfile,'r').read()
            lines = f.split('\n')
            lines.pop()
            family = []
            fam_num = 0
            marker_n=0
            for line in lines:
                tmp = []
                for i in line.split():                          #get every element from line
                    try:
                        tmp.append(int(i))
                    except ValueError:
                        try:
                            tmp.append(float(i))
                        except ValueError: #if things like '?' occurred
                            tmp.append(0)
                fid = tmp[0]                                   #family id
                marker_n = int((len(tmp)-6)/2)
                if fid >len(family):
                    family.append([])                          #family=[[[fam1member1][fam1member2]...],[fam2member1][]...]
                family[fid-1].append(tmp[1:])
            markers_to_analyze=[]
            if args.null:
                try:
                    m = markername.index(args.mname)
                    markers_to_analyze.append(m)
                except ValueError:
                    if args.snv:
                        continue
                    else:
                        markers_to_analyze=range(marker_n)
            else:
                markers_to_analyze=range(marker_n)
            if not args.snv:
                #get the freq file for CHP markers
                CHP_freq_byfam={}
                CHP_freqfile = '{}/'.format(cache_pdirname)+'{}.{}.freq'.format(rep,chr_id)
                try:
                    # read from freq file to get freq information for CHP alleles
		    # freq=[Freq(major allele), Freq(minor allele)...]
                    freqfile = file(CHP_freqfile,'r').read()        
                    freqlines = freqfile.split('\n')[:-1]
                    for element in freqlines:
                        CHP_freq=[]
                        tmp = element.split()
                        if tmp[1] == args.mname:
                            tmp_freqs=[float(x) for x in tmp[2:]]
                            sum_freq=sum(tmp_freqs)
                            if sum_freq == 1:
                                #when the sum of freqs is 1, stating that only one RV locus present here
                                CHP_freq=tmp_freqs
                            elif sum_freq < 1:
                                #when the sum is smaller than 1
                                #stating that there are more than one RV loci present in family
                                #adjust the freqs so that the sum is 1
                                CHP_freq=[freq/sum_freq for freq in tmp_freqs]
                            else:
                                #when the sum is larger than 1
                                #stating that there is a missing allele denoted as '0' 
                                #resulted from imputation step in SEQLinkage
                                #use the freqs of alleles other than '0'
                                CHP_freq=[freq/sum(tmp_freqs[:-1]) for freq in tmp_freqs[:-1]]
                            CHP_freq_byfam[int(tmp[0])]=CHP_freq
                    #print CHP_freq_byfam
                except:
                    pass
            fam_num = len(family)
            numerator = {}                          #marker:[nom1,nom2...]               
            denominator = {}
            asymp_pv = {}
            sum_null_nd = {}
            Q_all = {}
            for m in markers_to_analyze:
                print markername[m]
                sum_null_nd[markername[m]]=[]
                family_to_analyze = range(fam_num)
                for fid in range(fam_num):
                    print "processing family %d ..."%(fid+1)
                    try:
                        fam_null_nd[fid]
                    except IndexError:
                        fam_null_nd.append({})
                    fam = Family(m)
                    if args.snv:
                        tmp_maf=maf_by_marker[markername[m]]
                        fam.mp_freq=[tmp_maf,1-tmp_maf]
                    elif CHP_freq_byfam.has_key(fid+1):
                        fam.mp_freq=CHP_freq_byfam[fid+1]
                    fam.setdict(family[fid])  #set family data structure as a dictionary
                    fam.classify_affect(m) #classify affected individuals into relative pairs and calculate IBD
                    #set traits if initial trait values are 0
                    if fam.traitflag == 1:   
                        try:
                            traits[fid]
                        except IndexError:
                            traits.append(fam.setTrait(mod=0))
                        for iid in fam.fam_dict.keys():
                            #make the same family of all replicates share the same traits
                            fam.fam_dict[iid]['trait'] = traits[fid][iid]
                    if fam.err == 1:
                        family_to_analyze[fid]=None
                        fam.clean()
                        continue
                    else:
			#get the prior IBD product for each pair of affected relative pairs
			#fam.prior=sum(prob*pi(i,j)*pi(k,l))
                        if fam.famstruct not in famstruct:
                            #print "new structure"
                            famstruct.append(fam.famstruct)
                            fam.null_inv()
                            famstruct_null[famstruct.index(fam.famstruct)]=fam.prior
                        else:
                            idx = famstruct.index(fam.famstruct)
                            fam.prior = famstruct_null[idx]
			#get the null distribution of numerator and denominator via permutations 
                        if not fam_null_nd[fid].has_key(fam.info):
                            fam.nulldist()
                            fam_null_nd[fid][fam.info]=fam.null_nd
                        else:
                            #same family between replicates share same traits and \
                            #having same founder informativity also share the same null_nd
                            fam.null_nd = fam_null_nd[fid][fam.info]
                        sum_null_nd[markername[m]].append(fam.null_nd)
                        fam.execute(m)
                        try:
                            numerator[m].append(fam.numerator)
                            denominator[m].append(fam.denominator)
                        except KeyError:
                            numerator[m] = []
                            denominator[m] = []
                            numerator[m].append(fam.numerator)
                            denominator[m].append(fam.denominator)
                        #print fam.null_nd
                        fam.clean()
                Q = sum(numerator[m])/sum(denominator[m])
                t = Q**2*sum(denominator[m]) if Q>0 else 0
                p_asymp = stats.chi2.sf(t,1)/2
                Q_all[markername[m]]=Q
                asymp_pv[markername[m]]=p_asymp
            p_emp_min=1
            for ele in sorted(asymp_pv.items(), key=lambda x:x[1])[:5]:
                apv=ele[1]
                mname=ele[0]
                #if apv<0.5:
		#random sample from null distribution of (numerator,denominator) of families to get the null distribution of Q
                cdist = randomsample(sum_null_nd[mname])
                Q=Q_all[mname]
                p_emp=len([q for q in cdist if q >= Q])/len(cdist)
                rtext += '%s\t%.4f\t%.7f\t%.7f\n'%(mname,Q,apv,p_emp)
                p_emp_min = p_emp if p_emp<p_emp_min else p_emp_min
            result.write(rtext)
            pv.write('%s\t%.7f\n'%(rep,p_emp_min))
            pv.close()
            result.close()
