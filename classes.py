import numpy as np
import os
import logging
import argparse

def max_H(dir,formula):
        max = 0
        for key in formula.keys():
            if key not in ('H','A2_residue','A1_residue','A_residue'):
                max = max+(dir.atomic_bond[dir.symbol2index[key]]-2)*formula[key]
        return max+2
    
class Dir():
    def __init__(self,number,atomic,atomic_weight,isotopic_pattern,atomic_valence):
        self.number = number
        self.atomic = atomic
        self.atomic_weight = atomic_weight
        self.isotopic_pattern = isotopic_pattern
        self.atomic_valence = atomic_valence
        self.symbol2index = dict()
        self.initialize_dict()
        
    def initialize_dict(self):
        for i in range(len(self.atomic)):
            self.symbol2index[self.atomic[i]] = i
            

class Peak():
    def __init__(self,A,A1_rel,A2_rel):
        self.A=round(A)
        self.A1_rel = A1_rel
        self.A2_rel = A2_rel
        self.A2_formula_list = list()
        self.A1_formula_list = list()
        self.A_formula_list = list()
        self.formula = dict()
        self.error_list = list()
        self.combination = dict()
        self.normalized_error_list = list()
        self.normalized_formula_list = list()
    
    def symbol2index(self,symbol,dir):
        for i in range(len(dir.atomic)):
            if dir.atomic[i]==symbol:
                return i
    
    def A2_search(self,dir,topk,formula,residue,key_index):
        if dir.A2_atomic[key_index]==0 or residue<0:
            formula['A2_residue'] = residue
            A1_residue = self.A1_rel
            for key in formula.keys():
                if key not in ('A1_residue','A2_residue','A_residue'):
                    A1_residue = A1_residue-formula[key]*dir.A1_atomic[dir.symbol2index[key]]
            formula['A1_residue'] = A1_residue
            self.A2_formula_list.append(formula.copy())
            return 
        max_number = min(max(int(residue / dir.A2_atomic[key_index]) + 1,2), int(self.A / dir.atomic_weight[key_index]) + 1 )
        for i in range(max_number):
            formula[dir.atomic[key_index]] = i
            self.A2_search(dir,topk,formula,residue-i*dir.A2_atomic[key_index],key_index+1)
    
    def A1_search(self,dir,topk,formula,residue,key_index):
        if (dir.A2_atomic[key_index]!=0 or dir.atomic[key_index]=='H'):
            self.A1_search(dir,topk,formula,residue,key_index+1)
            return
        if dir.A1_atomic[key_index]==0:
            formula['A1_residue'] = residue
            A_residue = self.A
            for key in formula.keys():
                if key not in ('A1_residue','A2_residue','A_residue'):
                    A_residue = A_residue-formula[key]*dir.atomic_weight[dir.symbol2index[key]]
            formula['A_residue'] = A_residue
            self.A1_formula_list.append(formula.copy())
            return 
        max_number = min(max(int(residue / dir.A1_atomic[key_index]) + 1,2), int(self.A / dir.atomic_weight[key_index]) + 1 )
    #    print(f'max_number{dir.atomic[key_index]}={max_number}')
        for i in range(max_number):
            formula[dir.atomic[key_index]] = i
            self.A1_search(dir,topk,formula.copy(),residue-i*dir.A1_atomic[key_index],key_index+1)
    
    def A_search(self,dir,topk,formula,residue,key_index):
        if dir.atomic[key_index]=='H':
            H = max_H(dir,formula)
        #    logging.info(f'H={H}')
        #    logging.info(f'residue={residue}')
            if residue<-0.5:
                formula['H']=0
                formula['A_residue'] = residue
            elif residue>H:
                formula['H'] = H
                formula['A_residue'] = residue-H
            else:
                formula['H'] = round(residue)
                formula['A_residue'] = residue-formula['H']
            if formula not in self.A_formula_list:
            #    logging.info(formula)
                self.A_formula_list.append(formula.copy())
            return
        if dir.A2_atomic[key_index] !=0 or dir.A1_atomic[key_index]!=0 or(residue<0 and dir.atomic[key_index]!='H'):
            self.A_search(dir,topk,formula,residue,key_index+1)
            return
        max_number = int(residue/dir.atomic_weight[key_index]) +1
        for i in range(max_number):
            formula[dir.atomic[key_index]] = i
            self.A_search(dir,topk,formula.copy(),residue-i*dir.atomic_weight[key_index],key_index+1)
        
    def error_calculate(self,dir):
        eps = 1e-4
        for formula in self.A_formula_list:
            formula['A2_residue'] = formula['A2_residue'] - formula['H'] * dir.A2_atomic[dir.symbol2index['H']]
            error_A = (100*formula['A_residue'] / (self.A+eps))  ** 2
            error_A1 = (formula['A1_residue'] / (self.A1_rel+eps)) ** 2
            error_A2 = (formula['A2_residue'] / (self.A2_rel+eps)) ** 2
            error = np.average((error_A,error_A1,error_A2))
            self.error_list.append(error)
        combined = list(zip(self.error_list,self.A_formula_list))
        combined.sort(key=lambda x:x[0])
        self.error_list,self.A_formula_list = map(list,zip(*combined))
        return
    
    def formula_search(self,dir,topk):
        formula = dict()
        for atomic in dir.atomic:
            formula[atomic]=0
            
        dir.A2_sorted()
        self.A2_search(dir,topk,formula,self.A2_rel,0)
        
        dir.A1_sorted()
        for i in range(len(self.A2_formula_list)):
            formula = self.A2_formula_list[i].copy()
            residue = formula['A1_residue']
        #    print(self.A2_formula_list[i])
            self.A1_search(dir,topk,formula,residue,0)
        
        dir.A_sorted()
        for i in range(len(self.A1_formula_list)):
            formula = self.A1_formula_list[i].copy()
            residue = formula['A_residue']
            #logging.info(self.A1_formula_list[i])
            self.A_search(dir,topk,formula,residue,0)
            
        self.error_calculate(dir)
        
        return

    def topk_output(self,topk):
        logging.info(f'Peak{self.A}:')
        for i in range(min(topk,len(self.A_formula_list))):
            logging.info(self.A_formula_list[i])
            for (key,value) in self.A_formula_list[i].items():
                if key not in ('A1_residue','A2_residue','A_residue') and value>0:
                    logging.info(f'{key}{self.A_formula_list[i][key]}'.rstrip())
            logging.info(f' error={self.error_list[i]}')
            formula = self.A_formula_list[i]
            eps = 1e-3
        #    error_A = (10*formula['A_residue'] / (self.A+eps))  ** 2
        #    error_A1 = (formula['A1_residue'] / (self.A1_rel+eps)) ** 2
        #    error_A2 = (formula['A2_residue'] / (self.A2_rel+eps)) ** 2
        #    logging.info(f' error_A={error_A}')
        #    logging.info(f' error_A1={error_A1}')
        #    logging.info(f' error_A2={error_A2}')
        return
    
    def combined(self,topk):
        number = min(len(self.A_formula_list),topk)
        self.normalized_error_list = [1 / x for x in self.error_list]
        self.normalized_error_list = (self.normalized_error_list[:number]) / (sum(self.normalized_error_list[:number]))
    #    print(self.normalized_error_list)
        for i in range(number):
            for (key,value) in self.A_formula_list[i].items():
                if key in self.combination.keys():
                    self.combination[key] = self.combination[key]+ self.normalized_error_list[i]*value
                else:
                    self.combination[key] = self.normalized_error_list[i]*value
        return
            