import numpy as np
import os
import logging
import argparse
import random
import math
from classes import *

MAX = 9999.0
#只考虑最大峰值，并默认为分子离子峰
#先根据分子离子峰结果找出所有候选分子式，再模拟出同位素分布，与谱图进行匹配
#利用经验规则筛选：原子总数的两倍减去一必须大于或等于原子的价数之和（化合价之和小于等于2*（atom-1））
def reading_spectrum(file_name):
    MS_data_list= list()
    with open(file_name,'r',encoding='utf-8') as file:
        MS_data_position= list()
        MS_data_height = list()
        for line in file:
            if line=='\n':
                continue
            data = line.replace('\n','').split('\t')
            if not data[0][0].isdigit() and len(MS_data_position)!=0:
                MS_data_list.append((MS_data_position,MS_data_height))
                MS_data_position= list()
                MS_data_height = list()
            if data[0][0].isdigit():
                MS_data_position.append(float(data[0]))
                MS_data_height.append(float(data[2]))
        MS_data_list.append((MS_data_position,MS_data_height))
    return MS_data_list

def reading_dir(dir_file):
    with open(dir_file,'r',encoding='utf-8') as file:
        number = int(file.readline().replace('\n',''))
        atomic = (file.readline().strip().replace('\n','').split(' '))
        atomic_weight = list(map(lambda x:float(x),(file.readline().strip().replace('\n','').split(' '))))
        isotopic_pattern = list()
        for i in range(number):
            pattern = list(map(lambda x:float(x)/100.0,(file.readline().strip().replace('\n','').split(' '))))
            isotopic_pattern.append(pattern)
        atomic_valence = list(map(lambda x:int(x),(file.readline().strip().replace('\n','').split(' '))))
        return Dir(number,atomic,atomic_weight,isotopic_pattern,atomic_valence)


def formula_decomposition(molecular_mass,mass,dir,candidate_list,err,index,candidate):
    if mass < 0 or index==-1:
        
        mass_sum = sum([candidate[i]*dir.atomic_weight[i] for i in range(dir.number)])
        atomic_sum = sum(candidate)
        valence_sum = sum([candidate[i]*dir.atomic_valence[i] for i in range(dir.number)])
        if (abs((mass_sum-molecular_mass)/molecular_mass) <= err) and valence_sum>=2*atomic_sum-1 and valence_sum%2==0:
            candidate_list.append(candidate.copy())
        return
    max_number = round(mass / dir.atomic_weight[index])
    for i in range(max_number+1):
        candidate[index] = i
        formula_decomposition(molecular_mass,mass-dir.atomic_weight[index]*i,dir,candidate_list,err,index-1,candidate)
    candidate[index] = 0
    
def formula_selection(candidate_formulas,spectrum,spectrum_index,dir,topk):
    #由于设备精度限制，只考虑A+1,...,A+4的匹配程度
    intensity = np.zeros(4)
    score = np.zeros(len(candidate_formulas))
    
    for i in range(1,5):
        if spectrum[0][spectrum_index+i]-spectrum[0][spectrum_index]<=i+0.005 and spectrum[0][spectrum_index+i]-spectrum[0][spectrum_index]>=i-0.005:
            intensity[i-1] = spectrum[1][spectrum_index+i] / spectrum[1][spectrum_index]
    for j in range(len(candidate_formulas)):
        simulation_intensity = np.zeros(4)
        for k in range(1,5):
            possibility = 1
            index = 0
            isotopic_simulation(candidate_formulas[j],simulation_intensity,k,k,index,possibility,dir)
        score[j] = rmse(intensity,simulation_intensity)
    #    print(f'simulation_intensity:{simulation_intensity}')
    
    # print(f'intensity:{intensity}')
    
    combined = [(candidate_formulas[i],score[i]) for i in range(len(score))]
    combined.sort(key=lambda x:x[1])
    candidate_formulas,score = map(list,zip(*combined))
    
    for i in range(min(topk,len(score))):
        for j in range(dir.number):
            if candidate_formulas[i][j]!=0:
                logging.info(f'{dir.atomic[j]}{candidate_formulas[i][j]}')
        logging.info(score[i])
        logging.info('\n')

def rmse(original_intensity,simulation_intensity):
    score = 0
    for i in range(len(original_intensity)):
        score = score + (original_intensity[i]-simulation_intensity[i])**2
    return np.sqrt(score/4)
            
def isotopic_simulation(formula,simulation_intensity,intensity_index,k,index,possibility,dir):
    if k < 0 or (index==len(formula) and k>0):
        return
    if k==0:
        simulation_intensity[intensity_index-1] = simulation_intensity[intensity_index-1] + possibility
        return
    if len(dir.isotopic_pattern[index])>0:
        for j in range(1,len(dir.isotopic_pattern[index])):
            for x in range(int(formula[index])+1):
                factor = (dir.isotopic_pattern[index][j]**x)*math.factorial(int(formula[index]))/ (math.factorial(x)*math.factorial(int(formula[index]-x)))                                                                           
                possibility = possibility*factor
                
                isotopic_simulation(formula,simulation_intensity,intensity_index,k-j*x,index+1,possibility,dir)
                possibility = possibility / factor
    else:
        isotopic_simulation(formula,simulation_intensity,intensity_index,k,index+1,possibility,dir)
                
                
  
def solve_molecular_formula(spectrum,dir,topk):
    molecular_mass = 0
    spectrum_index = 0
    for i in range(len(spectrum[0])-1,1,-1):
        position = spectrum[0][i]
        height = spectrum[1][i]
        if height > 1 and spectrum[1][i-1]<height:
            molecular_mass=position
            spectrum_index = i
            break
    candidate_list = list()
    candidate = np.zeros(dir.number)
    
    formula_decomposition(molecular_mass,molecular_mass,dir,candidate_list,args.error,index=dir.number-1,candidate=candidate)
#    for i in range(len(candidate_list)):
#        logging.info(candidate_list[i])
    candidate = formula_selection(candidate_list,spectrum,spectrum_index,dir,topk)
    
    
    
    
                     
    
                
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='MassSpectrom with argsparse')
    parser.add_argument('--spectrum','-s',type=str,default='MS_spectrum.txt',help='Relative position of spectral file')
    parser.add_argument('--dir','-d',type=str,default='dir.txt',help='Relative position of dir file')
    parser.add_argument('--output','-o',type=str,default='output.txt',help='Relative position of output file')
    parser.add_argument('--topk','-k',type=int,default=5,help='Number of candidate results to output for each file')
    parser.add_argument('--error',type=float,default=0.005,help='Maximal error for mass')
    args = parser.parse_args()
    
    spectrum_file = args.spectrum
    dir = args.dir
    spectrum = reading_spectrum(spectrum_file)
    dir = reading_dir(dir)
    logging.basicConfig(filename=args.output, level=logging.INFO, format='%(message)s',filemode='w')
    for i in range(len(spectrum)):
        logging.info(f'Number:{i}\n')
        solve_molecular_formula(spectrum[i],dir,args.topk)
