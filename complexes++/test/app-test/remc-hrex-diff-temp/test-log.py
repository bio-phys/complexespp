#!/usr/bin/env python
# Copyright (c) 2018 the complexes++ development team and contributors
# (see the file AUTHORS for the full list of names)
#
# This file is part of complexes++.
#
# complexes++ is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# complexes++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with complexes++.  If not, see <https://www.gnu.org/licenses/>
import re
import numpy as np
import pandas as pd
import sys
import os
from glob import glob

##########################################################

def parse_sweep(line):
    return int(line.split(':')[-1])


def parse_delta_energy(line):
    return float(line.split(':')[-1].split()[0])


def parse_attempted_exchange(line):
    return line.split(':')[-1].lstrip()
    #return tuple(np.int64(line.split(':')[-1].split('-')))


def parse_exchanges(line):
    line = line.split(':')[-1].lstrip()
    changes = dict(c.split('=') for c in line.split(','))
    exchanges = {}
    for k, v in changes.iteritems():
        if v == 'x':
            v = 1
        elif v == 'o':
            v = 0
        else:
            v = float('NaN')
        exchanges[k.lstrip().strip()] = v

    return exchanges


def parse_probabilities(line):
    line = line.split(':')[-1].lstrip()
    changes = dict(c.split('=') for c in line.split(','))
    probabilities = {}
    for k, v in changes.iteritems():
        try:
            number = float(v)
        except ValueError:
            number = float('NaN')
        probabilities[k.lstrip().strip()] = number
    return probabilities


converter = {'sweep': parse_sweep,
             'delta_energy': parse_delta_energy,
             'attempted_exchange': parse_attempted_exchange,
             'exchanges': parse_exchanges,
             'probabilities': parse_probabilities}

def parse_exchange_log(fname):
    with open(fname) as f:
        log = f.readlines()
    attempt_lines = np.where(['Replica_Attempt' in line for line in log])[0]

    attempts = []
    for attempt in attempt_lines:
        d = {}
        # attemp record is maximal 5 lines, add 1 for header line
        for line in log[attempt + 1: attempt + 6]:
            for key, func in converter.iteritems():
                if key in line:
                    d[key] = func(line.replace('\n', ''))
            # we reached end of record
            if 'LOG' in line:
                break

        attempts.append(d)
    return attempts

##########################################################
def parse_get_transition_matrix(fname, matdim):
    with open(fname) as f:
        log = f.readlines()
    attempt_lines = np.where(['empirical_transition_matrix' in line for line in log])[0]

    if len(attempt_lines) < 1 :
        raise ValueError('There are {} transition matrices found in the file {}'.format(len(attempt_lines), fname))

    textmat = "".join([line for line in log[attempt_lines[0] + 1: attempt_lines[0] + 1 + matdim]])

    # use the first one
    return np.matrix(textmat).reshape(matdim, matdim)

##########################################################
def results_are_similar(res1, res2):
    # simply compare the resulting string
    if str(res1) != str(res2) :
        return False    

    return True

##########################################################
verboseOutput = True if os.getenv('VERBOSE', "0") == "1" else False
if verboseOutput :          
    print('Verbose is turned on' )

all_ok=True
# compare main logs
remc_main = parse_exchange_log('remc.log.out') 
mat_remc_main = parse_get_transition_matrix('remc.log.out', 4)
hrex_main = parse_exchange_log('hrex.log.out')        
mat_hrex_main = parse_get_transition_matrix('hrex.log.out', 4)

if not results_are_similar(remc_main,hrex_main) :
    if verboseOutput:
        print('Main log are not equal')
    all_ok=False

if not np.allclose(mat_remc_main, mat_hrex_main, rtol=0, atol=.05)  :
    if verboseOutput:
        print('Main log matrices are not equal')
    all_ok=False

# compare each dir log
for i in range(0,3):
    hrexlog = parse_exchange_log('test' + str(i) + '/complexes-hrex.log')
    mat_hrex = parse_get_transition_matrix('test' + str(i) + '/complexes-hrex.log', 4)

    remclog = parse_exchange_log('test' + str(i) + '/complexes-remc.log')    
    mat_remc = parse_get_transition_matrix('test' + str(i) + '/complexes-remc.log', 4)

    if not results_are_similar(remclog,hrexlog) :
        if verboseOutput:
            print('test' + str(i) + '/ is not correct')
        all_ok=False

    if not np.allclose(mat_hrex, mat_hrex_main, rtol=0, atol=.05) :
        if verboseOutput:
            print('test' + str(i) + '/complexes-hrex.log matrix is not correct')
        all_ok=False

    if not np.allclose(mat_remc_main, mat_remc, rtol=0, atol=.05) :
        if verboseOutput:
            print('test' + str(i) + '/complexes-remc.log matrix is not correct')
        all_ok=False

# ensure transition matrix coherency
nb_exchanges = np.zeros(shape=(4,4))
nb_total_attemps = 0

for attempt in remc_main:
    for exchanges in attempt['exchanges'] :
        idx1, idx2 = exchanges.split('-')
        result = attempt['exchanges'][exchanges]
        if result == 1 :
            nb_exchanges[idx1, idx2] += 1
            nb_exchanges[idx2, idx1] += 1
        nb_total_attemps += 1

exchange_proba = np.zeros(shape=(4,4))
for row in range(0, 4):         
    for col in range(0, 4):       
        if row != col:  
            coef = nb_exchanges[row, col] / nb_total_attemps
            exchange_proba[row, col] = coef
            exchange_proba[row, row] -= coef
            exchange_proba[col, col] -= coef
        else :
            exchange_proba[col, col] += 1
            
if not np.allclose(mat_hrex, mat_hrex_main, rtol=0, atol=.05) :
        if verboseOutput:
            print('transition matrix is wrong')
            print('expected transition matrix')
            print(exchange_proba)
        all_ok=False


sys.exit(not all_ok)
