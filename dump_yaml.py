#!/usr/bin/env python3

sys.path.append( './malt' )
from malt.algorithm.topology import Fragmentation

import pickle
import sys
import yaml
import numpy as np
from configparser import ConfigParser
from fireworks import LaunchPad
from glom import glom
from tqdm import tqdm

def create_launchpad(db_config_file):
    """use to create a FW launchpad using mongodb creds from file"""
    config = ConfigParser()
    config.read(db_config_file)
    db = config['db']

    lpad = LaunchPad(
        host=db['host'],
        port=int(db['port']),
        name=db['name'],
        username=db['username'],
        password=db['password'])

    return lpad

def get_family(fw):
    parent = glom(fw, 'launches.0.action.stored_data.parent')
    children = glom(fw, 'launches.0.action.stored_data.protonated_children')

    return parent, children


def get_record(parent):
    record = {}
    try:
        record['symbol'] = [str(elem[0]) for elem in parent['atom_list']]
        record['position'] = parent['coords'][-1]
    except:
        record = None
       
    return record

def get_children_results(children):       
    children_data = []

    for child in children:
        record = {}
        try:
            symbols = [str(elem[0]) for elem in child['atom_list']]
            record['bond_break_prediction'] = child['likely_to_break']
            record['site'] = [str(child['element']), child['atom_index']]

            last = child['coords'][-1]
            lap = frag.get_laplacian(np.asarray(last), symbols, longrange)
            eigvals = np.linalg.eigvalsh(lap)
            if len([eig for eig in eigvals if eig < 10**-10]) > 1:
                record['fragmented'] = True
            else:
                record['fragmented'] = False
        except:
            record = None

        if record:
            children_data.append(record)

    return children_data

if __name__ == '__main__':
	LOCAL_DB_CONFIG = '/home/bkrull/.fireworks/local_db.ini'
	lpad = create_launchpad(LOCAL_DB_CONFIG)

	frag = Fragmentation()
	longrange = frag.weights_longrange


	base = '/scratch/users/bkrull/yaml/'

	completed = lpad.get_fw_ids({'state': 'COMPLETED'})
		                    
	for id in tqdm(completed, total=len(completed)):
	    data = lpad.get_fw_dict_by_id(id)
	    name = data['name'].split('/')[-1]
	    
	    parent, children = get_family(data)
	    record = get_record(parent)        
	    record['protonation_tests'] = get_children_results(children)
		
	    if record['protonation_tests']:
		with open(base+'{}.yaml'.format(name), 'w') as f:
		    yaml.dump(record, f)
