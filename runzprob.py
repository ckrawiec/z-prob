usage = """
Calculate the Bayesian probability that
fluxes of input galaxies belong to 
certain redshift ranges using photo-z's.

usage: python runzprob.py <config file>

outputs a binary fits table and will overwrite existing file
"""

import numpy as np
import zprobability as zprob
import ConfigParser
import json
import time
import sys
import os
from astropy.io import fits

def parseconfig(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    params['template_file'] = config.get('I/O','template_file')
    params['target_file'] = config.get('I/O','target_file')
    params['output_file'] = config.get('I/O','output_file')
    if os.path.exists(params['output_file']):
        sys.exit("Error: {}\nOutput file already exists. Delete before running again.\n    Exiting.".format(params['output_file']))
    
    params['filters'] = config.get('parameters','filters')
    params['redshift_ranges'] = json.loads(config.get('parameters','redshift_ranges'))

    params['num_threads'] = config.getint('parameters','num_threads')
    params['integration'] = config.get('parameters','integration')
    global query_radius
    query_radius = None
    if params['integration'] == 'tree':
        query_radius = config.getfloat('parameters','query_radius')
        
    params['template_id_column'] = config.get('data','template_id_column')
    params['template_data_column'] = config.get('data','template_data_column')
    params['redshift_column'] = config.get('data','redshift_column')
    
    params['target_id_column'] = config.get('data','target_id_column')
    params['target_data_column'] = config.get('data','target_data_column')
    params['target_error_column'] = config.get('data','target_error_column')

    #if start/end indices not specified, use all targets
    if config.has_option('data','target_start_index'):
        params['target_start_index'] = config.getint('data','target_start_index')
    else:
        params['target_start_index'] = 0

    if config.has_option('data','target_end_index'):
        params['target_end_index'] = config.getint('data','target_end_index')
    else:
        params['target_end_index'] = None
            
    params['debug'] = config.getboolean('debugging','debug')
    if params['debug'] == True:
        params['N_debug'] = config.getint('debugging','N_debug')

    return params

def printtime():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now

def writetofile(params, Pdict, targets):
    col_defs = [fits.Column(name=params['target_id_column'], format='K', array=targets.ids)]
        
    P_norm = np.zeros(len(targets.data))
    for k in Pdict.keys():
        P_norm += Pdict[k]
    for k in Pdict.keys():
        col_defs.append(fits.Column(name='P'+k, format='D', array=Pdict[k]/P_norm))

    col_defs.append(fits.Column(name='PNORM', format='D', array=P_norm))

    pri_hdr = fits.Header()
    tb_hdr = fits.Header()
    tb_hdr['COMMENT'] = "Bayesian redshift probabilities for data in {} using photo-zs of templates from {}. \
                         Data vectors were comprised of {} for bands in {} with matching errors. \
                         Columns reported here are \'P[zmin, zmax]\'".format(params['target_file'],
                                                                             params['template_file'],
                                                                             params['target_data_column'],
                                                                             params['filters'])
                                                                        
    if params['integration']=='tree':
        tb_hdr['RADIUS'] = str(query_radius)

    pri_hdu = fits.PrimaryHDU(header=pri_hdr)
    tb_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_defs), nrows=len(targets.data), header=tb_hdr)
    hdu_list = fits.HDUList([pri_hdu, tb_hdu])
    hdu_list.writeto(params['output_file'], clobber=True)

def main(args):
    params = parseconfig(args[1])
    
    printtime()

    setup_start = time.time()

    templates = zprob.Templates(params['template_file'],
                                params['template_id_column'],
                                params['template_data_column'],
                                params['redshift_column'],
                                params['filters'])

    targets = zprob.Targets(params['target_file'],
                            params['target_id_column'],
                            params['target_data_column'],
                            params['target_error_column'],
                            params['target_start_index'],
                            params['target_end_index'],
                            params['filters'])
    
    setup_end = time.time()
    print "Loaded data in {} s".format(setup_end-setup_start)
    
    if params['debug']==True:
        targets.debug(templates,
                      params['N_debug'],
                      params['redshift_ranges'],
                      params['num_threads'],
                      params['output_file'])
        printtime()
        sys.exit()

    print "Working on {} galaxies (table indices {}-{}) ...".format(len(targets.data),
                                                                    params['target_start_index'],
                                                                    params['target_end_index'])

    #get Bayesian probabilities using templates and redshift ranges    
    P_dict = targets.calcProbabilities(templates,
                                       params['redshift_ranges'],
                                       params['num_threads'],
                                       params['integration'],
                                       query_radius)
        
    #write results to fits file
    writetofile(params, P_dict, targets)
   
    printtime()
        
if __name__=="__main__":
    main(sys.argv)
