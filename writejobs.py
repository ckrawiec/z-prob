run_script = '/home/ckrawiec/git/z-prob/runzprob.py'
num_threads = 8
mem_per_thread = '4G'
n_groups = 1

name = 'y1a1_gold_cosmosd04_cosmosdfull_matched_bothflagBRmasked_auto_griz_full_gauss_z4_6bins+stars'

galaxy_template_file = '/home/ckrawiec/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBRmasked_dfullflagBRmasked.fits'

star_template_file = '/home/ckrawiec/DES/data/y3_gold_v1.0_stars_flux_000001.fits'

target_file = '/home/ckrawiec/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBRmasked_dfullflagBRmasked.fits'

config_dict = {}

for i in range(1, n_groups+1):

    tab_num = str(i).zfill(2)
    config = '/home/ckrawiec/git/z-prob/config/{}.config'.format(name).format(tab_num)
    oe = '/home/ckrawiec/jobscripts/output/zprob_{}.oe'.format(name).format(tab_num)
    job = '/home/ckrawiec/jobscripts/zprob_{}.sh'.format(name).format(tab_num)

    config_dict['tab_num'] = tab_num
    config_dict['num_threads'] = num_threads
    config_dict['galaxy_template_file'] = galaxy_template_file
    config_dict['star_template_file'] = star_template_file
    config_dict['target_file'] = target_file.format(tab_num)
    config_dict['output_file'] = '/home/ckrawiec/DES/magnification/lbgselect/zproboutput/{}.fits'.format(name).format(tab_num)

    f = open(config, 'w')
    f.write("""#config parameters

[I/O]
#file names
galaxy_template_file = %(galaxy_template_file)s
star_template_file = %(star_template_file)s
target_file = %(target_file)s
output_file = %(output_file)s

[debugging]
#debug parameters
debug = False
N_debug = 500

[parameters]
#bands to use, same case as column names require
filters = GRIZ

#redshift groups for which to calculate probability
#of target membership
redshift_ranges = [[0.001, 0.5], [0.5, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 9.9]]

#number of threads for multiprocessing
num_threads = %(num_threads)d
#'tree' or 'full' integration over templates
integration = full
#query_radius = 5.

[data]
#target column names, case sensitive
target_id_column = COADD_OBJECTS_ID_d04
target_data_column = FLUX_AUTO_{}_d04
target_error_column = FLUXERR_AUTO_{}_d04

#template column names, case sensitive
galaxy_template_id_column = COADD_OBJECTS_ID_dfull
galaxy_template_data_column = FLUX_AUTO_{}_dfull
redshift_column = zminchi2_dfull

star_template_id_column = COADD_OBJECT_ID
star_template_data_column = FLUX_AUTO_{}

#specify indices to use slice of target data
#if not specified, all are used
target_start_index = 0
target_end_index = 1000
""" % config_dict)
    f.close()

    f = open(job,'w')
    f.write('#$ -pe omp {}\n'.format(num_threads))
    f.write('#$ -l h_vmem={}\n'.format(mem_per_thread))
    f.write('#$ -o {}\n'.format(oe))
    f.write('#$ -e {}\n'.format(oe))
    f.write('#$ -N job'+tab_num)
    f.write('\n\n')
    f.write('echo python {} {}\n'.format(run_script, config))
    f.write('python {} {}'.format(run_script, config))
    f.close()


