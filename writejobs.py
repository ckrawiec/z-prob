run_script = '/home/ckrawiec/git/z-prob/runzprob.py'
num_threads = 8
mem_per_thread = '3900M'
n_tab_groups = 11
n_jobs_per_tab = 1250
n_targets_per_job = 10000 #None=do all
config_dir = '/home/ckrawiec/git/z-prob-config/balrog/'

name = 'svbalrogtab{}_y1dfullcosmos_flagBRmasked_auto_griz_full_gauss_z4_6bins+GalaxSV3903sig'

galaxy_template_file = '/home/ckrawiec/DES/data/y1a1_gold_dfull_cosmos_flagBRmasked.fits'

star_template_file = '/home/ckrawiec/DES/data/Galaxia_SV390deg_desflux_radec.fits'

target_file = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM_zp_corr_fluxes.fits'
#y1a1_gold_d10_dfull_cosmos_matched_d10flagBRmasked_dfullflagBRmasked.fits'

config_dict = {}

for i in range(1, n_tab_groups+1):
    tab_num = str(i).zfill(2)

    if n_targets_per_job:
        start = 0
        end = n_targets_per_job

    for j in range(n_jobs_per_tab):
        if n_targets_per_job:
            inds = '{}-{}'.format(start, end)
            config_dict['start'] = start
            config_dict['end'] = end
        else:
            inds = 'all'
            config_dict['start'] = ''
            config_dict['end'] = ''

        config = '{}{}_{}.config'.format(config_dir, name, inds).format(tab_num)
        oe = '/home/ckrawiec/jobscripts/output/zprob_{}_{}.oe'.format(name, inds).format(tab_num)
        job = '/home/ckrawiec/jobscripts/zprob_{}_{}.sh'.format(name, inds).format(tab_num)

        config_dict['tab_num'] = tab_num

        config_dict['num_threads'] = num_threads
        config_dict['galaxy_template_file'] = galaxy_template_file
        config_dict['star_template_file'] = star_template_file
        config_dict['target_file'] = target_file.format(tab_num)
        config_dict['output_file'] = '/home/ckrawiec/DES/magnification/lbgselect/zproboutput/{}_{}.fits'.format(name, inds).format(tab_num)

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
target_id_column = BALROG_INDEX
target_data_column = FLUX_AUTO_{}
target_error_column = FLUXERR_AUTO_{}

#template column names, case sensitive
galaxy_template_id_column = COADD_OBJECTS_ID
galaxy_template_data_column = FLUX_AUTO_{}
#sq. deg.
galaxy_area = 2.
redshift_column = zminchi2

star_template_id_column = satid
star_template_data_column = flux_des_{}
#sq. deg.
star_area = 390.
#specify indices to use slice of target data
#if not specified, all are used
\n""" % config_dict)

        if n_targets_per_job:
            f.write("""
target_start_index = %(start)s
target_end_index = %(end)s
""" % config_dict)
            f.close()
        else:
            f.close()

        f = open(job,'w')
        f.write('#$ -pe omp {}\n'.format(num_threads))
        f.write('#$ -l h_vmem={}\n'.format(mem_per_thread))
        f.write('#$ -o {}\n'.format(oe))
        f.write('#$ -e {}\n'.format(oe))
        f.write('#$ -N j{}_{}'.format(tab_num, j))
        f.write('\n\n')
        f.write('echo python {} {}\n'.format(run_script, config))
        f.write('python {} {}'.format(run_script, config))
        f.close()

        if n_targets_per_job:
            start += n_targets_per_job
            end += n_targets_per_job
