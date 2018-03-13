import itertools
import time
import random
import sys
import numpy as np
from scipy.cluster.vq import kmeans, vq
from multiprocessing import Pool
from astropy.io import fits
from scipy.spatial import ckdtree


def likelihoods(val, err, truevals):
    """
    returns gaussian likelihood from each trueval given val & err
    """
    covI = 1./err**2.
    A = 1./np.sqrt( (2.*np.pi)**len(val) * np.prod(err**2.) )
    
    diff = val-truevals
    B = -0.5 * np.sum(diff**2. * covI, axis=1)
    return A * np.exp(B)
                
def P(vals, errs, truevals, nchunks=500, ntruechunks=1000):
    """
    sum the gaussian likelihoods L(vals|truevals) over truevals using the errs on vals
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """

    out = np.array([])

    #break data into chunks
    chunks = itertools.izip([vals[i:i+nchunks] for i in xrange(0, len(vals), nchunks)], 
                            [errs[i:i+nchunks] for i in xrange(0, len(vals), nchunks)])
    
    for chunk, errchunk in chunks:
        trueout = np.zeros(len(chunk))
    
        covIs = 1./errchunk**2.
        A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errchunk**2., axis=1 ) )

        #sum over all templates in chunks
        truechunks = (truevals[i:i+ntruechunks] for i in xrange(0, len(truevals), ntruechunks))
        for truechunk in truechunks:
            diff = chunk[:,np.newaxis,:]-truechunk[np.newaxis,:,:]

            B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)
            C = A[:,np.newaxis] * np.exp(B)
            
            trueout += np.sum(C, axis=1)
        
        #add target chunk to end of list
        out = np.concatenate((out, trueout))
    
    return out

def Ptree(vals, errs, truevals, treeunit=None, query_radius=3., thresh=1000):
    out = []
    lencheck = 0

    #build truth tree 
    start_tree = time.time()
    truetree = ckdtree.cKDTree(truevals/treeunit, balanced_tree=False, compact_nodes=False)
    end_tree = time.time()
    sys.stderr.write('tree created in {}s\n'.format(end_tree-start_tree))

    for val, err in zip(vals, errs):
        #get nearest neighbors within query_radius to target
        inear = truetree.query_ball_point(val/treeunit, r=query_radius)
        factor = 1.
        if len(inear)>thresh:
            lencheck+=1
            inear = random.sample(inear, thresh)
            factor = len(inear)/float(thresh)

        #data of nearest neighbors
        truearr = truetree.data[inear] * treeunit
        
        start_Ls = time.time()
        Ls = likelihoods(val, err, truearr)
        end_Ls = time.time()

        #sum likelihoods
        out.append(np.sum(Ls) * factor)

    if lencheck>0:
        sys.stderr.write("Warning: ball query returned >1m points for {} targets ({}%).\n".format(lencheck, 
                                                                                                  float(lencheck)/len(vals) * 100.))
    return np.array(out)

def Pwrapper(args):
    if args[0]=='tree':
        return Ptree(*args[1:])
    elif args[0]=='full':
        #exclude last argument
        return P(*args[1:-1])
    else:
        raise ValueError('Choose integration=\'full\' or \'tree\'')


class Templates:
       
    def __init__(self, filename, idcolumn, datacolumn, filters, zcolumn=None):
        #read in data and save columns
        data = fits.open(filename)[1].data

        self.filename = filename
        self.data = np.array( zip( *[data[datacolumn.format(f)] for f in filters] ) )
        self.ids = data[idcolumn]
        if zcolumn:
            self.redshifts = data[zcolumn]

        del data
        
    def getMask(self, zrange):
        z_mask = (self.redshifts > np.min(zrange)) & (self.redshifts < np.max(zrange)) 
        return z_mask

        
class Targets:

    def __init__(self, filename, idcolumn, datacolumn, errcolumn, start, end, filters):
        #read in data and save columns
        data = fits.open(filename)[1].data
        
        self.id = idcolumn
        self.filename = filename
        self.data = np.array( zip( *[ data[datacolumn.format(f)][start:end]
                                      for f in filters ] ) )
        self.errors = np.array( zip( *[ data[errcolumn.format(f)][start:end]
                                        for f in filters ] ) )
        self.ids = data[idcolumn][start:end]
        
        del data

    def calcProbabilities(self, galaxies, stars, zranges, numthreads, galaxyintegration):
        N = len(self.data)

        P_dict = {}
        P_dict[self.id] = self.ids

        if galaxyintegration=='tree':
            index_chunks, data_chunks, error_chunks, sigmas = self.kSeparate(numthreads)

            print "#Median errors used in galaxy tree:"
            for sigma in sigmas:
                print sigma, "\n"

        elif galaxyintegration=='full':
            n_per_process = int( np.ceil( len(self.data) / float(numthreads) ) )
            data_chunks = [ self.data[i:i+n_per_process] for i in xrange(0, N, n_per_process) ]
            error_chunks = [ self.errors[i:i+n_per_process] for i in xrange(0, N, n_per_process) ]

            sigmas = None

        else:
            raise ValueError('Choose integration=\'full\' or \'tree\'')

        #multiprocessing
        sys.stderr.write("#Multiprocessing checks:\n")
        sys.stderr.write("    Number of total targets: {}\n".format(N))
        sys.stderr.write("    Number of processes: {}\n".format(numthreads))
        sys.stderr.write("    Number of targets per process: {}\n".format(n_per_process))
        sys.stderr.write("    Number of target chunks, error chunks: {}, {}\n".format(len(data_chunks), len(error_chunks)))

        pool = Pool(processes=numthreads)
        
        for z_range in zranges:
            galaxy_template_data = galaxies.data[ galaxies.getMask(z_range) ]
            sys.stderr.write("Number of templates in redshift range {}: {}\n".format( str(z_range), len(galaxy_template_data) ))

            start_work = time.time()

            results = pool.map(Pwrapper, itertools.izip( itertools.repeat(galaxyintegration), data_chunks, error_chunks, 
                                                         itertools.repeat(galaxy_template_data), itertools.repeat(sigmas)) )
        
            chain_results = np.concatenate(results)
            
            if galaxyintegration=='tree':
                #kSeparate rearranged the order, put it back to match id column
                sorted_back = np.argsort(np.concatenate(index_chunks))
                if len(sorted_back) > len(chain_results):
                    final_results = np.concatenate((chain_results, np.array([np.nan]*(len(sorted_back)-len(chain_results)))))
                    P_dict[str(z_range)] = final_results[sorted_back]
                else:
                    P_dict[str(z_range)] = chain_results
                
            elif galaxyintegration=='full':
                P_dict[str(z_range)] = chain_results

            work_time = time.time() - start_work
            sys.stderr.write("    Work completed in {} s\n".format(work_time))

        pool.close()

        index_chunks, data_chunks, error_chunks, sigmas = self.kSeparate(numthreads)

        print "#Median errors used in tree:"
        for sigma in sigmas:
            print sigma, "\n"
        
        #stars
        pool = Pool(processes=numthreads)

        starintegration = 'tree'
        results = pool.map(Pwrapper, itertools.izip( itertools.repeat(starintegration), data_chunks, error_chunks,
                                                     itertools.repeat(stars.data), sigmas ))
        
        #kSeparate rearranged the order, put it back to match id column
        sorted_back = np.argsort(np.concatenate(index_chunks))
        chain_results = np.concatenate(results)

        if len(sorted_back) > len(chain_results):
            print len(chain_results), len(sorted_back)
            final_results = np.concatenate((chain_results, np.array([np.nan]*(len(sorted_back)-len(chain_results)))))
            P_dict['STAR'] = final_results[sorted_back]
        else:
            P_dict['STAR'] = chain_results

        pool.close()

        return P_dict

    def debug(self, templates, ntest, zranges, numthreads, output):
        debug_ids = self.ids[:ntest]

        #save likelihoods for each template
        z_range_tot = [np.min(np.min(zranges)), np.max(np.max(zranges))]
        z_mask = templates.getMask(z_range_tot)
        template_zs = templates.redshifts[z_mask]
        template_ids = templates.ids[z_mask]
        template_data = templates.data[z_mask]

        col_defs = [fits.Column(name='template_z', format='D', array=template_zs),
                    fits.Column(name='template_id', format='K', array=template_ids)]
        
        for n in range(ntest):
            template_Ls = likelihoods(self.data[n], self.errors[n], template_data)
            col_defs.append(fits.Column(name=str(int(debug_ids[n])), format='D', array=template_Ls))

        pri_hdr = fits.Header()
        tb1_hdr = fits.Header()
        tb1_hdr['COMMENT'] = "Gaussian likelihoods for data in {} \
                             using photo-zs of templates from {}.".format(self.filename, 
                                                                          templates.filename)

        tb1_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_defs), nrows=len(template_data), header=tb1_hdr)

        #save results of different methods
        id_col = fits.Column(name=self.id, format='K', array=debug_ids)

        ##full integration (optimized summing)
        Pfull_dict = self.calcProbabilities(templates, zranges, numthreads, 'full')
        
        ##tree integration
        Ptree_dict = self.calcProbabilities(templates, zranges, numthreads, 'tree', 10000)

        ##sum of likelihoods (no optimization)
        Lsum_dict = {}
        for z_range in zranges:
            template_data = templates.data[templates.getMask(z_range)]
            Lsum_dict[str(z_range)] = np.array([np.sum(likelihoods(self.data[n], self.errors[n], template_data))])

        col_defs = [id_col]
        for k in Pfull_dict.keys():
            col_defs.append(fits.Column(name='Pfull'+k, format='D', array=Pfull_dict[k]))
        for k in Ptree_dict.keys():
            col_defs.append(fits.Column(name='Ptree'+k, format='D', array=Ptree_dict[k]))
        for k in Lsum_dict.keys():
            col_defs.append(fits.Column(name='Lsum'+k, format='D', array=Lsum_dict[k]))
        
        tb2_hdr = fits.Header()
        tb2_hdr['COMMENT'] = "Results from different methods."
        tb2_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_defs), nrows=ntest, header=tb2_hdr)
            
        pri_hdu = fits.PrimaryHDU(header=pri_hdr)
       
        hdu_list = fits.HDUList([pri_hdu, tb1_hdu, tb2_hdu])
        hdu_list.writeto(output, clobber=True)

    def kSeparate(self, numthreads):    
        #kmeans does not take NaNs
        not_nan = ~np.isnan(self.errors)
        zero_errors = self.errors>0.
        e_mask = (np.all(not_nan, axis=1) & np.all(zero_errors, axis=1))
        not_e_mask = ~e_mask

        #separate data by closest errors
        centers, _ = kmeans(self.errors[e_mask], numthreads)
        k_indices, _ = vq(self.errors[e_mask], centers)

        #save indices of chunk pieces so can rearrange later
        index = np.arange(len(self.ids))
        
        index_chunks = [index[e_mask][k_indices==i] for i in range(numthreads)]
        data_chunks = [self.data[e_mask][k_indices==i] for i in range(numthreads)]
        error_chunks = [self.errors[e_mask][k_indices==i] for i in range(numthreads)]

        n_per_process = [len(data_chunk) for data_chunk in data_chunks]

        #median of errors along each filter axis
        sigmas = [np.median(error_chunk.T, axis=1) for error_chunk in error_chunks]

        #account for NaNs
        index_chunks.append(list(index[not_e_mask]))
        
        return index_chunks, data_chunks, error_chunks, sigmas


