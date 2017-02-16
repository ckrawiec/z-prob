import itertools
import time
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
from scipy.spatial import ckdtree

k_near, integrate = None, None

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

    chunks = itertools.izip([vals[i:i+nchunks] for i in xrange(0, len(vals), nchunks)], 
                            [errs[i:i+nchunks] for i in xrange(0, len(vals), nchunks)])
    
    for chunk, errchunk in chunks:
        trueout = np.zeros(len(chunk))
    
        covIs = 1./errchunk**2.
        A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errchunk**2., axis=1 ) )

        truechunks = (truevals[i:i+ntruechunks] for i in xrange(0, len(truevals), ntruechunks))
        for truechunk in truechunks:
            diff = chunk[:,np.newaxis,:]-truechunk[np.newaxis,:,:]

            B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)
            C = A[:,np.newaxis] * np.exp(B)
            
            trueout += np.sum(C, axis=1)

        out = np.concatenate((out, trueout))
    
    return out

def Ptree(vals, errs, truevals):
    
    #median of errors along each filter axis
    sigma = np.median(errs.T, axis=1)

    #build truth tree in units of sigma
    truetree = ckdtree.cKDTree(truevals)#/sigma)

    out = []
    for val, err in zip(vals, errs):

        #get knear nearest neighbors to target
        dnear, inear = truetree.query(val, k=k_near)#/sigma, k=k_near)

        #data of nearest neighbors
        truearr = truetree.data[inear] #* sigma

        Ls = likelihoods(val, err, truearr)

        #sum likelihoods and compensate for tree integral
        out.append(np.sum(Ls) * float(len(truevals))/len(truearr))
            
    return np.array(out)

def Pwrapper(args):
    if integrate=='tree':
        #if truth array shorter than k_near, just do full integration
        if len(args[-1]) > k_near:
            return Ptree(*args)
        else:
            print "Length of truth array ({}) <= {}, using brute force integration instead of tree".format(len(args[-1]), k_near)
            return P(*args)
    elif integrate=='full':
        return P(*args)
    else:
        raise ValueError('Choose integration=\'full\' or \'tree\'')
    

class Templates:
       
    def __init__(self, filename, idcolumn, datacolumn, zcolumn, filters):
        #read in data and save columns
        data = fits.open(filename)[1].data

        self.filename = filename
        self.data = np.array( zip( *[data[datacolumn.format(f)] for f in filters] ) )
        self.ids = data[idcolumn]
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

    def calcProbabilities(self, templates, zranges, numthreads, integration, knear=None):
        global integrate
        integrate = integration
        global k_near
        k_near = knear
    
        N = len(self.data)

        n_per_process = int( np.ceil( len(self.data) / float(numthreads) ) )
        data_chunks = [ self.data[i:i+n_per_process] for i in xrange(0, N, n_per_process) ]
        error_chunks = [ self.errors[i:i+n_per_process] for i in xrange(0, N, n_per_process) ]

        #multiprocessing
        print "#Multiprocessing checks:"
        print "    Number of processes: ", numthreads
        print "    Number of targets per process: ", n_per_process
        print "    Number of target chunks, error chunks: ", len(data_chunks), len(error_chunks)
        print "    Number of total targets: {}\n".format(N)

        P_dict = {}
        pool = Pool(processes=numthreads)
        
        for z_range in zranges:
            template_data = templates.data[ templates.getMask(z_range) ]
            print "Number of templates in redshift range {}: {}".format( str(z_range), len(template_data) )

            start_work = time.time()

            results = pool.map(Pwrapper, itertools.izip( data_chunks, error_chunks, itertools.repeat(template_data) ) )
        
            P_dict[str(z_range)] = np.concatenate(results)

            work_time = time.time() - start_work
            print "    Work completed in {} s".format(work_time)

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
