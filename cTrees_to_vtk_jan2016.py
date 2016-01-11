""" Consistent trees merger tree to vtk visualization file """

import numpy as np
from pandas import *


def link_halos(ctrees_df):
    " get hierarchy information "
    
    #get absolute halo and descendant IDs
    halo_ids = ctrees_df['id']
    desc_ids = ctrees_df['desc_id']
    """ make 'zeropoint' halo and descendant ID lists
    since absolute IDs are arbitrary and all we care about
    is hierarchy information """
    zeropoint_halo_ids = np.arange(len(halo_ids))
    zeropoint_desc_ids = np.zeros(len(desc_ids))

    for i,halo_id in enumerate(halo_ids):
        # loop over all halo IDs and find all parents of that halo
        parents = np.where(desc_ids==halo_id)[0]
        #print i,parents
        for parent in parents:
            dummy_desc_ids[parent]=i


colnames=['scale','id','desc_scale','desc_id','num_prog','pid','upid',
          'desc_pid','phantom','sam_mvir','mvir','rvir','rs','vrms',
          'mmp?','scale_of_last_MM','vmax','x','y','z','vx','vy','vz',
          'Jx','Jy','Jz','Spin','Breadth_first_ID','Depth_first_ID',
          'Tree_root_ID','Orig_halo_ID','Snap_num',
          'Next_coprogenitor_depthfirst_ID','Last_progenitor_depthfirst_ID',
          'Rs_Klypin','Mvir_all','M200b','M200c','M500c','M2500c','Xoff','Voff',
          'Spin_Bullock','b_to_a','c_to_a','A[x]','A[y]','A[z]',
          'b_to_a(500c)','c_to_a(500c)','A[x](500c)','A[y](500c)','A[z](500c)',
          'T/U','M_pe_Behroozi','M_pe_Diemer','Type','SM','Gas','BH_Mass']


# testing below
df = read_csv('newtree.dat',header=None,comment='#',sep='\s*',engine='python')
df.columns = colnames
link_halos(df)
