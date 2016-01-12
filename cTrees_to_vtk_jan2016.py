""" Consistent trees merger tree to vtk visualization file """
from __future__ import division
import numpy as np
from pandas import *

""" define column names from consistent trees Version 0.99.9.3 files """
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

def link_halos(ctrees_df):
    """ Get hierarchy information.
    Expecting a pandas dataframe with
    one consistent trees merger tree in it. """
    
    """ get absolute halo and descendant IDs """
    halo_ids = ctrees_df['id']
    desc_ids = ctrees_df['desc_id']
    """ make 'zeropoint' halo and descendant ID lists
    since absolute IDs are arbitrary and all we care about
    is hierarchy information """
    zeropoint_halo_ids = np.arange(len(halo_ids))
    zeropoint_desc_ids = np.zeros(len(desc_ids))

    for i,halo_id in enumerate(halo_ids):
        """ loop over all halo IDs and find all parents of that halo """
        parents = np.where(desc_ids==halo_id)[0]
        #print i,parents
        for parent in parents:
            zeropoint_desc_ids[parent]=i
  
    print 'Halos linked!'
    return zeropoint_halo_ids,zeropoint_desc_ids

def vtk_write(ctrees_df,zeropoint_halo_ids,zeropoint_desc_ids,ofpath,params):
    """ Write vtk file with all necessary hierarchy + parameter information. """
    
    # definitions
    npoints = ctrees_df.shape[0]
    nparams = len(params)

    #===================================================================
    # bare minimum for vtk file
    #===================================================================
    # open output file
    ofile = open(ofpath,'w+')
    # write necessary header info
    ofile.write('# vtk DataFile Version 3.0\n')
    ofile.write('vtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n')
    
    # x,y,z coordinates
    ofile.write('POINTS %d float\n' %npoints)
    for irow in range(npoints):
        # get and write spatial coordinates
        row = ctrees_df.loc[irow]
        x,y,z = row['x'],row['y'],row['z']
        ofile.write('%f %f %f\n' %(x,y,z))
    
    # hierarchy info
    ofile.write('CELLS %d %d\n' %(npoints,npoints*3))
    for irow in range(npoints):
        # get and write halo and desc IDs
        halo_id,desc_id = zeropoint_halo_ids[irow],zeropoint_desc_ids[irow]
        ofile.write('2 %d %d\n' %(halo_id,desc_id))

    # vtk cell types parameter (easiest to just leave at 3 for our purposes)
    ofile.write('CELL_TYPES %d\n' %npoints)
    for irow in range(npoints):
        ofile.write('3\n')
    ofile.write('\n')

    #===================================================================
    # from here down, we write the optional parameters to the vtk file
    #===================================================================
    ofile.write('POINT_DATA %d\n' %npoints)
    ofile.write('FIELD FieldData %d\n' %nparams)
    
    # redshift
    if 'redshift' in params:
        # get redshift from scale factor using the familiar formula
        z = (1/ctrees_df['scale']) - 1
        ofile.write('redshift 1 %d double\n' %npoints)
        for irow in range(npoints):
            ofile.write('%f\n' %z[irow])
    
    # virial mass
    if 'mvir' in params:
        # get log(M_vir)
        logmvir = np.log10(ctrees_df['mvir'])
        ofile.write('log(M_vir) 1 %d double\n' %npoints)
        for irow in range(npoints):
            ofile.write('%f \n' %logmvir[irow])

    # virial radius
    if 'rvir' in params:
        # get log(R_vir)
        logrvir = np.log10(ctrees_df['rvir'])
        ofile.write('log(R_vir) 1 %d double\n' %npoints)
        for irow in range(npoints):
            ofile.write('%f \n' %logrvir[irow])

    # velocity dispersion
    if 'vrms' in params:
        ofile.write('V_rms 1 %d double\n' %npoints)
        for irow in range(npoints):
            ofile.write('%f\n' %ctrees_df.loc[irow]['vrms'])

    # spin parameter
    if 'Spin' in params:
        # get log(Spin)
        logspin = np.log10(ctrees_df['Spin'])
        ofile.write('log(spin) 1 %d double\n' %npoints)
        for irow in range(npoints):
            ofile.write('%f\n' %logspin[irow])
    
    # most massive progenitor?
    if 'mmp?' in params:
        ofile.write('MMP? 1 %d int\n' %npoints)
        for irow in range(npoints):
            ofile.write('%i\n' %ctrees_df.loc[irow]['mmp?'])

    # phantom?
    if 'phantom' in params:
        ofile.write('phantom? 1 %d int\n' %npoints)
        for irow in range(npoints):
            ofile.write('%i\n' %ctrees_df.loc[irow]['phantom'])
    
    # velocity vectors
    if 'velocity' in params:
        ofile.write('velocity 3 %d double\n' %npoints)
        for irow in range(npoints):
            # get velocity components
            row = ctrees_df.loc[irow]
            vx,vy,vz = row['vx'],row['vy'],row['vz']
            ofile.write('%f %f %f\n' %(vx,vy,vz))
    
    # angular momentum
    if 'angular_momentum' in params:
        ofile.write('angular_momentum 3 %d double\n' %npoints)
        for irow in range(npoints):
            # get ang_mom components
            row = ctrees_df.loc[irow]
            jx,jy,jz = row['Jx'],row['Jy'],row['Jz']
            ofile.write('%f %f %f\n' %(jx,jy,jz))
    
    # ellipsoid axes
    if 'ellipsoid_shape' in params:
        ofile.write('ellipsoid_axes 3 %d double\n' %npoints)
        for irow in range(npoints):
            # get geometric components
            row = ctrees_df.loc[irow]
            Ax,Ay,Az = row['A[x]'],row['A[y]'],row['A[z]']
            ofile.write('%f %f %f\n' %(Ax,Ay,Az))


    ofile.close()

if __name__ == '__main__':
    # testing below
    df = read_csv('newtree.dat',header=None,comment='#',sep='\s*',engine='python')
    df.columns = colnames
    zp_hids,zp_dids = link_halos(df)
    params=['mvir','velocity','redshift','rvir','angular_momentum',
            'vrms','mmp?','ellipsoid_shape','phantom','Spin']
    vtk_write(df,zp_hids,zp_dids,'test_vtk_jan2016.vtk',params)


