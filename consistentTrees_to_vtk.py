"""
consistent_trees_to_vtk_v1.0
cTrees_to_vtk_v1.py
Written by Yotam Cohen in python 2.7
Must have numpy installed.
This script takes as input an output file from Consistent Trees Version 0.99.9.3
and makes a .vtk file out of it. vtk files can then be visualized using paraview
or any other visualization software that accepts vtk files.
Last update: 3/29/15
"""

from sys import argv, exit
import sys
import numpy as np
from pandas import *

def linkHaloes(datafile):
    hids, dids = [], []
    #for line in datafile[1:]:
    for line in datafile:
        hids.append(line.split()[1])
        dids.append(line.split()[3])
    hids = np.array(hids)
    dids = np.array(dids)

    zphids = np.arange(len(hids))
    zpdids = np.zeros(len(dids))
    
    for i,hid in enumerate(hids):
        matches = np.where(dids==hid)[0]
        for x in matches:
            zpdids[x]=i
    return zphids, zpdids

#def vtkprint(datafile,hids,dids):
def vtkprint(datafile):
    print '# vtk DataFile Version 3.0'
    print 'vtk output'
    print 'ASCII'
    print 'DATASET UNSTRUCTURED_GRID'

    print 'POINTS %d float' %(npoints) # points info
    for line in datafile:
        l = map(float,line.split())
        print l[17], l[18], l[19] # print x,y,z
    
    #print 'CELLS %d %d' %( npoints-1, (npoints-1)*3 ) # cells info
    print 'CELLS %d %d' %( npoints, npoints*3 ) # cells info
    for i,line in enumerate(datafile): # ignore first line (final single halo)
    #for line in datafile[1:]: # ignore first line (final single halo)
        l = map(float,line.split()) 
        hid, did = zphids[i], zpdids[i]
        print '2', int(hid), int(did) # print 2, halo_id, desc_id
        #try:
            #hid, did = zphids[i], zpdids[i]
            #print '2', int(hid), int(did) # print 2, halo_id, desc_id
        #except NameError:
            #print '2', int(l[1]), int(l[3]) # print 2, halo_id, desc_id

    #print 'CELL_TYPES %d' %(npoints-1) # cell types
    print 'CELL_TYPES %d' %(npoints) # cell types
    for line in datafile: # ignore first line (final single halo)
    #for line in datafile[1:]: # ignore first line (final single halo)
        print '3' # cell type (3 = line)

    # From here on down, all information is optional
    print ''
    print 'POINT_DATA %d' %(npoints)
    print 'FIELD FieldData %d' %(nfields)

    # do this next part for all relevant fields
    # print 'parameterName numComponents nPoints dataType' %(1 or 3, numpoints)
    print 'LogHaloMass 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print np.log10(l[10])

    print 'Velocities 3 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[20], l[21], l[22] # vx, vy, vz
        
    print 'Redshift 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        z = (1./l[0])-1.
        print z # z

    print 'Phantom 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[8] # bool

    print 'LogRvir 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print np.log10(l[11]) # rvir
    
    print 'MMP 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[14] # bool

    print 'OrigHaloID 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[30]

    print 'Spin 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[26]
    
    """
    print 'NumProgenitors 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[4]
    
    print 'ScaleFactor 1 %d double' %(npoints)
    for line in datafile:
        l = map(float,line.split())
        print l[0]
    """

if __name__ == '__main__':

    try:
        fname = argv[1]
    except IndexError:
        try:
            fname = raw_input('consisent trees file: ')
            f = open(fname)
        except IOError:
            print '%s does not exist in this directory' %fname
            sys.exit(0)
    
    with open(fname) as f:
        foo = f.readlines()
    del f
    treeNums, treeLines, headerLines, data = [],[],[],[]
    for i,line in enumerate(foo):
        if line.startswith('#') and line.startswith('#tree')==False:
            headerLines.append(line)
        elif line.startswith('#tree'):
            treeNums.append(line.split()[1])
            treeLines.append(i)
        elif len(map(float,line.split())) != 1:
            data.append(line)
    treeSizes = [ float(treeLines[i+1])-float(treeLines[i]) for i in xrange( len(treeLines[:-1]) ) ]
    treeSizes = map(int,treeSizes)

    if len(treeNums) > 1:
        if len(treeNums) > 10:
            listall = raw_input('%d trees in this file. List all? (y/n): ' %len(treeNums))
            if listall=='y':
                print 'tree ID \t number of haloes'
                for i,num in enumerate(treeNums): print treeNums[i],'\t','\t',treeSizes[i]
            else:
                howmany = int(raw_input('List how many?: '))
                print 'tree ID \t number of haloes'
                for i in xrange(howmany): print treeNums[i],'\t','\t',treeSizes[i]
        else:
            print 'tree_ID number_of_haloes'
            for i,num in enumerate(treeNums): print treeNums[i],'\t','\t',treeSizes[i]

        whichtree = raw_input('Tree number to visualize (0 for largest/parent tree): ')
        if whichtree not in treeNums:
            print 'Not a valid tree... exiting...'
            exit()
        elif whichtree != 0:
            dexnum = treeNums.index(whichtree)
            startline = treeLines[dexnum]+1
            endline = treeLines[dexnum+1]
            datalines = foo[startline:endline]
        elif whichtree == 0:
            startline = treeLines[0]+1
            endline = treeLines[1]
            datalines = foo[startline:endline]    
    else:
         whichtree = treeNums[0]
         datalines = data

    try:
        ofname = argv[2]
    except IndexError:
        #ofname = 'tree_%s.vtk' %whichtree
        ofname = 'tree'
        print 'Creating output file %s' %ofname
     
    dl = datalines
    with open('lines.txt','w+') as f:
        for i,line in enumerate(dl):
            f.write(line)
            #print i
    df = read_csv('lines.txt',header=None,sep='\s*',engine='python')
    print df.shape
    timesteps = np.unique(df[0])
    for i,step in enumerate(timesteps):
        dfsub = df[df[0]==step]
        temp = np.array(dfsub)
        temp = temp.tolist()
        newdl = []
        for line in temp:
            line = map(str,line)
            newl =  ' '.join(line)
            newdl.append(newl)
        print len(newdl)
        #off = '%s-s%d' %(ofname,i)
        off = '%s_%d.vtk' %(ofname,i)
        print off
        sys.stdout = open(off, 'w+')
        nfields = 8
        npoints = len(newdl)
        zphids, zpdids = linkHaloes(newdl)
        vtkprint(newdl)
        sys.stdout = sys.__stdout__
        """
        npoints = len(datalines)
        zphids, zpdids = linkHaloes(datalines)
        vtkprint(datalines)
        """
