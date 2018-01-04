#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2010-2015 CERN and Helsinki Institute of Physics.
# This software is distributed under the terms of the
# GNU General Public License version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md. In applying this
# license, CERN does not waive the privileges and immunities granted to it
# by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.
#
# Project website: http://arcpic.web.cern.ch/
# Developers: Helga Timko, Kyrre Sjobak
#
# 2Dpic_fieldMap.py:
# Plots field in gap as function of r, z at different time steps
#


import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt

import h5py

from matplotlib import rcParams
rcParams.update({'text.usetex': True})

#Get input
if len(sys.argv) != 4 and len(sys.argv) != 7:
    print "Usage: ./fieldMap.py <mintime> <maxtime> <every nth frame to analyse> [fieldVminR:fieldVmaxR fieldVminZ:fieldVmaxZ fieldVminAbs:fieldVmaxAbs]"
    print "Set min:max to '-:-' to use automatic/per-frame scaling"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

fieldVminR   = None
fieldVmaxR   = None
fieldVminZ   = None
fieldVmaxZ   = None
fieldVminAbs = None
fieldVmaxAbs = None
if len(sys.argv) == 7:
    minmax = sys.argv[4].split(":")
    if len(minmax) != 2:
        print "didn't understand fieldVminR:fieldVmaxR, got '" + sys.argv[4] + "'"
        exit(1)
    if minmax[0] == minmax[1] and minmax[0]=="-":
        fieldVminR = None
        fieldVmaxR = None
    else:
        fieldVminR = float(minmax[0])
        fieldVmaxR = float(minmax[1])

    minmax = sys.argv[5].split(":")
    if len(minmax) != 2:
        print "didn't understand fieldVminZ:fieldVmaxZ, got '" + sys.argv[5] + "'"
        exit(1)
    if minmax[0] == minmax[1] and minmax[0]=="-":
        fieldVminZ = None
        fieldVmaxZ = None
    else:
        fieldVminZ = float(minmax[0])
        fieldVmaxZ = float(minmax[1])

    minmax = sys.argv[6].split(":")
    if len(minmax) != 2:
        print "didn't understand fieldVminAbs:fieldVmaxAbs, got '" + sys.argv[6] + "'"
        exit(1)
    if minmax[0] == minmax[1] and minmax[0]=="-":
        fieldVminAbs = None
        fieldVmaxAbs = None
    else:
        fieldVminAbs = float(minmax[0])
        fieldVmaxAbs = float(minmax[1])

#Get the scaling setup
basefile = h5py.File("../out/output_00000000.h5",'r')

Ldb   = basefile['/METADATA/CALCULATED'].attrs['Ldb']
O_pe  = basefile['/METADATA/CALCULATED'].attrs['O_pe']
dZ    = basefile['/METADATA/CALCULATED'].attrs['dZ']
R     = basefile['/METADATA/CALCULATED'].attrs['R']
Z     = basefile['/METADATA/CALCULATED'].attrs['Z']
nr    = basefile['/METADATA/INPUTFILE'].attrs['nr']
nz    = basefile['/METADATA/INPUTFILE'].attrs['nz']

basefile.close()

rFactor   = 1e4*Ldb #output units -> um
fieldFactor = 510.998928e3/2.99792458e10**2 * Ldb * O_pe**2 * 1e-4 #output units -> MV/m

#Create output folder
ofoldername = "pngs/fieldMap"
if os.path.exists(ofoldername):
    if os.path.isdir(ofoldername):
        print "Removing " + ofoldername
        shutil.rmtree(ofoldername)
    else:
        print "Path '" + ofoldername + "' exists, but is not a directory. Aborting!"
        exit(1)
os.makedirs(ofoldername)
print "Created directory '%s'" % (ofoldername,)

print "Got options:"
print " - mintime     =", mintime
print " - maxtime     =", maxtime
print " - skipFrame   =", skipFrame
print
print " - fieldVminR     =", fieldVminR
print " - fieldVmaxR     =", fieldVmaxR
print
print " - fieldVminZ     =", fieldVminZ
print " - fieldVmaxZ     =", fieldVmaxZ
print
print " - fieldVminAbs   =", fieldVminAbs
print " - fieldVmaxAbs   =", fieldVmaxAbs
print


#Get list of output files and timestamp
stepNum   = []
timestamp = [] #[ns]
timeIdx = open("../out/timeIndex.dat")
timeIdx.readline() #Skip first line
if mintime == 0:
    stepNum.append("00000000")
    timestamp.append(0.0)
for line in timeIdx:
    ls = line.split()
    stepNum.append(ls[0])
    timestamp.append(float(ls[1]))
timeIdx.close()

QuiverScale      = None
QuiverScale_zoom = None

### Time step loop ###
pCyclic = 0; #Used for skipping
outIdx = 0;
firstTime = None
finalTime = None
for i in xrange(len(stepNum)):
    ts = int(stepNum[i])
    
    if ts < mintime:
        continue
    elif ts > maxtime:
        print "ts =", ts, "reached maxtime =", maxtime
        break
    elif pCyclic % skipFrame != 0:
        print "skipFrame", ts
        pCyclic += 1
        continue
    fname = "../out/output_" + stepNum[i] + ".h5"
    print "Plotting ts=", ts, ",", timestamp[i], "[ns], fname='" + fname + "'"
    
    #Read data for this timestep
    def readFieldFile(fname):
        datafile = h5py.File(fname, 'r')
        fieldData = datafile['EMFIELD/EFIELD']

        rData = fieldData[:,:,0]*rFactor
        zData = fieldData[:,:,1]*rFactor
        
        fieldList_z = fieldData[:,:,2]*fieldFactor
        fieldList_r = fieldData[:,:,3]*fieldFactor
        
        return (rData, zData, fieldList_z, fieldList_r)

    (rList, zList, fieldList_z, fieldList_r) = readFieldFile(fname)
    
    #assert rList_z.all() == rList_r.all()
    #assert zList_z.all() == zList_r.all()
    
    #plot Ez
    plt.contourf(zList,rList,fieldList_z,100, vmin=fieldVminZ, vmax=fieldVmaxZ);
    plt.colorbar()
    
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    plt.title("Z-field [MV/m], time = %.3f ns" % (timestamp[i],))

    plt.savefig(ofoldername + "/Ez_%08d.png" %(outIdx,),dpi=300)
    plt.clf()

    #plot Er
    plt.contourf(zList,rList,fieldList_r,100, vmin=fieldVminR, vmax=fieldVmaxR);
    plt.colorbar()
    
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    plt.title("R-field [MV/m], time = %.3f ns" % (timestamp[i],))

    plt.savefig(ofoldername + "/Er_%08d.png" %(outIdx,),dpi=300)
    plt.clf()

    #Field strength
    field = np.sqrt(fieldList_z**2 + fieldList_r**2)
    plt.contourf(zList,rList,field,100, vmin=fieldVminAbs, vmax=fieldVmaxAbs);
    plt.colorbar()
    
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    plt.title("Field magnitude [MV/m], time = %.3f ns" % (timestamp[i],))

    plt.savefig(ofoldername + "/Eabs_%08d.png" %(outIdx,),dpi=300)
    plt.clf()

    #Quiver plot
    Q = plt.quiver(zList, rList, fieldList_z, fieldList_r, field, scale=QuiverScale,
                   clim=(fieldVminAbs, fieldVmaxAbs))

    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    plt.title("Field, time = %.3f ns" % (timestamp[i],))
    
    plt.colorbar().set_label("Magnitude [MV/m]")

    plt.xlim(0.0-dZ*1e4,Z*1e4+dZ*1e4)
    plt.ylim(0.0-dZ*1e4,R*1e4+dZ*1e4)

    plt.savefig(ofoldername + "/quiver_%08d.png" %(outIdx,),dpi=300)
    if QuiverScale == None:
        QuiverScale = Q.scale
    plt.clf()
    
    #Quiver (zoomed)
    maxZ = Z*1e4/4.0
    maxR = R*1e4/3.0
    # find maximal indices to use
    iZ=0
    for j in xrange(nz+1):
        if zList[0,j] > maxZ:
            iZ = j
            break
    assert iZ != 0
    iR = 0
    for j in xrange(nr+1):
        if rList[j,0] > maxR:
            iR = j
            break
    assert iR != 0

    Q = plt.quiver(zList[:iR,:iZ], rList[:iR,:iZ],
                   fieldList_z[:iR,:iZ], fieldList_r[:iR,:iZ], field[:iR,:iZ],
                   scale=QuiverScale_zoom, clim=(fieldVminAbs, fieldVmaxAbs))
    plt.xlim(0.0-3*dZ*1e4,maxZ+dZ*1e4)
    plt.ylim(0.0-3*dZ*1e4,maxR+dZ*1e4)

    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    plt.title("Field, time = %.3f ns" % (timestamp[i],))
 
    plt.colorbar().set_label("Magnitude [MV/m]")

    plt.savefig(ofoldername + "/quiver_zoom_%08d.png" %(outIdx,),dpi=300)
    if QuiverScale_zoom == None:
        QuiverScale_zoom = Q.scale
    plt.clf()
    
    pCyclic += 1
    outIdx += 1    

#Movie assembly
def movieAssembly(movieName, imageBasename, speed):
    movieFileName = ofoldername + "/" + os.path.basename(movieName) + ".mp4"
    if mintime>0:
        out_dt = (timestamp[2]-timestamp[1])*skipFrame
    else:
        out_dt = (timestamp[1]-timestamp[0])*skipFrame

    fps = speed/out_dt
    fps = int(fps)
    if fps == 0:
        fps = 1
    print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
    
    ffmpegCommand = "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/" + imageBasename + "_%08d.png " + movieFileName
    print "Command: '" + ffmpegCommand
    
    print "(skipping actually running the command)"
    #os.system("rm " + movieFileName)
    #os.system(ffmpegCommand)

    #Make some space after ffmpeg output
    print
    print
speed = 0.25 #ns/s

#movieAssembly("Ez","Ez", speed)
#movieAssembly("Er","Er", speed)
#movieAssembly("Eabs","Eabs", speed)
#movieAssembly("quiver","quiver", speed)
#movieAssembly("quiver_zoom","quiver_zoom", speed)
