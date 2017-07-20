# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

import sys
import argparse
import pickle

from Geometry import *
from Configuration import *
from Tesselation import *
from Defects import *
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("--contractile", action='store_true',default=False, help="Adds contractile stresses to calculation")
parser.add_argument("-a", "--alpha", type=float, default=0.0, help="Prefactor of the contractile term")
parser.add_argument("--delaunay", action='store_true',default=False, help="Use Delaunay triangulation for tesselation.")
parser.add_argument("--closeHoles", action='store_true', default=False, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
parser.add_argument("-m", "--mult",type = float, default=1.0, help="initial Multiplier for tesselation neighbour radius")
parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
parser.add_argument("--getStatsBasic",action='store_true',default=False, help="Output basic stats (v2av, packing fraction, energy, pressure)")
parser.add_argument("--prefix",type=str,default='frame',help="prefix of vtp output files and pickle output")
args = parser.parse_args()


params = Param(args.directory+args.conffile)
print params

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]

if len(files) == 0:
  files = sorted(glob(args.directory + args.input+'*.dat.gz'))[args.skip:]

if args.writeD:
	defects_n_out=[[] for u in range(len(files))]
	defects_v_out=[[] for u in range(len(files))]
	numdefects_n_out=np.zeros(len(files))
	numdefects_v_out=np.zeros(len(files))
	
if args.getStatsBasic:
	vel2av=np.zeros(len(files))
	phival=np.zeros(len(files))
	ndensity=np.zeros(len(files))
	pressure=np.zeros(len(files))
	fmoment=np.zeros(len(files))
	energy=np.zeros(len(files))
	energytot=np.zeros(len(files))
	zav=np.zeros(len(files))

u=0
for f in files:
	print f
	# Ignore here is to simply calculate interactions of multiple types of particles (they have the same potential)
	conf = Configuration(params,f,True)
	if args.contractile:
		writeme = Writer(args.nematic,args.alpha)
	else:
		writeme = Writer(args.nematic)
	if args.writeP:
		outparticles = args.output + '/'+ args.prefix + '%06d_particles.vtp' % u  # + str(u) + '_particles.vtp'
		print outparticles
		writeme.writeConfigurationVTK(conf,outparticles)
	#plt.show()
	if args.getStatsBasic:
		vel2av[u], phival[u],ndensity[u], pressure[u],fmoment[u],energy[u],energytot[u],zav[u]=conf.getStatsBasic()
		#vel2av[u], phival[u],pressav[u],energy[u]= conf.getStatsBasic()
		print "Mean square velocity: " + str(vel2av[u])
		print "Packing fraction: " + str(phival[u])
		print "Number density: " + str(ndensity[u])
		print "Pressure: " + str(pressure[u])
		print "Mean force  moment: " + str(fmoment[u])
		print "Energy per particle: " + str(energy[u])
		print "Total energy: " + str(energytot[u])
		print "Contact number: " + str(zav[u])
	if args.writeD or args.writeT:
		conf.getTangentBundle()
		tess = Tesselation(conf)
		print "initialized tesselation"
		if args.delaunay:
                        LoopList,Ival,Jval = tess.findLoopDelaunay()
                else:
                        LoopList,Ival,Jval = tess.findLoop(args.closeHoles,args.mult,1.1)
		print "found loops"
		#print LoopList
		if args.writeD:
			#print "Still to be done ..."
			outdefects = args.output + '/' + args.prefix + '%06d_defects.vtp' % u # + str(u) + '_defects.vtp'	
			print outdefects
			defects = Defects(tess,conf)
			# Look for nematic defects in the director field, but do not look for velocity defects (since it's a mess)
			if args.nematic:
				defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('nematic',False)
			else:
				defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('polar')
			print "found defects"
			defects_n_out[u]=defects_n
			defects_v_out[u]=defects_v
			numdefects_n_out[u]=numdefect_n
			numdefects_v_out[u]=numdefect_v
			#defects.PlotDefects()
			writeme.writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outdefects)
		if args.writeT:
			outpatches = args.output + '/' + args.prefix + '%06d_patches.vtp' % u #+ str(u) + '_patches.vtp'
			print outpatches
			if args.makeEdges:
				tess.makeEdges(3.0)   
			tess.OrderPatches()
			print "ordered patches"
			if args.contractile:
				writeme.writePatches(tess,outpatches,True)
			else:
				writeme.writePatches(tess,outpatches,True)
	
	u+=1
#This is so that visualisation only doesn't overwrite previously analyzed data
if ((args.writeD) or (args.getStatsBasic)):
	try:
		data={'J':params.J,'v':params.v0,'k':params.pot_params['k'],'pot_params':params.pot_params,'population':params.population,'pop_params':params.pop_params}
	except:
                try:
                    data={'J':params.J,'v':params.v0,'k':params.pot_params['k'],'pot_params':params.pot_params}
                except:
                    data={'pot_params':params.pot_params}
	if args.writeD:
		dataD={'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out,'numdefects_v':numdefects_v_out}
		data.update(dataD)
	if args.getStatsBasic:
		dataS={'vel2av':vel2av,'phival':phival,'ndensity':ndensity,'pressure':pressure,'fmoment':fmoment,'energy':energy,'energytot':energytot,'zav':zav}
		data.update(dataS)
	if args.writeD:	
		if args.prefix=='frame':
			outpickle=args.directory+'defect_data.p'
		else:
			outpickle=args.directory+args.prefix+'_defect_data.p'
	else:
		if args.prefix=='frame':
			outpickle=args.directory+'configuration_data.p'
		else:
			outpickle=args.directory+args.prefix+'_configuration_data.p'
	print outpickle
	pickle.dump(data,open(outpickle,'wb'))
	
	
