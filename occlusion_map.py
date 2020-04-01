#!/bin/python

# the goal of this code is to find slices of space and reconstruct into 3d shape.
# apply along axis



import numpy as np
import math  
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
np.set_printoptions(threshold=sys.maxsize)

exclude = 10 # set this value to create an exclusion filter for coordinates that are these many Angstroms away from the N-termini

########## NEXT BLOCK OF CODE ##########

PDB=sys.argv[1]
L_fab=70 # length of proposed Fab used to access Nterm 
W_fab=55 # width of proposed Fab used to access Nterm 
avgVDW=1.6 # estimated length of peptide (avg 3.5 Angstroms per residue)
resolution=1 # math.degrees(math.atan2((W_fab/2),L_fab))=~20degr which is a lot
preArray=[]



# define functions for changing cartesian coordinates to spherical coordinates and vice versa
def sph2cart(sphCo): # [radius, polar, azimuth] in degr
    cPts = np.empty(sphCo.shape)
    pol=np.radians(sphCo[:,1])
    azmt=np.radians(sphCo[:,2])
    cPts[:,0] = sphCo[:,0]*np.sin(pol)*np.cos(azmt) #x
    cPts[:,1] = sphCo[:,0]*np.sin(pol)*np.sin(azmt) # y
    cPts[:,2] = sphCo[:,0]*np.cos(pol) # z
    return cPts

def cart2sph(cartCo): 
    sphPts = np.empty(cartCo.shape)
    xy = cartCo[:,0]**2 + cartCo[:,1]**2
    sphPts[:,0] = np.sqrt(xy + cartCo[:,2]**2) # radius (r)
    sphPts[:,1] = np.degrees(np.arctan2(np.sqrt(xy), cartCo[:,2])) # polar (z) (theta)
    sphPts[:,2] = np.degrees(np.arctan2(cartCo[:,1], cartCo[:,0])) # azimuth (x-y) (phi)
    sphPts=np.where(abs(sphPts) > 1e-4, sphPts,0)
    sphPts[(np.where(sphPts[:,0]==0),1)]=0 # set angles to be 0 if radius is 0
    sphPts[(np.where(sphPts[:,0]==0),2)]=0 
    return sphPts # [radius, polar, azimuth] in degr

# pass in cartesian coordinates of N-terminus and polar coordinates about crtN; effectively a 3-element list and a large mx3 array
def cart2Sph2(crtN, new_origin): 
	var2=np.subtract(crtN,new_origin)
	sph2C=cart2sph(var2)
	return sph2C

# calculate radius of each slice
def sliceR (angle_value):
	R_slice=[] # distances of slices from sphere center
	for i in range(angle_value): # for each number of slices
		R_slice.append(math.sqrt((angle_value)**2 - i**2))
	R_slice.append(0)
	return R_slice # list of diameters of each slice from center to tip

# define the area occluded by nearby atoms
def circleblock(rad, ang):
	theta=ang[0]
	phi=ang[1]
	t=[]
	p=[]
	arctan=int(round(math.degrees(math.atan2(avgVDW,rad)),0))
	
	if arctan % 2 == 1:
		j = 0
	else:
		j = 1
	for i in range(-arctan,arctan+1):
		t.append(theta+i)
		p.append(phi+i)

	R_span_list=sliceR(arctan)
	R_span_bk=np.rint(np.array(R_span_list)+arctan+j).astype(int)
	
	R_span_list.pop(0)
	R_span_list.reverse()
	R_span_fw=np.rint(np.array(R_span_list)+arctan+j).astype(int)
	frq=np.concatenate((R_span_fw, R_span_bk))
	values_t=[]
	values_p=[]
	
	for k in range(len(frq)):
		for m in range(frq[k]):
			values_t.append(t[k])
			values_p.append(p[((len(p)-frq[k])/2)+m])
	theta=(np.array(values_t)).reshape([len(values_t),1])
	phi=np.array(values_p).reshape([len(values_p),1])

	return np.hstack((theta,phi))

# to keep all angles within mapped bounds of 180 and 360
def withinrange(rawangles):
	processed=np.empty(np.shape(rawangles))
	tlist=rawangles[:,0]
	plist=rawangles[:,1]
	tlist[tlist > 179] -=360
	tlist[tlist < 0] *= -1
	tlist[tlist == 180] = 0
	plist[plist > 359] -= 360
	plist[plist < 0] += 360
	processed[:,0]=tlist
	processed[:,1]=plist
	return processed

# map-making of occluded space from POV of specified N-termini

print "making figures..."

with open("<folders>/orgCoord.list") as f:
	coord=f.readlines() 

with open("<folders>blacklist.list") as g:
	blacklist=np.loadtxt(fname = "<folders>blacklist.list")

with open("<folders>residues.list") as h:
	h.readline()
	residues=(h.readline()).split()

# make combinations of 0-360 degr in both theta and phi angles
box = np.mgrid[0:180:resolution,0:360:resolution]
alldegr=(np.dstack(box)).reshape(64800, 2)
xcdn=np.loadtxt(fname = "<folders>/XCOR.txt")
ycdn=np.loadtxt(fname = "<folders>/YCOR.txt")
zcdn=np.loadtxt(fname = "<folders>/ZCOR.txt")
cartCo=np.empty([len(xcdn),3])
cartCo[:,0]=xcdn
cartCo[:,1]=ycdn
cartCo[:,2]=zcdn
cartCo=[i for i in cartCo if i not in blacklist]


# create arrays of coordinates for each PDB structure
for line in coord:
	if line.strip(): # if line is not empty
		if "/" not in line:
			xyz=[];
			xyz.append(float(line.split()[0]))
			xyz.append(float(line.split()[1]))
			xyz.append(float(line.split()[2]))
			preArray.append(xyz)		 
	else:
		# array of all Nterm coordinates in specific PDB
		cartOrgs=np.array(preArray)
		
		# open the file with all the spherical coordinates and store in array
		for x in range(len(cartOrgs)):
			total=(alldegr).tolist()

			# find radius values withi L_fab away from Nterm's distance from sph1 origin
			#take these coordinates and switch them to sph2 FoR
			co_sph2=cart2Sph2(cartCo,cartOrgs[x]) #masked values to sph2 FoR #len=47850

			#check first column and filter for r values that are w/in L_fab of Nterm and return entire row
			mask1=co_sph2[:,0] > L_fab	
			m1=np.column_stack((mask1,mask1,mask1))
			m_L_fab=np.ma.compressed(np.ma.array(co_sph2, mask=m1))
			data_sph2_L_fab=(np.reshape(m_L_fab, [len(m_L_fab)/3,3]))

			if exclude > 0:
				mask2=data_sph2_L_fab[:,0] < exclude
				m2=np.column_stack((mask2,mask2,mask2))
				m_exclude=np.ma.compressed(np.ma.array(data_sph2_L_fab,mask=m2))
				data_sph2_L_fab=(np.reshape(m_exclude,[len(m_exclude)/3,3]))
			
			radii=data_sph2_L_fab[:, 0] 

			rawangles=np.empty([len(data_sph2_L_fab),2])
			rawangles[:,0]=np.rint(data_sph2_L_fab[:, 1]) 
			rawangles[:,1]=np.rint(data_sph2_L_fab[:, 2])
			rawangles[rawangles < 0] +=360

			
			# collect all points that are occluded and map. 
			for k in range(len(rawangles)):
				sizeblock=circleblock(radii[k], rawangles[k])
				rawangles=np.concatenate((rawangles, sizeblock), axis=0)
			
			processed=withinrange(rawangles)

			blocked=[tuple(row) for row in rawangles]
			inaccessible = np.unique(blocked)
			inaccessible=(inaccessible).tolist()

			for pair in inaccessible:
				total.remove(pair)	
			possible=np.array(total)	

			print "space available: " + str((np.shape(possible)[0])/648)+ "%"

			########## NEXT BLOCK OF CODE ##########
			
			thetax=possible[:,0]
			phiy=possible[:,1]
			
			plt.scatter(phiy, thetax, c='w', s=20, edgecolors='w')
			axes=plt.gca()
			axes.set_xlim([0,360])
			axes.set_ylim([180,0])
			axes.set_axis_off()
			plt.margins(0,0)
			axes.xaxis.set_major_locator(tick.NullLocator())
			axes.yaxis.set_major_locator(tick.NullLocator())


			plotname=('figures/{0}_{1}.png'.format(PDB,residues[x]))
			fig=plt.gcf()
			fig.set_size_inches(20, 10)
			fig.savefig(plotname, dpi=600, bbox_inches='tight', pad_inches=0.0, facecolor="black", edgecolor='k')

			plt.close()
					
			print("saving Figure " + str(x+1) + " of " + str(len(cartOrgs)))
			
			########## NEXT BLOCK OF CODE ##########
			
			

		print("success")

