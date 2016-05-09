import pandas as pd
from astropy.io import ascii as astropyascii
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams

######
# M-dwarf FITS files in './lightcurvedata/' + (9-digit Kepid)+'/'
# Synthetic FITS files in './Mdwarfsynthetics/'
######

# Amount of time before/after TCE midtransit for 2-panel (length=2) or 3-panel (lengh=3) figures
buffertimes=[0.5,1.5,5]

writefigures=True        # If true, writes .png files
writetcefiles=True       # If true, writes the TCE files with TCE labels
weighting=True           # If true, uses the user weighted data
plot=True                # If true and 'writefigures==False', plot pops up
if weighting==True:      # Picks which TCE filename to choose
	filewritename='weighted'
elif weighting==False:
	filewritename='nonweighted'

# List of M-dwarf (real data) TCEs
mdwarftces=astropyascii.read('allmdwarf'+filewritename+'tce.dat',header_start=0,data_start=1,delimiter=' ')
mdwarftces.sort('userxmid')
mdwarftces.sort('kepid')

# List of synthetic TCEs
syntces=astropyascii.read('syn'+filewritename+'tce.dat',header_start=0,data_start=1,delimiter=' ')
syntces.sort('userxmid')
syntces.sort('kepid')

tceid=0  # Same numbering system applies to synthetics and M-dwarf TCEs
alltces=[mdwarftces,syntces]  # Pseudo-combines them into one list
for i in range(len(alltces)):
#for i in [0]:
	tces=alltces[i]  # Pick either M-dwarf or synthetic TCE array
	tceidcolumn=[]   # List to add to TCE array to write output including TCE ID
	if i==0:         # Which folder name to put output
		datatype='mdwarfs'
	elif i==1:
		datatype='synthetics'
	for j in range(len(tces)):  # Going through each TCE
	#for j in [0]:
		# Assiging a TCE ID to the data
		tceid+=1
		tceidcolumn.append(tceid)
		tce=tces[j]
		# Extracting transit information
		kepid=tce['kepid']
		userxmid=tce['userxmid']
		userxmin=tce['userxmin']
		userxmax=tce['userxmax']
		if i==0:    # M-dwarf TCEs
			quarter=tce['datalocation'].split('_')[-1].split('.')[0]
			syntheticid='N/A'
			locationsplit=tce['fits'].split('/')
			fitslocation='./lightcurvedata/'+locationsplit[-2]+'/'+locationsplit[-1]
		elif i==1:  # Synthetic TCEs
			quarter=tce['datalocation'].split('_')[-1].split('-')[0]
			syntheticid=str(tce['syntheticid'])
			fitslocation='./Mdwarfsynthetics/synthetic_'+syntheticid+'.fits'
		# Extracting FITs times and fluxes
		ffits=fits.open(fitslocation)
		alltime=ffits[1].data['TIME']
		allflux=ffits[1].data['PDCSAP_FLUX']
		ffits.close()
		# Setting up figure depending on whether to have 2 or 3 zoom panels
		if len(buffertimes)==3:
			fig,(ax1,ax2,ax3) = plt.subplots(3,1)
			axes=[ax1,ax2,ax3]
			if plot==False:
				fig.set_size_inches(10,15.0)
				fig.subplots_adjust(right=0.98,top=0.98,left=0.1,bottom=0.04)
		if len(buffertimes)==2:
			fig,(ax1,ax2,) = plt.subplots(2,1)
			axes=[ax1,ax2]
			if plot==False:
				fig.set_size_inches(10,10.0)
				fig.subplots_adjust(right=0.98,top=0.97,left=0.1,bottom=0.05)
		#axes[0].set_title('Kepid = '+str(kepid)+'. Synthetic ID = '+syntheticid+'. Quarter = '+quarter+'.')	
		yaxisformatter = ScalarFormatter(useOffset=False)
		for k in range(len(buffertimes)):
		# Cutting out region of transit +/- buffertime
			buffertime=buffertimes[k]
			ax=axes[k]
			time=24.0*(alltime[np.where((alltime>(userxmin-buffertime)) & (alltime<(userxmax+buffertime)))]-userxmid)
			origflux=allflux[np.where((alltime>(userxmin-buffertime)) & (alltime<(userxmax+buffertime)))]
			flux=origflux/np.nanmedian(origflux)  # Normalizing fluxes to ~1.0
			ax.plot(time,flux,marker='o',color='b',ls='-')
			# Plotting extra edge space in x-axis
			xlim=np.nanmax([abs(np.nanmin(time)),abs(np.nanmax(time))])+2.0
			ax.set_xlim(-xlim,xlim)
			# Plotting extra edge space in y-axis
			ydiff=np.nanmax(flux)-np.nanmin(flux)
			ydifffactor=0.05
			ylims=(np.nanmin(flux)-ydiff*ydifffactor,np.nanmax(flux)+ydiff*ydifffactor)
			ax.set_ylim(ylims)
			ax.yaxis.set_major_formatter(yaxisformatter)
			# Drawing vertical lines around transit
			normuserxmin,normuserxmax=24.0*(userxmin-userxmid),24.0*(userxmax-userxmid)
			ax.plot([0,0],ylims,marker='None',color='k',ls='--')
			ax.plot([normuserxmin,normuserxmin],ylims,marker='None',color='k',ls='-')
			ax.plot([normuserxmax,normuserxmax],ylims,marker='None',color='k',ls='-')
			ax.set_ylabel('Normalized flux')
			ax.set_xlabel('Time (hours from midtransit)')
		# Writing figures or plotting
		if writefigures==True:
			fig.savefig('./round2cutouts/'+datatype+'/'+str(len(buffertimes))+'panel/tce'+str(tceid).zfill(5)+'.png')
			plt.close(fig)
		elif plot==True:
			plt.show()
		print '{0:0.2f}'.format(100.0*float(tceid)/float(len(mdwarftces)+len(syntces)))+'% completed. Datatype = '+datatype+'. j = '+str(j)+'. TCE ID = '+str(tceid).zfill(5)+'.'
	if writetcefiles==True:  # Writing out TCE array with TCE IDs
		tces['tceid']=tceidcolumn
		astropyascii.write(tces,'./round2cutouts/'+datatype+'tces.dat')







