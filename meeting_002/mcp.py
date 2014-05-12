import os
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
import scipy.stats as stats
from matplotlib import font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import mcmath
from matplotlib.pylab import date2num,num2date
from mcmath import n2d,d2n
#from mpl_toolkits.basemap import interp as bminterp

global plt,date2num,num2date,n2d,d2n

class NoVarError(Exception):
    def __init__(self,valus):
        self.value = valus
    def __str__(self):
        return repr(self.value)

class MCPlotError(Exception):
	def __init__(self,valus):
		self.value = valus
	def __str__(self):
		return repr(self.value)

##
##### BEGIN plot1

def plot1(x,y=None,clrfirst=True,newfig=False,fylim=None,**kwargs):
	if y is not None:
		ymin = np.nanmin(y) - 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
		ymax = np.nanmax(y) + 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
	else:
		ymin = np.nanmin(x) - 0.02 * (abs(np.nanmax(x) - np.nanmin(x)))
		ymax = np.nanmax(x) + 0.02 * (abs(np.nanmax(x) - np.nanmin(x)))
	if newfig is False:
		fig= plt.gcf()
	else:
		fig = plt.figure()
	if clrfirst is True:
		plt.clf()
		ax = fig.add_subplot(111,ylim=(ymin,ymax))
	else:
		ax = plt.gca()
	if (np.sign(ymax) != np.sign(ymin)):
		ax.axhline(0.0,c='black',lw=1.2,alpha=0.7)
	if y is None:
		y = x.copy()
		x = np.arange(y.size)

	ax.plot(x,y,'+-',**kwargs)

	ax.grid(True)
#	plt.show()
	##ylabel('temperature anomaly')
	plt.figure(fig.number)
	return None

### BEGIN dateplt1

def dateplt1(t,y,add=None,yl=5,fxlim=None,fylim=None,fighold=False,nomth=False,**kwargs):
	from matplotlib.dates import YearLocator,DateFormatter,MonthLocator
	if add is None:
		ymin = np.nanmin(y) - 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
		ymax = np.nanmax(y) + 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
		if isinstance(y,np.ma.MaskedArray):
			data_extent = np.where(y.mask == False)
			if len(data_extent) == 0:
				data_start = 0
				data_end = y.size - 1
			else:
				data_start = data_extent[0][0]
				data_end = data_extent[0][-1]
		else:
			data_start = 0
			data_end = y.size - 1
		if (isinstance(t[0],np.float64) is True) or (isinstance(t[0],np.float32) is True):
			datemin = num2date(t[data_start])
			datemax = num2date(t[data_end])
		elif isinstance(t[0],dt):
			datemin = t[data_start]
			datemax = t[data_end]
		else:
			return None
		datediff = datemax - datemin
		datediff = datediff / 50
		datemin = datemin - datediff
		datemax = datemax + datediff
		if fighold is True:
			plt.hold(False)
			fig = plt.gcf()
			plt.clf()
		else:
			fig = plt.figure()
		if (fylim is not None):
			if type(fylim) is not tuple:
				print "option <fylim> must be a tuple"
				return None
			(ymin,ymax) = fylim
		ax = fig.add_subplot(111,autoscale_on=False,ylim=(ymin,ymax))
		if (np.sign(ymax) != np.sign(ymin)):
			ax.axhline(0.0,c='black',lw=1.2,alpha=0.7)
		line1 = ax.plot(t,y,'+-',**kwargs)
		years = YearLocator(yl)
		yearsFmt = DateFormatter('%Y')
		if nomth == False:
			months = MonthLocator(1)

		ax.set_xlim(datemin,datemax)

		#fig.autofmt_xdate()
		ax.xaxis.set_major_locator(years)
		ax.xaxis.set_major_formatter(yearsFmt)
		if nomth == False:
			ax.xaxis.set_minor_locator(months)

		ax.grid(True)
		my_fmt_xdate(ax)
		##plt.show()
		##ylabel('temperature anomaly')
		plt.figure(fig.number)
		return line1

	else:
		xchange = 0
		ychange = 0
		fig = plt.gcf()
		ax = plt.gca()
		years = YearLocator(yl)
		yearsFmt = DateFormatter('%Y')
		months = MonthLocator(1)
		line1 = ax.plot(t,y,'+-',**kwargs)
		if isinstance(y,np.ma.MaskedArray):
			data_extent = np.where(y.mask == False)
			if len(data_extent) == 0:
				data_start = 0
				data_end = y.size - 1
			else:
				data_start = data_extent[0][0]
				data_end = data_extent[0][-1]
		else:
			data_start = 0
			data_end = y.size - 1
		if (isinstance(t[0],np.float64) is True) or (isinstance(t[0],np.float32) is True):
			datemin = num2date(t[data_start])
			datemax = num2date(t[data_end])
		elif isinstance(t[0],dt):
			datemin = t[data_start]
			datemax = t[data_end]
		else:
			return None

		datediff = datemax - datemin
		datediff = datediff / 50
		datemin = datemin - datediff
		datemax = datemax + datediff
		ymin = np.nanmin(y) - 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
		ymax = np.nanmax(y) + 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
		(oldymin,oldymax) = ax.get_ylim()
		(oldxmin,oldxmax) = ax.get_xlim()

		if (fylim is not None):
			if type(fylim) is not tuple:
				print "option <fylim> must be a tuple"
				return None
			(ymin,ymax) = fylim
		if datemin.toordinal() < oldxmin:
			xchange = 1
			oldxmin = datemin.toordinal()
		if datemax.toordinal() > oldxmax:
			xchange = 1
			oldxmax = datemax.toordinal()
		if ymin < oldymin:
			ychange = 1
			oldymin = ymin
		if ymax > oldymax:
			ychange = 1
			oldymax = ymax
		if xchange:
			plt.setp(ax,xlim=ax.set_xlim(oldxmin,oldxmax))
		if ychange:
			plt.setp(ax,ylim=ax.set_ylim(oldymin,oldymax))
		ax.xaxis.set_major_locator(years)
		ax.grid(True)
		my_fmt_xdate(ax)
		plt.figure(fig.number)
		return line1

### END dateplt1

### BEGIN dateplt_m

def dateplt_m(t,y,yl=5,fxlim=None,fylim=None,ticksize=10,nomth=False,fmt='+-',**kwargs):
	from matplotlib.dates import YearLocator,DateFormatter,MonthLocator
	xchange = 0
	ychange = 0
	fig = plt.gcf()
	ax = plt.gca()
	years = YearLocator(yl)
	yearsFmt = DateFormatter('%Y')
	if nomth == False:
		months = MonthLocator(1)
	
	if isinstance(y,np.ma.MaskedArray):
		data_extent = np.where(y.mask == False)
		if len(data_extent) == 0:
			data_start = 0
			data_end = y.size - 1
		else:
			data_start = data_extent[0][0]
			data_end = data_extent[0][-1]
	else:
		data_start = 0
		data_end = y.size - 1
	if (isinstance(t[0],np.float64) is True) or (isinstance(t[0],np.float32) is True):
		datemin = num2date(t[data_start])
		datemax = num2date(t[data_end])
	elif isinstance(t[0],dt):
		datemin = t[data_start]
		datemax = t[data_end]
	else:
		return None

	datediff = datemax - datemin
	datediff = datediff / 50
	datemin = datemin - datediff
	datemax = datemax + datediff
	ymin = np.nanmin(y) - 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
	ymax = np.nanmax(y) + 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
	(oldymin,oldymax) = ax.get_ylim()
	(oldxmin,oldxmax) = ax.get_xlim()
	

	if (np.sign(ymax) != np.sign(ymin)):
		ax.axhline(0.0,c='black',lw=1.2,alpha=0.7)
	line1 = ax.plot(t,y,fmt,**kwargs)
#	ax.set_xlim(datemin,datemax)
	ax.grid(True)
	

	if (datemin.toordinal() < oldxmin) | (oldxmin == 0.):
		xchange = 1
		oldxmin = datemin.toordinal()
	if datemax.toordinal() > oldxmax:
		xchange = 1
		oldxmax = datemax.toordinal()
	if ymin < oldymin:
		ychange = 1
		oldymin = ymin
	if ymax > oldymax:
		ychange = 1
		oldymax = ymax

	if (fylim is not None):
		if type(fylim) is not tuple:
			print "option <fylim> must be a tuple"
			return None
		(ymin,ymax) = fylim
		plt.setp(ax,ylim=ax.set_ylim(ymin,ymax))
	else:
		if xchange:
			plt.setp(ax,xlim=ax.set_xlim(oldxmin,oldxmax))
		if ychange:
			plt.setp(ax,ylim=ax.set_ylim(oldymin,oldymax))
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	if nomth == False:
		ax.xaxis.set_minor_locator(months)

	##plt.show()

	my_fmt_xdate(ax,rot=45,hal='right')

	xtl = ax.get_xticklabels()
	ytl = ax.get_yticklabels()
	plt.setp(xtl,'size',ticksize)
	plt.setp(ytl,'size',ticksize)

	plt.figure(fig.number)
	return line1

### END dateplt_m

### BEGIN dateplt_m2

def dateplt_m2(t,y,yl=5,fxlim=None,fylim=None,ticksize=10,nomth=False,fmt='+-',**kwargs):
	"""Different from dateplt_m in that it's just a 2-digit year for x tick marks"""
	
	from matplotlib.dates import YearLocator,DateFormatter,MonthLocator
	xchange = 0
	ychange = 0
	fig = plt.gcf()
	ax = plt.gca()
	years = YearLocator(yl)
	yearsFmt = DateFormatter("'%y")
	if nomth == False:
		months = MonthLocator(1)
	
	if isinstance(y,np.ma.MaskedArray):
		data_extent = np.where(y.mask == False)
		if len(data_extent) == 0:
			data_start = 0
			data_end = y.size - 1
		else:
			data_start = data_extent[0][0]
			data_end = data_extent[0][-1]
	else:
		data_start = 0
		data_end = y.size - 1
	if (isinstance(t[0],np.float64) is True) or (isinstance(t[0],np.float32) is True):
		datemin = num2date(t[data_start])
		datemax = num2date(t[data_end])
	elif isinstance(t[0],dt):
		datemin = t[data_start]
		datemax = t[data_end]
	else:
		return None

	datediff = datemax - datemin
	datediff = datediff / 50
	datemin = datemin - datediff
	datemax = datemax + datediff
	ymin = np.nanmin(y) - 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
	ymax = np.nanmax(y) + 0.02 * (abs(np.nanmax(y) - np.nanmin(y)))
	(oldymin,oldymax) = ax.get_ylim()
	(oldxmin,oldxmax) = ax.get_xlim()
	

	if (np.sign(ymax) != np.sign(ymin)):
		ax.axhline(0.0,c='black',lw=1.2,alpha=0.7)
	line1 = ax.plot(t,y,fmt,**kwargs)
#	ax.set_xlim(datemin,datemax)
	ax.grid(True)
	

	if (datemin.toordinal() < oldxmin) | (oldxmin == 0.):
		xchange = 1
		oldxmin = datemin.toordinal()
	if datemax.toordinal() > oldxmax:
		xchange = 1
		oldxmax = datemax.toordinal()
	if ymin < oldymin:
		ychange = 1
		oldymin = ymin
	if ymax > oldymax:
		ychange = 1
		oldymax = ymax

	if (fylim is not None):
		if type(fylim) is not tuple:
			print "option <fylim> must be a tuple"
			return None
		(ymin,ymax) = fylim
		plt.setp(ax,ylim=ax.set_ylim(ymin,ymax))
	else:
		if xchange:
			plt.setp(ax,xlim=ax.set_xlim(oldxmin,oldxmax))
		if ychange:
			plt.setp(ax,ylim=ax.set_ylim(oldymin,oldymax))
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	if nomth == False:
		ax.xaxis.set_minor_locator(months)

	xtl = ax.get_xticklabels()
	ytl = ax.get_yticklabels()
	plt.setp(xtl,'size',ticksize)
	plt.setp(ytl,'size',ticksize)

	##plt.show()

	##my_fmt_xdate(ax,rot=45,hal='right')
	
	plt.figure(fig.number)
	return line1

### END dateplt_m2

##
#### BEGIN longts

def longts(t,x,yl=5,fxlim=None,fylim=None,ticksize=10,nomth=False,**kwargs):
	fig = plt.figure(figsize=(12,5))
	ax = fig.add_axes((0.07,0.1,0.88,0.8),autoscale_on=False,ylim=(-0.01,0.01))
	line1 = dateplt_m(t,x,yl=yl,fxlim=fxlim,fylim=fylim,ticksize=ticksize,nomth=nomth,**kwargs)
	return line1

def barplt1(x,y,yl=1,width=1,mth=7):
	from matplotlib.dates import YearLocator,DateFormatter,MonthLocator
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.bar(x,y,align='center',width=width)
	years = YearLocator(yl,month=1,day=1)
	yearsFmt = DateFormatter('%Y')
	#months = MonthLocator(mth)
	datemin = x[0] + td(-20)
	datemax = x[len(x) - 1] + td(20)
	ax.set_xlim(datemin,datemax)
	ax.grid(True)
	fig.autofmt_xdate()
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	#ax.xaxis.set_minor_locator(months)
	
	##plt.show()
	plt.figure(fig.number)
	return None

def my_fmt_xdate(ax=None,rot=30,hal='right'):
	if ax is None:
		ax = plt.gca()

	for i in ax.get_xmajorticklabels():
		i.set_rotation(rot)
		i.set_ha(hal)

	return None

def my_tax_spacing(yl=5,ax=None):
	from matplotlib.dates import YearLocator,DateFormatter,MonthLocator
	if ax is None:
		ax = plt.gca()
	years = YearLocator(yl,month=1)
	yearsFmt = DateFormatter('%Y')
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	fig = plt.gcf()
	plt.figure(fig.number)
	return None


def latlonlab(lat=None,lon=None):
	retstr = []
	if lat is not None:
		if lat < 0:
			latstr = "%d$^{\circ}$ S" % abs(int(lat))
		elif lat == 0:
			latstr = "Equator"
		else:
			latstr = "%d$^{\circ}$ N" % abs(int(lat))
		retstr.append(latstr)
	if lon is not None:
		if lon < 0:
			lonstr = "%d$^{\circ}$ W" % abs(int(lon))
		elif lon >= 360:
			lonstr = "%d$^{\circ}$ E" % abs(int(lon - 360))
		elif lon > 180:
			lonstr = "%d$^{\circ}$ W" % abs(int(360 - lon))
		else:
			lonstr = "%d$^{\circ}$ E" % abs(int(lon))
		retstr.append(lonstr)

	if len(retstr) > 1:
		return retstr
	elif len(retstr) == 1:
		return retstr[0]
	else:
		return None


def latlontxt(lat=None,lon=None):
	retstr = []
	if lat is not None:
		if lat < 0:
			latstr = "%d S" % abs(int(lat))
		elif lat == 0:
			latstr = "Equator"
		else:
			latstr = "%d N" % abs(int(lat))
		retstr.append(latstr)
	if lon is not None:
		if lon < 0:
			lonstr = "%d W" % abs(int(lon))
		elif lon >= 360:
			lonstr = "%d E" % abs(int(lon - 360))
		elif lon > 180:
			lonstr = "%d W" % abs(int(360 - lon))
		else:
			lonstr = "%d E" % abs(int(lon))
		retstr.append(lonstr)

	if len(retstr) > 1:
		return retstr
	elif len(retstr) == 1:
		return retstr[0]
	else:
		return None



def figtitle(s):
	fig = plt.gcf()
	fig.text(0.5,0.98,s,ha='center',va='top',size='large',weight='bold')
	plt.figure(fig.number)
	return None

####
######## BEGIN shade_coord

def shade_coord(xin,yin,datain=None,lon0=None):
	from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
	from mcmath import edges

	# invert lat coords if they are descending instead of ascending
	if yin[-1] < yin[0]:
		yin = yin[::-1]
		if datain is not None:
			ydimloc = np.where(yin.size == np.array(datain.shape))[0]
			if len(ydimloc) == 0:
				raise MCPlotError("no dimension in 'datain' matches length of 'yin'")
			y_ind = []
			for i in xrange(datain.ndim):
				y_ind.append(slice(None))
			y_ind[ydimloc[0]] = slice(None,None,-1)
			datain = datain[y_ind]
	yedg = edges(yin)
	
	# convert xin to -180:180 range; roll and addcyclic as needed to have lon0 as central lon
	xedg = edges(xin)
	xspan = xedg[-1] - xedg[0]
	if ((xspan < 365.) & (xspan > 355.)):
		xin = np.where(xin >= 180,xin-360,xin)
		roll_ind = np.where(xin < np.roll(xin,1))[0]
		if len(roll_ind) == 0:
			raise MCPlotError("can't find pivot between +180 and -180 in xin")
		xin = np.roll(xin,-roll_ind)
		if datain is not None:
			xdimloc = np.where(xin.size == np.array(datain.shape))[-1]
			if len(xdimloc) == 0:
				raise MCPlotError("no dimension in 'datain' matches length of 'xin'")
			datain = np.roll(datain,-roll_ind,axis=xdimloc)

			if (lon0 is not None):
				(datain,xin) = addcyclic(datain,xin)
				(datain,xin) = shiftgrid(lon0-180,datain,xin)

			xedg = edges(xin)

			## Fix for keeping the squares near the edges of global maps filled just to the edge (map boundary)
			if (lon0 is not None):
				if xedg[0] < (lon0-180):
					xedg[0] = lon0-180
				if xedg[-1] > (lon0+180):
					xedg[-1] = lon0+180

			return [xedg,yedg,datain,xin,yin]
		else:
			xedg = edges(xin)

	return [xedg,yedg,datain,xin,yin]


######## END shade_coord
####


####
######## BEGIN shade

def shade(xin,yin,valsin,proj='cyl',lon0=None,res='c',blat=30,lat0=None,fz=(11.6,5.8),add=None,cbticksize=None,cbticks=None,nocb=None,ex=None,sx=None,xe=None,cllw=1.0,cbnds=None,cbext='both',hatch=None,hatch1=None,hatchmask=None,**kwargs):
	"""An increasingly complicated function for generating decent graphs of shaded contour plots
	on the surface of the earth.  This version is really only designed to work with rectilinear
	grids.  Placement of colorbars in multi-panel figures has proven to be a big problem, and more
	parameters now exist to try to manage this better in multiple circumstances.
	Params: xin - longtitude array of data (1-d)
	        yin - latitude array of data (1-d)
		valsin - data to shade
		proj - projection to employ,
		       currently only the following are enabled: cyl, moll, cass, tmerc, npl, spl, laea, robin
		lon0 - central longitude meridian for plot
		res - same as in Basemap
		blat - bounding latitude for polar projections: npl,spl
		fz - tuple for figure size
		add - add the shade plot to current axes instead of making a new figure
		cbticksize, cbticks - self-explanatory, hopefully
		nocb - don't auto-add a colorbar
		COLORBAR management:
		ex - end of x-interval in which to create a colorbar (colorbar should not step over this (0,1)-axis value)
		sx - spacer for colorbar (fraction of current figure 'x-axis' to draw colorbar from edge of last map

		hatch - use the built-in (new as of 1.2.0) hatching keyword arg; called using the symbol desired
		        e.g., hatch='x'
		hatch1 - use my hatching scheme which makes smaller hatching marks (just dots at the moment)
		         using this argument requires passing a 2-tuple, with the desired density of dots first,
		         and the desired size second; a suggested first attempt is "hatch1=(5,0.3)"
		hatchmask - required for either of the two hatching methods; needs to be a MaskedArray,
		            and not just a mask (not just a boolean array)
			"""
	from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
	from mcmath import edges

	if add is None:
		fig = plt.figure(figsize=fz)
		ax = fig.add_axes((0.05,0.05,0.8,0.9),autoscale_on=False)
	else:
		fig = plt.gcf()
		ax = plt.gca()

	xedg,yedg,dataout,xout,yout = shade_coord(xin,yin,valsin,lon0)
	
	if proj == 'moll':
		m1 = Basemap(projection=proj,lon_0=lon0,resolution=res)
	elif proj == 'robin':
		m1 = Basemap(projection=proj,lon_0=lon0,resolution=res)
##	elif proj == 'eck4':
##		m1 = Basemap(projection=proj,lon_0=xedg[0],resolution=res)
	elif proj == 'cyl':
		m1 = Basemap(projection=proj,llcrnrlat=yout[0],llcrnrlon=xout[0],
					 urcrnrlat=yout[-1],urcrnrlon=xout[-1],resolution=res)
	elif proj == 'cass':
		m1 = Basemap(llcrnrlat=yout[0],llcrnrlon=xout[0],urcrnrlat=yout[-1],
					 urcrnrlon=xout[-1],resolution=res,projection=proj,
					 lon_0=(xedg[0]+xedg[-1])/2,lat_0=(yedg[0]+yedg[-1])/2)
	elif proj == 'tmerc':
		m1 = Basemap(llcrnrlat=yout[0],llcrnrlon=xout[0],urcrnrlat=yout[-1],
					 urcrnrlon=xout[-1],resolution=res,projection=proj,
					 lon_0=(xout[0]+xout[-1])/2,lat_0=(yout[0]+yout[-1])/2)
	elif proj == 'npl':
		m1 = Basemap(projection='nplaea',boundinglat=blat,lon_0=lon0,resolution=res)
	elif proj == 'spl':
		m1 = Basemap(projection='splaea',boundinglat=blat,lon_0=lon0,resolution=res)
	elif proj == 'laea':
		if (lat0 is None):
			raise MCPlotError("central latitude not specified")
		m1 = Basemap(width=wid,height=ht,resolution=res,projection='laea',\
            lat_ts=lat0,lat_0=lat0,lon_0=lon0)
	else:
		m1 = Basemap(projection=proj,llcrnrlat=yout[0],llcrnrlon=xout[0],
					 urcrnrlat=yout[-1],urcrnrlon=xout[-1],resolution=res)


	## FIX: there is a problem with the splaea projection when one of the
	## returned latitudes from shade_coord is +90 deg (North).  Then the
	## Basemap object (here, m1) returns for some longitude map projection
	## coords 1e30.  To fix this, I scrap the last row of yedg if it's +90 (only for splaea).
	if proj == 'spl':
		if yedg[-1] == 90.:
			yedg = yedg[:-1]
			dataout = dataout[:-1]
	
	(x,y) = m1(*np.meshgrid(xedg,yedg))

	m1.drawmapboundary()

	cs = m1.pcolormesh(x,y,dataout,**kwargs)
	
	m1.drawcoastlines(linewidth=cllw)

	plt.draw()

	if not bool(nocb):
		if add is None:
			bbox1 = ax.get_position()
			endx = bbox1.get_points()[1,0]
			endy = bbox1.get_points()[1,1]
			begy = bbox1.get_points()[0,1]
			yext = endy - begy
			xleft = 1 - endx
			xspacer = 0.15*xleft
			xext = 0.83*xleft
			ax2 = fig.add_axes((endx+xspacer,0.1,0.02,0.8),autoscale_on=False)
			if cbticks is None:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,boundaries=cbnds)
			else:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks,boundaries=cbnds)
		else:
			bbox1 = plt.gca().get_position()
			endx = bbox1.get_points()[1,0]
			endy = bbox1.get_points()[1,1]
			begy = bbox1.get_points()[0,1]
			yext = endy - begy
			if ex is not None:
				xleft = ex - endx
			else:
				xleft = 1 - endx
			if sx is not None:
				xspacer = sx
			else:
				xspacer = 0.15*xleft
			if xe is not None:
				xext = xe
			else:
				xext = 0.83*xleft
			if xext > (0.08*yext):
				xext = 0.08*yext

			ax2 = fig.add_axes((endx+xspacer,begy,xext,yext),autoscale_on=False)
			if cbticks is None:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,boundaries=cbnds)
			else:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks,boundaries=cbnds)
		if cbticksize is not None:
			ytl = ax2.get_yticklabels()
			plt.setp(ytl,'size',cbticksize)

	if (hatch is not None) & (hatchmask is not None):
		xedg,yedg,hatchout,xout,yout = shade_coord(xin,yin,hatchmask,lon0)
		(x,y) = m1(*np.meshgrid(xedg,yedg[:-1]))
		m1.contourf(x,y,hatchout,hatches=hatch,alpha=0.0)

	if (hatch1 is not None) & (hatchmask is not None):
		xedg,yedg,hatchout,xout,yout = shade_coord(xin,yin,hatchmask,lon0)
		(x,y) = m1(*np.meshgrid(xedg,yedg))
		density = hatch1[0]
		a2 = hatchout * hatch1[1]
		(x,y) = m1(*np.meshgrid(xout,yout))
		m1.scatter(x[::density,::density],y[::density,::density],a2[::density,::density],ax=ax,marker=',')


	return m1,cs


##	return m1

######## END shade
####



####
######## BEGIN contourf

### WARNING!  Under construction!!
### like shade(), but for making plots using contourf instead of pcolormesh
### plan to add this functionality to shade() once it's finished.


def contourf(xin,yin,valsin,proj='cyl',lon0=None,res='c',blat=30,lat0=None,fz=(11.6,5.8),add=None,cbticksize=None,cbticks=None,nocb=None,ex=None,sx=None,xe=None,cllw=1.0,cbnds=None,cbext='both',hatch=None,hatch1=None,hatchmask=None,**kwargs):
	"""An increasingly complicated function for generating decent graphs of shaded contour plots
	on the surface of the earth.  This version is really only designed to work with rectangular
	grids.  Placement of colorbars in multi-panel figures has proven to be a big problem, and more
	parameters now exist to try to manage this better in multiple circumstances.
	Params: xin - longtitude array of data (1-d)
	        yin - latitude array of data (1-d)
		valsin - data to shade
		proj - projection to employ,
		       currently only the following are enabled: cyl, moll, cass, tmerc, npl, spl, laea, robin
		lon0 - central longitude meridian for plot
		res - same as in Basemap
		blat - bounding latitude for polar projections: npl,spl
		fz - tuple for figure size
		add - add the shade plot to current axes instead of making a new figure
		cbticksize, cbticks - self-explanatory, hopefully
		nocb - don't auto-add a colorbar
		COLORBAR management:
		ex - end of x-interval in which to create a colorbar (colorbar should not step over this (0,1)-axis value)
		sx - spacer for colorbar (fraction of current figure 'x-axis' to draw colorbar from edge of last map

		hatch - use the built-in (new as of 1.2.0) hatching keyword arg; called using the symbol desired
		        e.g., hatch='x'
		hatch1 - use my hatching scheme which makes smaller hatching marks (just dots at the moment)
		         using this argument requires passing a 2-tuple, with the desired density of dots first,
		         and the desired size second; a suggested first attempt is "hatch1=(5,0.3)"
		hatchmask - required for either of the two hatching methods; needs to be a MaskedArray,
		            and not just a mask (not just a boolean array)
			"""
	from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
	from mcmath import edges

	if add is None:
		fig = plt.figure(figsize=fz)
		ax = fig.add_axes((0.05,0.05,0.8,0.9),autoscale_on=False)
	else:
		fig = plt.gcf()
		ax = plt.gca()

	xedg,yedg,dataout,xout,yout = shade_coord(xin,yin,valsin,lon0)
	
	if proj == 'moll':
		m1 = Basemap(projection=proj,lon_0=xedg[0]+180,resolution=res)
	elif proj == 'robin':
		m1 = Basemap(projection=proj,lon_0=xedg[0]+180,resolution=res)
##	elif proj == 'eck4':
##		m1 = Basemap(projection=proj,lon_0=xedg[0]+180,resolution=res)
	elif proj == 'cyl':
		m1 = Basemap(projection=proj,llcrnrlat=yout[0],llcrnrlon=xout[0],
					 urcrnrlat=yout[-1],urcrnrlon=xout[-1],resolution=res)
	elif proj == 'cass':
		m1 = Basemap(llcrnrlat=yout[0],llcrnrlon=xout[0],urcrnrlat=yout[-1],
					 urcrnrlon=xout[-1],resolution=res,projection=proj,
					 lon_0=(xedg[0]+xedg[-1])/2,lat_0=(yedg[0]+yedg[-1])/2)
	elif proj == 'tmerc':
		m1 = Basemap(llcrnrlat=yout[0],llcrnrlon=xout[0],urcrnrlat=yout[-1],
					 urcrnrlon=xout[-1],resolution=res,projection=proj,
					 lon_0=(xout[0]+xout[-1])/2,lat_0=(yout[0]+yout[-1])/2)
	elif proj == 'npl':
		m1 = Basemap(projection='nplaea',boundinglat=blat,lon_0=lon0,resolution=res)
	elif proj == 'spl':
		m1 = Basemap(projection='splaea',boundinglat=blat,lon_0=lon0,resolution=res)
	elif proj == 'laea':
		if (lat0 is None):
			raise MCPlotError("central latitude not specified")
		m1 = Basemap(width=wid,height=ht,resolution=res,projection='laea',\
            lat_ts=lat0,lat_0=lat0,lon_0=lon0)
	else:
		m1 = Basemap(projection=proj,llcrnrlat=yout[0],llcrnrlon=xout[0],
					 urcrnrlat=yout[-1],urcrnrlon=xout[-1],resolution=res)
	(x,y) = m1(*np.meshgrid(xedg,yedg))

##	if (proj=='cyl') | (proj=='npl') | (proj=='spl') | (proj=='laea') | (proj=='tmerc') | (proj=='cass'):
	m1.drawmapboundary()

	cs = m1.contourf(x,y,dataout,**kwargs)

	m1.drawcoastlines(linewidth=cllw)

	plt.draw()

	if not bool(nocb):
		if add is None:
			bbox1 = ax.get_position()
			endx = bbox1.get_points()[1,0]
			endy = bbox1.get_points()[1,1]
			begy = bbox1.get_points()[0,1]
			yext = endy - begy
			xleft = 1 - endx
			xspacer = 0.15*xleft
			xext = 0.83*xleft
			ax2 = fig.add_axes((endx+xspacer,0.1,0.02,0.8),autoscale_on=False)
			if cbticks is None:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,boundaries=cbnds)
			else:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks,boundaries=cbnds)
		else:
			bbox1 = plt.gca().get_position()
			endx = bbox1.get_points()[1,0]
			endy = bbox1.get_points()[1,1]
			begy = bbox1.get_points()[0,1]
			yext = endy - begy
			if ex is not None:
				xleft = ex - endx
			else:
				xleft = 1 - endx
			if sx is not None:
				xspacer = sx
			else:
				xspacer = 0.15*xleft
			if xe is not None:
				xext = xe
			else:
				xext = 0.83*xleft
			if xext > (0.08*yext):
				xext = 0.08*yext

			ax2 = fig.add_axes((endx+xspacer,begy,xext,yext),autoscale_on=False)
			if cbticks is None:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,boundaries=cbnds)
			else:
				if cbnds is None:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks)
				else:
					plt.colorbar(cax=ax2,ax=ax,extend=cbext,ticks=cbticks,boundaries=cbnds)
		if cbticksize is not None:
			ytl = ax2.get_yticklabels()
			plt.setp(ytl,'size',cbticksize)

	if (hatch is not None) & (hatchmask is not None):
		xedg,yedg,hatchout,xout,yout = shade_coord(xin,yin,hatchmask,lon0)
		(x,y) = m1(*np.meshgrid(xedg,yedg[:-1]))
		m1.contourf(x,y,hatchout,hatches=hatch,alpha=0.0)

	if (hatch1 is not None) & (hatchmask is not None):
		xedg,yedg,hatchout,xout,yout = shade_coord(xin,yin,hatchmask,lon0)
		(x,y) = m1(*np.meshgrid(xedg,yedg))
		density = hatch1[0]
		a2 = hatchout * hatch1[1]
		(x,y) = m1(*np.meshgrid(xout,yout))
		m1.scatter(x[::density,::density],y[::density,::density],a2[::density,::density],ax=ax,marker=',')


	return m1,cs


##	return m1

######## END contourf
####

####
######## BEGIN cvshade

def cvshade(clon,clat,valsin,proj='cyl',lon0=215,res='c',blat=30,fz=(11.6,5.8),add=None,cbticksize=12,cbticks=None,nocb=None,vlo=None,vhi=None,lat0=None,wid=12000000,ht=8000000,**kwargs):
	from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
	from mcmath import edges
	from matplotlib.patches import Polygon
	from matplotlib.collections import PatchCollection
	from ferr import use

##	dgrid = use('/scratch/local2/u241180/data/mill/MPI_lon-lat-pressure-points_ocean-grid.nc',silent=1)
##	clon = dgrid.v['grid_center_lon'][:]
##	clat = dgrid.v['grid_center_lat'][:]
##	clon_crnr = dgrid.v['grid_corner_lon'][:]
##	clat_crnr = dgrid.v['grid_corner_lat'][:]

	if clon.shape[0] == 4 and clat.shape[0] == 4:
		clat = np.transpose(clat,(1,2,0))
		clon = np.transpose(clon,(1,2,0))

	clonc = clon.reshape((-1,4))
	clatc = clat.reshape((-1,4))

	clonc1 = np.where(clonc < 35, clonc + 360, clonc)
	clonc2 = np.empty_like(clonc1)
	clatc2 = clatc.copy()
	lhs_x = np.empty((0,4),dtype=np.float32)
	lhs_y = np.empty((0,4),dtype=np.float32)
	lhs_val = np.empty(0,dtype=np.float32)
	val_flat = valsin.copy()
	if isinstance(valsin,np.ma.MaskedArray):
		val_flat[val_flat.mask] = np.nan
	val_flat = val_flat.flatten()

	for ii in xrange(clonc1.shape[0]):
		clonc2[ii] = clonc1[ii]
		if ((clonc1[ii] > 350).any() & (clonc1[ii] < 50).any()):
			clonc2[ii] = np.where(clonc1[ii] < 50, 395., clonc1[ii])
			lhs_x = np.append(lhs_x,np.where(clonc1[ii] > 50,35.,clonc1[ii]).reshape(1,-1),axis=0)
			lhs_y = np.append(lhs_y,clatc2[ii].reshape(1,-1),axis=0)
			lhs_val = np.append(lhs_val,val_flat[ii])

	val2 = np.append(val_flat,lhs_val)
	if isinstance(valsin,np.ma.MaskedArray):
		val2 = np.ma.masked_invalid(val2)
	clonc2 = np.append(clonc2,lhs_x,axis=0)
	clatc2 = np.append(clatc2,lhs_y,axis=0)

	if proj == 'cyl':
		m1 = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
					 llcrnrlon=35,urcrnrlon=395,resolution='c')
	elif proj == 'moll':
		m1 = Basemap(projection='moll',lon_0=lon0)
	elif proj == 'npl':
		m1 = Basemap(projection='nplaea',boundinglat=blat,lon_0=lon0,resolution='l')
	elif proj == 'spl':
		m1 = Basemap(projection='splaea',boundinglat=blat,lon_0=lon0,resolution='l')
	elif proj == 'laea':
		if (lat0 is None):
			raise MCPlotError("central latitude not specified")
		m1 = Basemap(width=wid,height=ht,resolution='l',projection='laea',\
            lat_ts=lat0,lat_0=lat0,lon_0=lon0)
	else:
		raise MCPlotError("no proj %s found!" % proj)

	if add is None:
		fig = plt.figure(figsize=fz)
		ax = fig.add_axes((0.05,0.05,0.8,0.9),autoscale_on=False)
	else:
		fig = plt.gcf()
		ax = plt.gca()

	patches = np.empty(clonc2.shape[0],dtype=object)

	for ii in xrange(clonc2.shape[0]):
		xx1 = clonc2[ii]
		yy1 = clatc2[ii]
		xx,yy = m1(xx1,yy1)
		indx = np.array(zip(xx,yy))
		patches[ii] = Polygon(indx, True)

	p1 = PatchCollection(patches, edgecolors='none', linewidths=0,**kwargs)

	p1.set_array(val2)
	if ((vlo is not None) & (vhi is not None)):
		p1.set_clim(vmin=vlo,vmax=vhi)
	p1.set_antialiased(False)
	m1.drawmapboundary()
##	ax = gca()
	ax.add_collection(p1)
	m1.drawcoastlines()
	
	if not bool(nocb):
		if add is None:
			ax2 = fig.add_axes((0.89,0.1,0.02,0.8),autoscale_on=False)
			if cbticks is None:
				plt.colorbar(p1,cax=ax2,ax=ax,extend='both')
			else:
				plt.colorbar(p1,cax=ax2,ax=ax,extend='both',ticks=cbticks)
		else:
			bbox1 = ax.get_position()
			endx = bbox1.get_points()[1,0]
			endy = bbox1.get_points()[1,1]
			begy = bbox1.get_points()[0,1]
			begx = bbox1.get_points()[0,0]
			yext = 0.9*(endy - begy)
			xleft = endx - begx
			xspacer = 0.005
			yspacer = 0.02
			xext = 0.03*xleft
##			if xext > (0.08*yext):
##				xext = 0.08*yext
			ax2 = fig.add_axes((endx+xspacer,begy+yspacer,xext,yext),autoscale_on=False)
			if cbticks is None:
				plt.colorbar(p1,cax=ax2,ax=ax,extend='both')
			else:
				plt.colorbar(p1,cax=ax2,ax=ax,extend='both',ticks=cbticks)
		if cbticksize is not None:
			ytl = ax2.get_yticklabels()
			plt.setp(ytl,'size',cbticksize)

	return m1,p1


######## END cvshade
####


####
######## BEGIN multiline

def multiline(x,o='v',**kwargs):
	if o == 'v':
		for i in x:
			plt.axvline(i,**kwargs)
	if o == 'h':
		for i in x:
			plt.axhline(i,**kwargs)

	return None

######## END multiline
####


####
######## BEGIN fftplt1

def fftplt1(x,Mch,dtrnd=True,demean=True):
	fs=1

	if demean:
		x = x - x.mean()

	if dtrnd:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_linear, noverlap=Mch/2, Fs=fs)
	else:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_none, noverlap=Mch/2, Fs=fs)

	xdtrnd = pylab.detrend_linear(x)
	xauto = mcmath.acorr_mlab(xdtrnd,2)
	rhoxauto = (xauto[1] + np.sqrt(abs(xauto[2])))/2
	R = mcmath.redspec(rhoxauto,np.arange(Mch/2),Mch)

	P = Pout[0][:R.size]
	F = Pout[1][:R.size]
	Psum = P.sum()
	Rsum = R.sum()
	PRratio = Psum/Rsum
	Rcmp = R*PRratio

	plt.figure()
	plt.plot(F,P)
	plt.plot(F,Rcmp)

	return (F,P,Rcmp)

######## END fftplt1
####

####
######## BEGIN fftplt2

def fftplt2(x,Mch,pval=0.1,dtrnd=True,demean=True,titlestr='',Dt=0.01):
	"""Differs from fftplt1 by returning a plot highlighting sig peaks, with their
	freq axis tick marks notated by their period.  Also allows adjustment of sig level and title string."""
	fs=1
	titl_plc = (1.,1.)

	if demean:
		x = x - x.mean()

	if dtrnd:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_linear, noverlap=Mch/2, Fs=fs)
	else:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_none, noverlap=Mch/2, Fs=fs)

	xdtrnd = pylab.detrend_linear(x)
	xauto = mcmath.acorr_mlab(xdtrnd,2)
	rhoxauto = (xauto[1] + np.sqrt(abs(xauto[2])))/2
	R = mcmath.redspec(rhoxauto,np.arange(Mch/2),Mch)

	P = Pout[0][:R.size]
	F = Pout[1][:R.size]
	Psum = P.sum()
	Rsum = R.sum()
	PRratio = Psum/Rsum
	Rcmp = R*PRratio

	dof = (x.size / (Mch/2.)) * 1.2
	Fval = stats.f.isf(pval,dof,dof)

	tst = P / Rcmp
	pass1 = np.where(tst > Fval)[0]
	maxs = mcmath.extrema_find(P,'max',t=F,dt=Dt)
	max_ind = np.intersect1d(maxs,pass1)
	Fmaxs = F[max_ind]
	Fmaxs2 = np.append(Fmaxs,0.1)
	Fmaxs2.sort()
	Tmaxs = 1/(Fmaxs2)
	Tmaxs = np.round(Tmaxs,2)

	fig = plt.figure(figsize=(7,5))
	ax = fig.add_axes((0.1,0.14,0.88,0.81))
	ax.plot(F,P,)
	ax.plot(F,Rcmp)
	
	ax.set_xticks(Fmaxs2)
	ax.set_xticklabels(Tmaxs)
	my_fmt_xdate(ax,rot=90,hal='center')
	multiline(Fmaxs,c='red',ls='--')

	xtl = ax.get_xticklabels()
	ytl = ax.get_yticklabels()
	plt.setp(xtl,'size',10)
	plt.setp(ytl,'size',9)

	ppct = int((1 - pval)*100)
	titl_str = '%s FFT (chunksize: %d) :: peak CL: %d%%' % (titlestr,Mch,ppct)
	plt.text(titl_plc[0],titl_plc[1],titl_str,ha='right',va='bottom',size=12,transform=ax.transAxes)
	plt.xlabel('Peak period (yrs)',size=11)


	return (F,P,Rcmp)

######## END fftplt2
####

####
######## BEGIN fftplt3

def fftplt3(x,Mch,pval=0.1,dtrnd=True,demean=True,titlestr='',Dt=0.01):
	"""Differs from fftplt3 solely by plotting the power spectral density
	in log-form, 10*log10(P)."""
	fs=1
	titl_plc = (1.,1.)

	if demean:
		x = x - x.mean()

	if dtrnd:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_linear, noverlap=Mch/2, Fs=fs)
	else:
		Pout = plt.psd(x, NFFT=Mch, detrend=pylab.detrend_none, noverlap=Mch/2, Fs=fs)

	xdtrnd = pylab.detrend_linear(x)
	xauto = mcmath.acorr_mlab(xdtrnd,2)
	rhoxauto = (xauto[1] + np.sqrt(abs(xauto[2])))/2
	R = mcmath.redspec(rhoxauto,np.arange(Mch/2),Mch)

	P = Pout[0][:R.size]
	F = Pout[1][:R.size]
	Psum = P.sum()
	Rsum = R.sum()
	PRratio = Psum/Rsum
	Rcmp = R*PRratio

	dof = (x.size / (Mch/2.)) * 1.2
	Fval = stats.f.isf(pval,dof,dof)

	tst = P / Rcmp
	pass1 = np.where(tst > Fval)[0]
	maxs = mcmath.extrema_find(P,'max',t=F,dt=Dt)
	max_ind = np.intersect1d(maxs,pass1)
	Fmaxs = F[max_ind]
	Fmaxs2 = np.append(Fmaxs,0.1)
	Fmaxs2.sort()
	Tmaxs = 1/(Fmaxs2)
	Tmaxs = np.round(Tmaxs,2)

	ax = plt.gca()
	ax.plot(F,10*np.log10(Rcmp))

	ax.set_xticks(Fmaxs2)
	ax.set_xticklabels(Tmaxs)
	my_fmt_xdate(ax,rot=90,hal='center')
	multiline(Fmaxs,c='red',ls='--')

	xtl = ax.get_xticklabels()
	ytl = ax.get_yticklabels()
	plt.setp(xtl,'size',10)
	plt.setp(ytl,'size',9)

	ppct = int((1 - pval)*100)
	titl_str = '%s FFT (chunksize: %d) :: peak CL: %d%%' % (titlestr,Mch,ppct)
	plt.text(titl_plc[0],titl_plc[1],titl_str,ha='right',va='bottom',size=12,transform=ax.transAxes)
	plt.xlabel('Peak period (yrs)',size=11)


	return (F,P,Rcmp)

######## END fftplt3
####

