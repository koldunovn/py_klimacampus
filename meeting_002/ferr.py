import netCDF4 as nc4
import re
import datetime
from datetime import datetime as dt
import numpy as np
import scipy as sp
import os
import scipy.stats as stats
import matplotlib.dates
from netCDF4 import date2num,num2date
import matplotlib as mpl
import mcmath
from mcmath import n2d,d2n
global nc4,re,datetime,d,v,blah,tulu,xulu,yulu,zulu,tufactor,monthnum,np,sp,stats,date2num,mcmath,mpl,dt,num2date
tulu={'hours':'hours','HRS': 'hours','hour':'hours','HOUR': 'hours','days': 'days','DAYS': 'days','day': 'days','DAY': 'days','seconds': 'seconds','SECONDS': 'seconds','second': 'seconds','SECOND': 'seconds','sec':'seconds','secs':'seconds',}
tufactor = {'hours': (24.), 'seconds': (86400.), 'days': 1.}
yulu=['lat','latitude','ecmwfy','y']
xulu=['lon','longitude','ecmwfx','x']
zulu=['depth','zlev','level','lv','height','z']
monthnum={'jan':1,'JAN':1,'feb':2,'FEB':2,'mar':3,'MAR':3,'apr':4,'APR':4,'may':5,'MAY':5,'jun':6,'JUN':6,'jul':7,'JUL':7,'aug':8,'AUG':8,'sep':9,'SEP':9,'oct':10,'OCT':10,'nov':11,'NOV':11,'dec':12,'DEC':12}


class FerrError(Exception):
    def __init__(self,valus):
        self.value = valus
    def __str__(self):
        return repr(self.value)

################
################
#### MEGA CLASS: USE
################
################

class use:

####
######## BEGIN __init__

	def __init__(self,nfn,silent=False,_append=False):
		#global Nio,Ngl,re,d,v,blah,tulu,xulu,yulu,zulu
		import re
		v = {}
		d = {}
		self.cv = {}
		self.vax={}
		if _append is False:
			blah=nc4.Dataset(nfn)
		elif _append is True:
			blah=nc4.Dataset(nfn,mode='a')
		dims=blah.dimensions.keys()
		varis = blah.variables.keys()
		self.name = nfn
		for x in dims:
			if x in varis:
				d[x.lower()]= blah.variables[x]
				if ('_bnd' not in x) and ('edge' not in x):

					if (self.attchk(d[x.lower()],'axis')):
						ax = self.attchk(d[x.lower()],'axis')
						if (getattr(d[x.lower()],ax).upper() == 'X'):
							self.cv[x.lower()] = 'xax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						if (getattr(d[x.lower()],ax).upper() == 'Y'):
							self.cv[x.lower()] = 'yax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						if (getattr(d[x.lower()],ax).upper() == 'Z'):
							self.cv[x.lower()] = 'zax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						if (getattr(d[x.lower()],ax).upper() == 'T'):
							self.cv[x.lower()] = 'tax'


					if ((hasattr(d[x.lower()],'units')) and (mpl.is_string_like(d[x.lower()].units))):
						if (('east' in d[x.lower()].units.lower()) and ('degree' in d[x.lower()].units.lower())):
							self.cv[x.lower()] = 'xax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						if (('north' in d[x.lower()].units.lower()) and ('degree' in d[x.lower()].units.lower())):
							self.cv[x.lower()] = 'yax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						if (('meter' in d[x.lower()].units.lower()) or ('m' is d[x.lower()].units.lower())):
							self.cv[x.lower()] = 'zax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if hasattr(getattr(d[x.lower()],att),'lower'):
									self.cv[xx][att] = getattr(d[x.lower()],att).lower()
								else:
									self.cv[xx][att] = getattr(d[x.lower()],att)
							continue
						
						units=d[x.lower()].units.split()
						if (units[0].lower() in tulu.keys()):
							if (x.lower() not in self.cv.keys()):
								self.cv[x.lower()] = 'tax'
							xx = x.lower() + '_atts'
							self.cv[xx] = {}
							self.cv[xx]['tunits']=tulu[units[0].lower()]
							self.cv[xx]['tufac']=tufactor[tulu[units[0].lower()]]
							if hasattr(d[x.lower()],'bounds'):
								self.cv[xx]['bounds'] = getattr(d[x.lower()],'bounds').lower()
							elif hasattr(d[x.lower()],'edges'):
								self.cv[xx]['edges'] = getattr(d[x.lower()],'edges').lower()
							r1=re.compile(r'([\d]{1,4}).([\d]{1,2}).([\d]{1,2})')
							date1=r1.search(d[x.lower()].units)
							if date1 is not None:
								if int(date1.group(1)) == 0:
									torg= date1.group(1) + '-' + date1.group(2) + '-' +  date1.group(3)
									self.cv[xx]['torg'] = torg
								else:
									torg=[date1.group(1),date1.group(2),date1.group(3)]
									r1=re.compile(r'([\d]{2}):([\d]{2})[:]{0,1}([\d]{0,2})')
									time1=r1.search(d[x.lower()].units)
									if time1 is not None:
										torg.extend([time1.group(1),time1.group(2)])
										self.cv[xx]['torg'] = dt(int(torg[0]),int(torg[1]),
													int(torg[2]),int(torg[3]),int(torg[4]))
									else:
										self.cv[xx]['torg'] = dt(int(torg[0]),int(torg[1]),int(torg[2]))
							else:
								self.cv[xx]['torg'] = 'climatology'
							atts = d[x.lower()].ncattrs()
							for att in atts:
								if att not in self.cv[xx].keys():
									if hasattr(getattr(d[x.lower()],att),'lower'):
										self.cv[xx][att] = getattr(d[x.lower()],att).lower()
									else:
										self.cv[xx][att] = getattr(d[x.lower()],att)

					if ((x.lower() in xulu) or ('lon' in x.lower())):
						self.cv[x.lower()] = 'xax'
						xx = x.lower() + '_atts'
						self.cv[xx] = {}
						atts = d[x.lower()].ncattrs()
						for att in atts:
							if hasattr(getattr(d[x.lower()],att),'lower'):
								self.cv[xx][att] = getattr(d[x.lower()],att).lower()
							else:
								self.cv[xx][att] = getattr(d[x.lower()],att)
						continue
					if ((x.lower() in yulu) or ('lat' in x.lower())):
						self.cv[x.lower()] = 'yax'
						xx = x.lower() + '_atts'
						self.cv[xx] = {}
						atts = d[x.lower()].ncattrs()
						for att in atts:
							if hasattr(getattr(d[x.lower()],att),'lower'):
								self.cv[xx][att] = getattr(d[x.lower()],att).lower()
							else:
								self.cv[xx][att] = getattr(d[x.lower()],att)
						continue
					if ((x.lower() in zulu) or ('depth' in x.lower())
						or ('lv' in x.lower()) or ('level' in x.lower())):
						self.cv[x] = 'zax'
						xx = x.lower() + '_atts'
						self.cv[xx] = {}
						atts = d[x.lower()].ncattrs()
						for att in atts:
							if hasattr(getattr(d[x.lower()],att),'lower'):
								self.cv[xx][att] = getattr(d[x.lower()],att).lower()
							else:
								self.cv[xx][att] = getattr(d[x.lower()],att)
						continue
					if ('time' in x.lower()):
						if (x.lower() not in self.cv.keys()):
							self.cv[x.lower()] = 'tax'

							# old method for torg attr: self.vax[x]['torg'] = torg
				
		for y in varis:
			if ('_bnd' in y):
				d[y.lower()]= blah.variables[y]
				bnds_re = re.compile(r'(.*)_bnd[s]{0,1}$')
				bres = bnds_re.search(y)
				if bres is not None:
					if len(bres.groups()) > 0:
						if len(bres.group(1)) > 0:
							name = bres.group(1)
							xx = name.lower() + '_atts'
							if xx not in self.cv.keys():
								self.cv[xx] = {}
							if 'bounds' not in self.cv[xx].keys():
								self.cv[xx]['bounds'] = y.lower()
			elif ('edge' in y):
				d[y.lower()]= blah.variables[y]
				edg_re = re.compile(r'(.*)edge[s]{0,1}$')
				eres = edg_re.search(y)
				if eres is not None:
					if len(eres.groups()) > 0:
						if len(eres.group(1)) > 0:
							name = eres.group(1)
							xx = name.lower() + '_atts'
							if xx not in self.cv.keys():
								self.cv[xx] = {}
							if 'edges' not in self.cv[xx].keys():
								self.cv[xx]['edges'] = y.lower()
			elif y not in dims:
				v[y.lower()]=blah.variables[y]
				self.vax[y.lower()] = {}
				v[y.lower()].set_auto_maskandscale(1)
		self.d = d
		self.v = v
		self.f = blah

		#try to get x,y,z,t dims (and dim order)
		for x in self.v.keys():
			vdims=self.v[x].dimensions
			self.vax[x.lower()]['dimord'] = []
			for y in vdims:
				if y.lower() in self.cv.keys():
					self.vax[x.lower()][self.cv[y.lower()]] = y.lower()
					self.vax[x.lower()]['dimord'].append(self.cv[y.lower()][0])
##					yy = y.lower() + '_atts'
##					if yy in self.cv.keys():
##						for kk in self.cv[yy].keys():
##							if ('bound' not in kk) and ('edge' not in kk):
##								self.vax[x.lower()][kk] = self.cv[yy][kk]
		#end of x,y,t,z axis-finding
		if bool(silent) is False:
			self.shda()
		return None

######## END __init__
####


	def shda(self):
		for x in self.v.keys():
			shp=self.v[x].shape
			dim=self.v[x].dimensions
			att=self.v[x].ncattrs()
			print "\n%s: " % (x.upper()),
			for i in range(len(shp)):
				print dim[i],"(%s" % shp[i],
				for y in self.vax[x].keys():
					if self.vax[x][y] == dim[i]:
						print ", %s" % (y),
						continue
				print ") : ",
			print "\n"
			for a in self.v[x].ncattrs():
				print "%s: %s :: " % (a,getattr(self.v[x],a)),
			print "\n"
		return None

	def shax(self):
		for x in self.d.keys():
			if 'bnd' in x or 'edge' in x:
				continue
			vals = self.d[x][:]
			print "\n%s:  [%5.2f %5.2f ... %5.2f %5.2f]" % (x,vals[0],vals[1],vals[-2],vals[-1])
			for a in self.d[x].ncattrs():
				b = getattr(self.d[x],a)
				print "%s: %s :: " % (a,b),
			print "\n--------------\n"
		return None

	def __repr__(self):
		s = "Ferr Dataset obj: " + self.name
		return s
		

	def varchk(self,varcall):
		global vmatch
		vmatch = None
		varcall = "^" + varcall + "$"
		var_re = re.compile(varcall,re.I)
		for vnames in self.v.keys():
			vmatch = var_re.match(vnames)
			if vmatch is not None:
				break
		if vmatch is None:
			for vnames in self.d.keys():
				vmatch = var_re.match(vnames)
				if vmatch is not None:
					break
		if vmatch is not None:
			varcall = vmatch.group()
			return varcall
		else:
			return None

	def attchk(self,vobj,attcall):
		global amatch
		amatch = None
		attcall = "^" + attcall + "$"
		att_re = re.compile(attcall,re.I)
		for anames in vobj.ncattrs():
			amatch = att_re.match(anames)
			if amatch is not None:
				break
		if amatch is not None:
			attcall = amatch.group()
			return attcall
		else:
			return None

##	def bndchk(self,vari,axstr='All'):
##		varcall = self.varchk(vari)

##		if varcall is None:
##			print "\nno variable called %s\n" % (vari)
##			return None

##		var1 = self.v[varcall]
##		var1ax = self.vax[varcall]

##		if 'All' in axstr:
##			all = {}
##			if 'xax' in var1ax.keys():
##				xax = self.d[var1ax['xax']][:]
##				all['x'] = xax

##			if 'yax' in var1ax.keys():
##				yax = self.d[var1ax['yax']][:]
##				all['y'] = yax

##			if 'zax' in var1ax.keys():
##				zax = self.d[var1ax['zax']][:]
##				all['z'] = zax

##			if 'tax' in var1ax.keys():
##				tax = self.d[var1ax['tax']][:]
##				all['t'] = tax

##			return all

##		if 'x' in axstr:
##			if 'xax' in var1ax.keys():
##				xax = self.d[var1ax['xax']][:]
##				return xax

##		if 'y' in axstr:
##			if 'yax' in var1ax.keys():
##				yax = self.d[var1ax['yax']][:]
##				return yax

##		if 'z' in axstr:
##			if 'zax' in var1ax.keys():
##				zax = self.d[var1ax['zax']][:]
##				return zax

##		if 't' in axstr:
##			if 'tax' in var1ax.keys():
##				tax = self.d[var1ax['tax']][:]
##				return tax
	
##		return None



	
##################################
## getax
##################################

	def getax(self,vari,axstr='All'):
		"""Get axis subroutine; gets 'All' by default.
		'All' gets all axes of the variable in question, in a dict,
		with 'x'-> vals structure (where 'x' is the	coordinate variable
		paired with its values).

		Otherwise, can ask for a single axis by putting one of 'x', 'y', 'z', 't' in axstr."""

		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None

		var1 = self.v[varcall]
		var1ax = self.vax[varcall]

		if 'All' in axstr:
			all = {}
			if 'xax' in var1ax.keys():
				xax = self.d[var1ax['xax']][:]
				all['x'] = xax

			if 'yax' in var1ax.keys():
				yax = self.d[var1ax['yax']][:]
				all['y'] = yax

			if 'zax' in var1ax.keys():
				zax = self.d[var1ax['zax']][:]
				all['z'] = zax

			if 'tax' in var1ax.keys():
				tax = self.d[var1ax['tax']][:]
				all['t'] = tax

			return all
			

		if 'x' in axstr:
			if 'xax' in var1ax.keys():
				xax = self.d[var1ax['xax']][:]
				return xax

		if 'y' in axstr:
			if 'yax' in var1ax.keys():
				yax = self.d[var1ax['yax']][:]
				return yax

		if 'z' in axstr:
			if 'zax' in var1ax.keys():
				zax = self.d[var1ax['zax']][:]
				return zax

		if 't' in axstr:
			if 'tax' in var1ax.keys():
				tax = self.d[var1ax['tax']][:]
				return tax
	
		return None


#### end of getax


##################################
## get_tunits
##################################

	def get_tunits(self,vari):
		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None
		
		var1ax = self.vax[varcall]
		if 'tax' not in var1ax.keys():
			print "variable %s has no time axis found by ferr\n" % (varcall)
			return None
		else:
			tax = var1ax['tax']
		tax_atts = tax + '_atts'
		if tax_atts not in self.cv.keys():
			print "time axis %s does not have any attrs assigned to cv dict\n" % tax
			return None
		if 'torg' not in self.cv[tax_atts].keys():
			print "time axis for %s does not have a time origin attribute\n" % (varcall)
			return None
		else:
			torg = self.cv[tax_atts]['torg']
		if 'tunits' not in self.cv[tax_atts].keys():
			print "time axis for %s does not have a time units attribute\n" % (varcall)
			return None
		else:
			tunits1 = self.cv[tax_atts]['tunits']
		tunits = tunits1 + ' since ' + str(torg)

		return tunits

#### end of get_tunits


##################################
## dt_vals
##################################

	def dt_vals(self,vari):

		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None

		if varcall in self.vax.keys():
			var1ax = self.vax[varcall]
			if 'tax' not in var1ax.keys():
				print "variable %s has no time axis found by ferr\n" % (varcall)
				return None
			else:
				tax = var1ax['tax']
				tatts = tax + '_atts'
			tvals1 = self.d[tax][:]
			tunits_str = self.get_tunits(varcall)
			if 'calendar' in self.cv[tatts].keys():
				cal_str = self.cv[tatts]['calendar']
			else:
				cal_str = 'standard'
		if varcall in self.cv.keys():
			var_atts = varcall + '_atts'
			if var_atts not in self.cv.keys():
				print "info for coord_var %s not found\n" % (varcall)
				return None
			tvals1 = self.d[varcall][:]
			if 'units' not in self.cv[var_atts]:
				print "tunits info for coord_var %s not found\n" % (varcall)
				return None
			tunits_str = self.cv[var_atts]['units']
			if 'calendar' in self.cv[var_atts].keys():
				cal_str = self.cv[var_atts]['calendar']
			else:
				cal_str = 'standard'
		if int(tunits_str.split()[2].split('-')[0]) == 0:
			tunits_str = tunits_str.split()[0] + ' since 0001-01-01'
			print "NOTE: the time axis is a climatological time axis\n"

		d0001 = n2d(tvals1,units=tunits_str,cal=cal_str)

## t0001 = date2num(torg) + tvals1 / tufac  ## these lines don't handle torgs of 0001-1-1 quite right
## d0001 = np.array(num2date(t0001))        ## whereas the netCDF4 num2date with proper units str does!

		return d0001


#### end of dt_vals


##################################
## axedges
##################################

	def axedges(self,vari,xyzt):
		"""Calculates axis edges for datapoints in variable 'varcall'.  Looks for axes corresponding
		to XYZT depending on contents of string 'xyzt'.  If xyzt = "xyt", then this function looks for
		an x-, y-, and t-axis for variable 'varcall'.  Failure to find one causes an exception.
		
		Returns a tuple of arrays with edges in order of XYZT for all axes requested."""

		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None

		var1 = self.v[varcall]
		var1ax = self.vax[varcall]
		xflag = 0
		yflag = 0
		zflag = 0
		tflag = 0
		if ('x' in xyzt.lower()):
			xflag = 1
			if ('xax' in var1ax.keys()):
				xname = var1ax['xax']
				xvals = self.d[xname][:]
			else:
				print "axedges error: No x-axis for %s\n"  % (varcall)
				return None
		if ('y' in xyzt.lower()):
			yflag = 1
			if ('yax' in var1ax.keys()):
				yname = var1ax['yax']
				yvals = self.d[yname][:]
				return None
			else:
				print "axedges error: No y-axis for %s\n"  % (varcall)
				return None
		if ('z' in xyzt.lower()):
			zflag = 1
			if ('zax' in var1ax.keys()):
				zname = var1ax['zax']
				zvals = self.d[zname][:]
			else:
				print "axedges error: No z-axis for %s\n"  % (varcall)
				return None
		if ('t' in xyzt.lower()):
			tflag = 1
			if ('tax' in var1ax.keys()):
				tname = var1ax['tax']
				tbndsname = None
				bnds_re = re.compile(r'(.*_bnds)')
				for tr in self.d.keys():
					if bnds_re.match(tr) != None:
						tbndsname = bnds_re.match(tr).group(1)
						tbnds1 = self.d[tbndsname][:]
					else:
						tvals = self.d[tname][:]
			else:
				print "axedges error: No t-axis for %s\n"  % (varcall)
				return None





		return None


#### end of axedges



##################################
## tax_intvl
##################################

	def tax_intvl(self,vari,timeloc,tvalsin=None):
		tspot = ''
		tbegin = ''
		tend = ''
		tstr_err = "Date string format not understood (e.g., '1-jan-1970:31-dec-1990')"
		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None

		
		var1ax = self.vax[varcall]
		if 'tax' not in var1ax.keys():
			print "variable %s has no time axies found by ferr\n" % (varcall)
			return None
		else:
			tax = var1ax['tax']
			tax_atts = tax + '_atts'
		
		d1900 = datetime.datetime(1900,1,1)
		torg = self.cv[tax_atts]['torg']
		ddiff = torg - d1900
		tufac = self.cv[tax_atts]['tufac']
		if tvalsin is None:
			tvals = self.dt_vals(varcall)
		else:
			tvals = tvalsin
			if not isinstance(tvals[0],dt):
				tunits_str = self.get_tunits(varcall)
				dt0 = nc4.num2date(tvals1,units=tunits_str)
				tvals = dt0
		
		date_re = re.compile(r'([\d]{1,2})-([a-zA-Z]{3})-([\d]{4})')
		tbegin = timeloc.rsplit(':')[0]
		tend = timeloc.rsplit(':')[1]
		if date_re.match(tbegin) != None:
			#put date conversion to index value stuff here
			tbeg_dt = datetime.datetime.strptime(tbegin,'%d-%b-%Y')
			tindex = np.argwhere(tvals >= tbeg_dt)[0]
			tbegin = int(tindex)
		else:
			print tstr_err
			return None
		
		if date_re.match(tend) != None:
			#put date conversion to index value stuff here
			tend_dt = datetime.datetime.strptime(tend,'%d-%b-%Y')
			tindex = max(np.argwhere(tvals <= tend_dt))
			tend = int(tindex)
		else:
			print tstr_err
			return None

		return (tbegin,tend)


#### end of tax_intvl

##################################
## gv
##################################

	def gv(self,vari,locstr='',verbose=None):
		"""replacement get_value function with time converted to
		user-friendly ferret-like values"""
		global tloc,zloc,yloc,xloc
		tloc = None
		zloc = None
		yloc = None
		xloc = None
		dtax = None
		
		tspot = ''
		tbegin = ''
		tend = ''
		xspot = 0.
		yspot = 0.
		zspot = 0.
		notax = False
		noxax = False
		noyax = False
		nozax = False
		tstr_err = 'Time string format following \"t=\" not understood'
		xstr_err = 'Longitude string format following \"x=\" not understood'
		ystr_err = 'Latitude string format following \"y=\" not understood'
		zstr_err = 'Depth / height string format following \"z=\" not understood'
		varcall = self.varchk(vari)

		if varcall is None:
			print "\nno variable called %s\n" % (vari)
			return None
		
		var1 = self.v[varcall]
		var1ax = self.vax[varcall]

		axs = self.getax(varcall)
		
		if 't' not in axs.keys():
			notax = 1
		else:
			dtax = self.dt_vals(varcall)
			if dtax is None:
				return None
		if 'x' not in axs.keys():
			noxax = 1
		if 'y' not in axs.keys():
			noyax = 1
		if 'z' not in axs.keys():
			nozax = 1

		t_re = re.compile(r't=([a-zA-Z0-9\-\:]+)')
		x_re = re.compile(r'x=([a-zA-Z0-9\-\:]+)')
		y_re = re.compile(r'y=([a-zA-Z0-9\-\:]+)')
		z_re = re.compile(r'z=([a-zA-Z0-9\-\:]+)')

		if t_re.search(locstr) != None:
			if notax:
				print "Error: no time axis for variable \'varcall\'"
				return None
			tvals = self.d[var1ax['tax']][:]
			timeloc = t_re.search(locstr).group(1)
			#print timeloc
			date_re = re.compile(r'([\d]{1,2})-([a-zA-Z]{3})-([\d]{4})')
			index_re = re.compile(r'([\d]+)')
			if timeloc.find(':') > -1:
				tbegin = timeloc.rsplit(':')[0]
				tend = timeloc.rsplit(':')[1]
				t1 = None
				t2 = None
				if date_re.match(tbegin) != None:
					tbeg_dt = datetime.datetime.strptime(tbegin,'%d-%b-%Y')
					t1 = dtax >= tbeg_dt
				elif index_re.search(tbegin) != None:
					tbegin = int(index_re.search(tbegin).group(1))
				else:
					print tstr_err
					return None

				if date_re.match(tend) != None:
					tend_dt = datetime.datetime.strptime(tend,'%d-%b-%Y')
					t2 = dtax <= tend_dt
				elif index_re.search(tend) != None:
					tend = int(index_re.search(tend).group(1))
				else:
					print tstr_err
					return None

				if (t1 is not None) and (t2 is not None):
					tloc = t1 & t2
				elif type(tbegin) is int and type(tend) is int:
					if tend > tbegin:
						if (tend - tbegin == 1):
							tloc = '%d' % (tbegin)
						else:
							tloc = '%d:%d' % (tbegin,tend)
					else:
						if (tbegin - tend == 1):
							tloc = '%d' % (tend)
						else:
							tloc = '%d:%d' % (tend,tbegin)
				#print 'tloc is ', tloc

			else:
				if date_re.match(timeloc) != None:
					#put axedges stuff in here to properly "place" the tspot_dt in index space
					tspot_dt = datetime.datetime.strptime(timeloc,'%d-%b-%Y')
					tdiff = abs(dtax - tspot_dt)
					tloc = np.zeros(dtax.shape)
					tloc = tloc.astype(bool)
					tloc[np.where(tdiff == np.min(tdiff))[0]] = True
					
				elif index_re.search(timeloc) != None:
					tloc = int((index_re.search(timeloc)).group(1))
				else:
					print tstr_err
					return None

				#print 'tloc is ', tloc


		if x_re.search(locstr) != None:
			if noxax:
				print "Error: no longitude (x) axis for variable \'varcall\'"
				return None

			xax = axs['x']
			xax_360 = np.zeros_like(xax)
			if np.where(xax > 360.)[0].size > 0:
				xax_360 = np.where(xax < 360, xax, xax - 360)
			if np.where(xax < 0.)[0].size > 0:
				xax_360 = np.where(xax < 360, xax, xax - 360)
			if xax_360.nonzero()[0].size == 0:
				xax_360 = xax
			lonloc = x_re.search(locstr).group(1)
			#print lonloc
			E_re = re.compile(r'([\d]+[\.]{0,1}[\d]*)[eE]')
			W_re = re.compile(r'([\d]+[\.]{0,1}[\d]*)[wW]')
			num_re = re.compile(r'([\d]+[\.]{0,1}[\d]*)$')
			index_re = re.compile(r'[\d]+\:{0,1}[\d]*$')
##			if index_re.match(lonloc) != None:
##				pass # this is still broken
			if lonloc.find(':') > -1:
				xbegin = lonloc.rsplit(':')[0]
				xend = lonloc.rsplit(':')[1]
				if E_re.search(xbegin) != None:
					xbegin = np.double(E_re.search(xbegin).group(1))
				elif W_re.search(xbegin) != None:
					xbegin = np.double(W_re.search(xbegin).group(1))
					xbegin = 360.0 - xbegin
				elif num_re.match(xbegin) != None:
					if xbegin >= 360:
						xbegin = np.double(xbegin) - 360
					else:
						xbegin = np.double(xbegin)
				else:
					print xstr_err
					return None

				if E_re.search(xend) != None:
					xend = np.double(E_re.search(xend).group(1))
				elif W_re.search(xend) != None:
					xend = np.double(W_re.search(xend).group(1))
					xend = 360.0 - xend
				elif num_re.match(xend) != None:
					if xend >= 360:
						xend = np.double(xend) - 360
					else:
						xend = np.double(xend)
				else:
					print xstr_err
					return None

				if xend < xbegin:
					xax_mod = np.where(xax_360 < xbegin, xax_360 + 360, xax_360)
					xend = xend + 360
					x1 = xax_mod >= xbegin
					x2 = xax_mod <= xend
					xloc = x1 & x2
				elif xbegin < xend:
					x1 = xax_360 >= xbegin
					x2 = xax_360 <= xend
					xloc = x1 & x2

				#print 'xloc is ', xloc

			else:
				if E_re.search(lonloc) != None:
					xspot = np.double(E_re.search(lonloc).group(1))
				elif W_re.search(lonloc) != None:
					xspot = np.double(W_re.search(lonloc).group(1))
					xspot = 360.0 - xspot
				elif num_re.match(lonloc) != None:
					xspot = np.double(lonloc)
					if xspot >= 360:
						xspot = xspot - 360
				else:
				    print xstr_err
				    return None
				xdiff = abs(xax - xspot)
				xloc = np.zeros(xax.shape).astype(bool)
				xloc[np.argmin(abs(xax_360 - xspot))] = True

				#print 'xloc is ', xloc



		if y_re.search(locstr) != None:
			if noyax:
				print "Error: no latitude axis for variable \'varcall\'"
				return None

			yax = axs['y']
			latloc = y_re.search(locstr).group(1)
			#print latloc  #debugging line, usu. keep commented
			N_re = re.compile(r'([\d]+[\.]{0,1}[\d]*)[nN]')
			S_re = re.compile(r'([\d]+[\.]{0,1}[\d]*)[sS]')
			num_re = re.compile(r'([-]{0,1}[\d]+[\.]{0,1}[\d]*)$')
			index_re = re.compile(r'i[\d]+\:{0,1}[\d]*$')
##			if index_re.match(latloc) != None:
##				pass # this is still broken
			if latloc.find(':') > -1:
				ybegin = latloc.rsplit(':')[0]
				yend = latloc.rsplit(':')[1]
				if N_re.search(ybegin) != None:
					ybegin = np.double(N_re.search(ybegin).group(1))
				elif S_re.search(ybegin) != None:
					ybegin = np.double(S_re.search(ybegin).group(1))
					ybegin = 0.0 - ybegin
				elif num_re.match(ybegin) != None:
					ybegin = np.double(ybegin)
				else:
					print ystr_err
					return None

				if N_re.search(yend) != None:
					yend = np.double(N_re.search(yend).group(1))
				elif S_re.search(yend) != None:
					yend = np.double(S_re.search(yend).group(1))
					yend = 0.0 - yend
				elif num_re.match(yend) != None:
					yend = np.double(yend)
				else:
					print ystr_err
					return None
				
				y1 = yax >= ybegin
				y2 = yax <= yend
				yloc = y1 & y2
				#print 'yloc is ', yloc

			else:
				if N_re.search(latloc) != None:
					yspot = np.double(N_re.search(latloc).group(1))
				elif S_re.search(latloc) != None:
					yspot = np.double(S_re.search(latloc).group(1))
					yspot = 0.0 - yspot
				elif num_re.match(latloc) != None:
					yspot = np.double(latloc)
				else:
					print ystr_err
					return None

				ydiff = abs(yax - yspot)
				yloc = np.zeros(yax.shape).astype(bool)
				yloc[np.argmin(abs(yax - yspot))] = True
				
				#print 'yloc is ', yloc

		if z_re.search(locstr) != None:
			if nozax:
				print "Error: no depth / height axis for variable \'varcall\'"
				return None

			zax = axs['z']
			deploc = z_re.search(locstr).group(1)
			#print latloc  #debugging line, usu. keep commented
			num_re = re.compile(r'([-]{0,1}[\d]+[\.]{0,1}[\d]*)$')
			index_re = re.compile(r'i[\d]+\:{0,1}[\d]*$')
##			if index_re.match(latloc) != None:
##				pass # this is still broken
			if deploc.find(':') > -1:
				zbegin = deploc.rsplit(':')[0]
				zend = deploc.rsplit(':')[1]
				if num_re.match(zbegin) != None:
					zbegin = np.double(zbegin)
				else:
					print zstr_err
					return None

				if num_re.match(zend) != None:
					zend = np.double(zend)
				else:
					print zstr_err
					return None
				
				z1 = zax >= zbegin
				z2 = zax <= zend
				zloc = z1 & z2
				#print 'zloc is ', zloc

			else:
				if num_re.match(deploc) != None:
					zspot = np.double(deploc)
				else:
					print zstr_err
					return None

				zdiff = abs(zax - zspot)
				zloc = np.zeros(zax.shape).astype(bool)
				zloc[np.argmin(abs(zax - zspot))] = True

				#print 'zloc is ', zloc

		def whichloc(dimstr):
			global tloc,zloc,yloc,xloc
			if dimstr == 't':
				if tloc is not None:
					return tloc
				else:
					tloc = slice(0,axs['t'].size)
					return tloc
			if dimstr == 'z':
				if zloc is not None:
					return zloc
				else:
					zloc = slice(0,axs['z'].size)
					return zloc
			if dimstr == 'y':
				if yloc is not None:
					return yloc
				else:
					yloc = slice(0,axs['y'].size)
					return yloc
			if dimstr == 'x':
				if xloc is not None:
					return xloc
				else:
					xloc = slice(0,axs['x'].size)
					return xloc
			return None

		loclist = []

		for i in var1ax['dimord']:
			loclist.append(whichloc(i))

		loctup = tuple(loclist)
		
		outdata = np.squeeze(var1[loctup])

		if isinstance(outdata,np.ma.MaskedArray):
			if outdata.mask.all() == False:
				for i in var1ax['dimord']:
					if (i is 't') & (dtax is not None):
						setattr(outdata,i,dtax[whichloc(i)])
					else:
						setattr(outdata,i,axs[i][whichloc(i)])
				return outdata
			else:
				return None
		elif isinstance(outdata,np.ndarray):
			if outdata.size > 0:
				outdata = np.ma.MaskedArray(outdata)
				for i in var1ax['dimord']:
					if (i is 't') & (dtax is not None):
						setattr(outdata,i,dtax[whichloc(i)])
					else:
						setattr(outdata,i,axs[i][whichloc(i)])
				return outdata
			else:
				return None
		else:
			return None


#### end of gv


################
################
################
##  END OF CLASS 'use'
################
################
################

####################################
##  START OF NETCDF SAVE FUNCTIONS
####################################

##################################
## save
##################################

def save(nfn,data,dataname,x=None,y=None,z=None,t=None,tunits=None,dtim=None,
			 xbnds=None,ybnds=None,zbnds=None,tbnds=None,dims=None,xname='Lon',
			 yname='Lat',zname='Depth',tname='Time',dhist=None,dunits=None,t_atts=None,append=False,app_in_t=False,silent=False):

	if dims is None:
		raise FerrError("Dimensions and order not specified!")
	ndims = data.ndim
	ct_dims = 0
	if x is not None:
		ct_dims = ct_dims + 1
		if not isinstance(x,np.ndarray):
			raise FerrError("x is not an ndarray; all coord var data must be in np.ndarrays")
	if y is not None:
		ct_dims = ct_dims + 1
		if not isinstance(y,np.ndarray):
			raise FerrError("y is not an ndarray; all coord var data must be in np.ndarrays")
	if z is not None:
		ct_dims = ct_dims + 1
		if not isinstance(z,np.ndarray):
			raise FerrError("z is not an ndarray; all coord var data must be in np.ndarrays")
	if (t is not None) or (dtim is not None):
		ct_dims = ct_dims + 1
		if t is not None:
			if not isinstance(t,np.ndarray):
				raise FerrError("t is not an ndarray; all coord var data must be in np.ndarrays")
		if dtim is not None:
			if not isinstance(dtim,np.ndarray):
				raise FerrError("dtim is not an ndarray; all coord var data must be in np.ndarrays")

	if bool(app_in_t) == False:
		##if ndims < ct_dims:
		##raise FerrError("Not enough dim info to go by.  Please call with some dimensional data")
		if ndims != len(dims):
			raise FerrError("Dimensions and order not fully specified!")
	if ct_dims < len(dims):
		raise FerrError("Not enough dim info to go by.  Please call with more coord vars data")

	T = tname.lower()
	Z = zname.lower()
	Y = yname.lower()
	X = xname.lower()

	global outf
	outf = None

	if append is False:
		if os.path.isfile(nfn) == True:
			ques_str = "output file \'%s\' exists.  Delete? (y/n): " % nfn
			ans = raw_input(ques_str)
			if 'y' in ans.lower():
				os.remove(nfn)

		h1 = "created by ferr.py on %s, using py netCDF4" % datetime.datetime.ctime(datetime.datetime.now())
		outf = nc4.Dataset(nfn, 'w', clobber=False, format='NETCDF3_CLASSIC')
		outf.history = h1

		if 't' in dims:
			outf.createDimension(T,None)
			time_var = outf.createVariable(tname.lower(),'d',(T,))
			time_var.long_name = 'Time'
			time_var.axis = 'T'
			if dtim is not None:
				time_var.units = "days since 0001-01-01 00:00:00"
				time_var.time_origin = "0001-01-01 00:00:00"
				time_var[:] = nc4.date2num(dtim,units="days since 0001-01-01 00:00:00")
			elif t is not None:
				if tunits is None:
					raise FerrError("No time units specified, and only t-values given!")
				time_var.units = tunits
				time_var[:] = t
			if tbnds is not None:
				time_var.bounds = tname.lower() + "_bnds"
			if t_atts is not None:
				for ii in t_atts.keys():
					setattr(time_var,ii,t_atts[ii])

		if 'x' in dims:
			outf.createDimension(X,x.size)
			lon_var = outf.createVariable(xname.lower(),'d',(X,))
			lon_var.long_name = 'Longitude'
			lon_var.axis = 'X'
			lon_var.units = 'degrees_east'
			if xbnds is not None:
				lon_var.point_spacing = 'uneven'
				lon_var.bounds = xname.lower() + "_bnds"
			else:
				lon_var.point_spacing = 'even'
			lon_var.modulo = np.array([360.])
			lon_var[:] = x[:]

		if 'y' in dims:
			outf.createDimension(Y,y.size)
			lat_var = outf.createVariable(yname.lower(),'d',(Y,))
			lat_var.long_name = 'Latitude'
			lat_var.axis = 'Y'
			lat_var.units = 'degrees_north'
			if ybnds is not None:
				lat_var.point_spacing = 'uneven'
				lat_var.bounds = yname.lower() + "_bnds"
			else:
				lat_var.point_spacing = 'even'
			lat_var[:] = y[:]

		if 'z' in dims:
			outf.createDimension(Z,z.size)
			depth_var = outf.createVariable(zname.lower(),'d',(Z,))
			depth_var.long_name = 'Depth'
			depth_var.axis = 'Z'
			depth_var.units = 'meters'
			depth_var.positive = 'down'
			if zbnds is not None:
				depth_var.point_spacing = 'uneven'
				depth_var.bounds = zname.lower() + "_bnds"
			else:
				depth_var.point_spacing = 'even'
			depth_var[:] = z[:]

		if (tbnds is not None) | (xbnds is not None) | (ybnds is not None) | (zbnds is not None):
			outf.createDimension('bnds',2)

		if tbnds is not None:
			tbndsname = tname.lower() + "_bnds"
			tbnds_var = outf.createVariable(tbndsname,'d',(T,'bnds'))
			if isinstance(tbnds[0,0],dt):
				tbnds_var[:] = nc4.date2num(tbnds,units="days since 0001-01-01 00:00:00")
			elif isinstance(tbnds[0][0],float):
				tbnds_var[:] = tbnds

		if xbnds is not None:
			xbndsname = xname.lower() + "_bnds"
			xbnds_var = outf.createVariable(xbndsname,'d',(X,'bnds'))
			xbnds_var[:] = xbnds

		if ybnds is not None:
			ybndsname = yname.lower() + "_bnds"
			ybnds_var = outf.createVariable(ybndsname,'d',(Y,'bnds'))
			ybnds_var[:] = ybnds

		if zbnds is not None:
			zbndsname = zname.lower() + "_bnds"
			zbnds_var = outf.createVariable(zbndsname,'d',(Z,'bnds'))
			zbnds_var[:] = zbnds

	elif append is True:
		tax_found = False
		zax_found = False
		yax_found = False
		xax_found = False
		Tsize = 0

		if os.path.isfile(nfn) is False:
			raise FerrError("file %s does not exist; can't append to it" % nfn)

		outf = use(nfn,silent=True,_append=True)
		indims = outf.d.keys()
		invars = outf.v.keys()
		cv1 = outf.cv
		cv1_keys = cv1.keys()
		fdims = outf.f.dimensions.keys()

		if (dataname in invars) & (app_in_t is False):
			raise FerrError("variable %s already exists; can't append it to file %s" % (dataname,nfn))

		if (tbnds is not None) | (xbnds is not None) | (ybnds is not None) | (zbnds is not None):
			if 'bnds' not in fdims:
				outf.f.createDimension('bnds',2)

		if ('t' in dims) and (app_in_t is False):
			if (t is not None):
				Tsize = t.size
				for i in cv1_keys:
					if cv1[i] is 'tax':
						tax1 = outf.d[i][:]
						cmp1 = t == tax1
						##another option: sp.special.array_equiv()
						if isinstance(cmp1,np.ndarray):
							cmp1 = cmp1.all()
						if cmp1 is True:
							tax_found = True
							for fd in fdims:
								if i.lower() == fd.lower():
									T = fd
							break
			if (dtim is not None):
				Tsize = dtim.size
				for i in cv1_keys:
					if cv1[i] is 'tax':
						tax1 = outf.dt_vals(i)
						cmp1 = dtim == tax1
						if isinstance(cmp1,np.ndarray):
							cmp1 = cmp1.all()
						if cmp1:
							tax_found = True
							for fd in fdims:
								if i.lower() == fd.lower():
									T = fd
							break

			if tax_found is False:
				while T in indims:
					if T[-1].isdigit():
						T = T[0:-1] + str(int(T[-1]) + 1)
					else:
						T = T + '1'

				outf.f.createDimension(T,Tsize)
				time_var = outf.f.createVariable(T,'d',(T,))
				time_var.long_name = 'Time'
				time_var.axis = 'T'
				if dtim is not None:
					time_var.units = "days since 0001-01-01 00:00:00"
					time_var.time_origin = "0001-01-01 00:00:00"
					time_var[:] = nc4.date2num(dtim,units="days since 0001-01-01 00:00:00")
				elif t is not None:
					if tunits is None:
						raise FerrError("No time units specified, and only t-values given!")
					time_var.units = tunits
					time_var[:] = t
				if tbnds is not None:
					time_var.bounds = T + "_bnds"
				if t_atts is not None:
					for ii in t_atts.keys():
						setattr(time_var,ii,t_atts[ii])

				if tbnds is not None:
					tbndsname = T + "_bnds"
					tbnds_var = outf.f.createVariable(tbndsname,'d',(T,'bnds'))
					if isinstance(tbnds[0,0],dt):
						tbnds_var[:] = nc4.date2num(tbnds,units="days since 0001-01-01 00:00:00")
					elif isinstance(tbnds[0][0],float):
						tbnds_var[:] = tbnds

		if('z' in dims):
			if (z is not None):
				for i in cv1_keys:
					if cv1[i] is 'zax':
						zax1 = outf.d[i][:]
						cmp1 = z == zax1
						if isinstance(cmp1,np.ndarray):
							cmp1 = cmp1.all()
						if cmp1:
							zax_found = True
							for fd in fdims:
								if i.lower() == fd.lower():
									Z = fd
							break
				if zax_found is False:
					while Z in indims:
						if Z[-1].isdigit():
							Z = Z[0:-1] + str(int(Z[-1]) + 1)
						else:
							Z = Z + '1'

					outf.f.createDimension(Z,z.size)
					depth_var = outf.f.createVariable(Z,'d',(Z,))
					depth_var.long_name = 'Depth'
					depth_var.axis = 'Z'
					depth_var.units = 'meters'
					depth_var.positive = 'down'
					if zbnds is not None:
						depth_var.point_spacing = 'uneven'
						depth_var.bounds = Z + "_bnds"
					else:
						depth_var.point_spacing = 'even'
					depth_var[:] = z[:]

					if zbnds is not None:
						zbndsname = Z + "_bnds"
						zbnds_var = outf.f.createVariable(zbndsname,'d',(Z,'bnds'))
						zbnds_var[:] = zbnds

		if 'y' in dims:
			if (y is not None):
				for i in cv1_keys:
					if cv1[i] is 'yax':
						yax1 = outf.d[i][:]
						cmp1 = y == yax1
						if isinstance(cmp1,np.ndarray):
							cmp1 = cmp1.all()
						if cmp1:
							yax_found = True
							for fd in fdims:
								if i.lower() == fd.lower():
									Y = fd
							break
				if yax_found is False:
					while Y in indims:
						if Y[-1].isdigit():
							Y = Y[0:-1] + str(int(Y[-1]) + 1)
						else:
							Y = Y + '1'

					outf.f.createDimension(Y,y.size)
					lat_var = outf.f.createVariable(Y,'d',(Y,))
					lat_var.long_name = 'Latitude'
					lat_var.axis = 'Y'
					lat_var.units = 'degrees_north'
					if ybnds is not None:
						lat_var.point_spacing = 'uneven'
						lat_var.bounds = Y + "_bnds"
					else:
						lat_var.point_spacing = 'even'
					lat_var[:] = y[:]

					if ybnds is not None:
						ybndsname = Y + "_bnds"
						ybnds_var = outf.f.createVariable(ybndsname,'d',(Y,'bnds'))
						ybnds_var[:] = ybnds

		if 'x' in dims:
			if (x is not None):
				for i in cv1_keys:
					if cv1[i] is 'xax':
						xax1 = outf.d[i][:]
						cmp1 = x == xax1
						if isinstance(cmp1,np.ndarray):
							cmp1 = cmp1.all()
						if cmp1:
							xax_found = True
							for fd in fdims:
								if i.lower() == fd.lower():
									X = fd
							break
				if xax_found is False:
					while X in indims:
						if X[-1].isdigit():
							X = X[0:-1] + str(int(X[-1]) + 1)
						else:
							X = X + '1'

					outf.f.createDimension(X,x.size)
					lon_var = outf.f.createVariable(X,'d',(X,))
					lon_var.long_name = 'Longitude'
					lon_var.axis = 'X'
					lon_var.units = 'degrees_east'
					if xbnds is not None:
						lon_var.point_spacing = 'uneven'
						lon_var.bounds = X + "_bnds"
					else:
						lon_var.point_spacing = 'even'
					lon_var.modulo = np.array([360.])
					lon_var[:] = x[:]

					if xbnds is not None:
						xbndsname = X + "_bnds"
						xbnds_var = outf.f.createVariable(xbndsname,'d',(X,'bnds'))
						xbnds_var[:] = xbnds
	else:
		raise FerrError("kw 'append' has to be True or False")


	vdims = []
	for i in xrange(len(dims)):
		if 'x' is dims.lower()[i]:
			vdims.append(X)

		if 'y' is dims.lower()[i]:
			vdims.append(Y)

		if 'z' is dims.lower()[i]:
			vdims.append(Z)

		if 't' is dims.lower()[i]:
			vdims.append(T)

	vdims = tuple(vdims)
	if append is False:
		dat_var = outf.createVariable(dataname,'f',vdims,fill_value=-1.e+34)
		dat_var.missing_value = -1.e+34
		dat_var.long_name = dataname
		if dhist is not None:
			dat_var.history = dhist
		data = np.float32(data)
		data.set_fill_value(-1.e+34)
		data.data[data.mask] = data.fill_value
		if app_in_t is True:
			tsize = time_var.size
			if (tsize > 1) & (ndims == len(dims)):
				dims_a = nc4.stringtoarr(dims,len(dims))
				data_tsize = data.shape[np.where(dims_a == 't')[0]]
				dat_var[tsize:tsize+data_tsize,...] = data
			else:
				dat_var[0,...] = data
		else:
			dat_var[:] = data

	else:
		if app_in_t is True:
			if dataname.lower() not in outf.v.keys():
				raise FerrError("append mode, but var '%s' not found in file '%s'!" % (dataname,nfn))
			dat_var = outf.v[dataname.lower()]
			tim_var = outf.d[T]
			tsize = tim_var.size
			data = np.float32(data)
			data.set_fill_value(-1.e+34)
			data.data[data.mask] = data.fill_value
			if t is not None:
				t_insize = t.size
				tim_var[tsize:tsize + t_insize] = t
			elif dtim is not None:
				t_insize = dtim.size
				tim_var[tsize:tsize + t_insize] = nc4.date2num(dtim,units="days since 0001-01-01 00:00:00")
			else:
				raise FerrError("no time data, though append in time option was True")
			if (tsize > 1) & (ndims == len(dims)):
				dims_a = nc4.stringtoarr(dims,len(dims))
				data_tsize = data.shape[np.where(dims_a == 't')[0]]
				dat_var[tsize:tsize+data_tsize,...] = data
			else:
				dat_var[tsize,...] = data
		else:
			dat_var = outf.f.createVariable(dataname,'f',vdims,fill_value=-1.e+34)
			dat_var.missing_value = -1.e+34
			dat_var.long_name = dataname
			if dhist is not None:
				dat_var.history = dhist
			data = np.float32(data)
			data.set_fill_value(-1.e+34)
			data.data[data.mask] = data.fill_value
			dat_var[:] = data

	if append is False:
		outf.close()
	else:
		outf.f.close()

	if bool(silent) is False:
		print "\ndata written to %s\n" % nfn

	return None

#### end of save

##################################
## mpi_save
##################################


def mpi_save(nfn,data,datname,dtim,dtbnds=None):
	"""Save NETCDF file of 'data' which uses the MPI MILLENNIUM lat/lon
	coordinates.  Must pass the time axis as well, in an array of
	datetime objects.
	Input: 'nfn' - new filename
		'data' - the 3-d or 4-d data (either t,y,x or t,z,y,x)
		'datname' - data variable name
		'dtim' - array of date time objects for the T-axis
	"""
	if (data.ndim > 4) | (data.ndim < 3):
		raise FerrError("data array is not 3-d or 4-d")

	if (data.shape[0] != dtim.size):
		raise FerrError("input time array is of different size than first dim of input data")
	
	dfile = nc4.Dataset(nfn,'w',clobber=True,format='NETCDF4')
	
	clon,clat,clon_crnr,clat_crnr = mcp.mpi_grid()
	y1,x1 = clon.shape
	
	dfile.createDimension('time',dtim.size)
	dfile.createDimension('x',x1)
	dfile.createDimension('y',y1)
	dfile.createDimension('nv4',4)

	if dtbnds is not None:
		dfile.createDimension('bnds',2)

	time_var = dfile.createVariable('time','d',('time',))
	time_var.long_name = 'Time'
	time_var.units = 'days since 0001-01-01 00:00:00'
	time_var.axis = 'T'
	time_var.calendar = 'proleptic_gregorian'
	time_var[:] = mcmath.d2n(dtim)

	if dtbnds is not None:
		time_var.bounds = 'time_bnds'
		tbnds_var = dfile.createVariable('time_bnds','d',('time','bnds'))
		tbnds_var.units = 'days since 0001-01-01 00:00:00'
		tbnds_var.calendar = 'proleptic_gregorian'
		tbnds_var[:] = mcmath.d2n(dtbnds)
		
	lon_var = dfile.createVariable('lon','d',('y','x'))
	lon_var.long_name = 'longitude'
	lon_var.units = 'degrees'
	lon_var._CoordinateAxisType = 'Lon'
	lon_var.bounds = 'lon_bnds'
	lon_var[:] = clon

	lat_var = dfile.createVariable('lat','d',('y','x'))
	lat_var.long_name = 'latitude'
	lat_var.units = 'degrees'
	lat_var._CoordinateAxisType = 'Lat'
	lat_var.bounds = 'lat_bnds'
	lat_var[:] = clat

	lonbnds_var = dfile.createVariable('lon_bnds','d',('y','x','nv4'))
	lonbnds_var[:] = clon_crnr

	latbnds_var = dfile.createVariable('lat_bnds','d',('y','x','nv4'))
	latbnds_var[:] = clat_crnr

	if data.ndim == 4:
		dfile.createDimension('z',40)
		depth_var = dfile.createVariable('z','d',('z',))
		depth_var.axis = 'Z'
		depth_var.long_name = 'Depth'
		depth_var.positive = 'down'
		depth_var.point_spacing = 'uneven'
		depth_var.units = 'm'
		depth_var[:] = np.array([    6.,    17.,    27.,    37.,    47.,    57.,    69.,    83.,
		100.,   123.,   150.,   183.,   220.,   263.,   310.,   363.,
		420.,   485.,   560.,   645.,   740.,   845.,   960.,  1085.,
		1220.,  1365.,  1525.,  1700.,  1885.,  2080.,  2290.,  2525.,
		2785.,  3070.,  3395.,  3770.,  4195.,  4670.,  5170.,  5720.])

		if isinstance(data,np.ma.MaskedArray):
			data_var = dfile.createVariable(datname,'f',('time','z','y','x'),fill_value=data.fill_value)
			data_var.missing_value = data.fill_value
			data[data.mask] = data.fill_value
		else:
			data_var = dfile.createVariable(datname,'f',('time','z','y','x'))
		data_var.long_name = datname
		
		data_var[:] = data
	else:
		if isinstance(data,np.ma.MaskedArray):
			data_var = dfile.createVariable(datname,'f',('time','y','x'),fill_value=data.fill_value)
			data_var.missing_value = data.fill_value
			data[data.mask] = data.fill_value
		else:
			data_var = dfile.createVariable(datname,'f',('time','y','x'))
		data_var.long_name = datname
		
		data_var[:] = data

	dfile.sync()
	dfile.close()

	return None


#### end of mpi_save

