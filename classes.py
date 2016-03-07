import numpy as np
import scipy as sp
from scipy import stats as spstats

import cdo; cdo = cdo.Cdo()
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import ceil
import datetime
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset, num2date, date2num
import os
import datetime
from pcolor import default_pcolor_args
from taylor import TaylorDiagram
        
class Plot(object):
    def __init__(self, **kwargs):
        self.ifile = kwargs.pop('ifile', '')
        self.runID = kwargs.pop('runID', None)
        self.experiment = kwargs.pop('experiment', 'historical')
        self.variable = kwargs.pop('variable', None)
        self.data_type = kwargs.pop('data_type', 'climatology')
        self.projection = kwargs.pop('projection', 'global_map')
        self.start_date = str(kwargs.pop('start_date', self._default_start_date()))
        self.end_date = str(kwargs.pop('end_date', self._default_end_date()))
        self.remap = kwargs.pop('remap', 'remapdis')
        self.remap_grid = kwargs.pop('remap_grid', 'r360x180')
        self.scale = kwargs.pop('scale', 1)
        self.shift = kwargs.pop('shift', 0)
        self.frequency = kwargs.pop('frequency', 'mon')
        self.idealdepth = kwargs.pop('depth', None)
        self.realm = kwargs.pop('realm', 'atmos')
        self.realization = kwargs.pop('realization', 'r1i1p1')
        self.divergent = kwargs.pop('divergent', False)
        self.pcolor_args = dict(kwargs.pop('pcolor_args', {}))
        self.ax_args = dict(kwargs.pop('ax_args', {}))
        self.fig = kwargs.pop('fig', None)
        self.ax = kwargs.pop('ax', None)
        self.seasons = kwargs.pop('seasons', ['DJF', 'MAM', 'JJA', 'SON'])
        self.rmse = kwargs.pop('rmse', False)
        self.depthneeded = kwargs.pop('depthneeded', None)
        self.cdostring = kwargs.pop('cdostring', None)
        self.png = kwargs.pop('png', False)
        self.pdf = kwargs.pop('pdf', True)
        self.cbaxis = kwargs.pop('cbaxis', None)
        self.alpha = kwargs.pop('alpha', None)
        self._flatdata = kwargs.pop('flatdata', None)
        self._dataset = kwargs.pop('dataset', None)
        self._data = kwargs.pop('data', None)
        self._full_data = kwargs.pop('fulldata', None)
        self._full_dataset = kwargs.pop('full_dataset', None)
        self._full_ncvar = kwargs.pop('full_ncvar', None)
        self._ffile = kwargs.pop('ffile', None)
        self._units = kwargs.pop('units', None)
        self._lat = kwargs.pop('lat', None)
        self._lon = kwargs.pop('lon', None)
        self._slat = kwargs.pop('slat', None)
        self._slon = kwargs.pop('slon', None)
        self._ncvar = kwargs.pop('ncvar', None)
        self._realm_cat = kwargs.pop('realm_cat', None)
        self._plotname = kwargs.pop('plotname', None)
        self._plotdata = kwargs.pop('plotdata', None)
        self._unscaleddata = kwargs.pop('unscaleddata', None)
        self._truedepth = kwargs.pop('truedepth', None)
        self._ofile = kwargs.pop('ofile', None)
        self.sister = kwargs.pop('sister', None)
        self._statstring = kwargs.pop('statstring', None)
        self._pvalues = kwargs.pop('pvalues', None)
        self._plotstd = kwargs.pop('plotstd', None)
        self._std = kwargs.pop('std', None)
        
        self._time_averaged = None
        self._fulldepths = None
        self._stats = None
        self._standard_deviation = None
        self._trendscale = None
        self._nobs = None

        self.plot_args = self._fill_plot_args(**kwargs.pop('plot_args', {}))     
    
    def _get_datetime(self, datestring):
        datelist = ['%Y-%M-%d', '%Y-%M', '%Y']
        for form in datelist:
            try:
                obj = datetime.datetime.strptime(datestring, form)
            except:
                continue
            else:
                return obj       
    
    @property
    def nobs(self):
        if self._nobs is None:
            sd = self._get_datetime(self.start_date)
            ed = self._get_datetime(self.end_date)
            number_of_days = (ed - sd).days
            if self.frequency == 'day':
                self._nobs =  number_of_days
            elif self.frequency == 'mon':
                self._nobs = number_of_days / 30
            else:
                self._nobs = number_of_days / 365
        return self._nobs
    
    @property
    def plotstd(self):
        if self._plotstd is None:
            self._plotstd = self._removedepth(self.std) * self.scale
        return self._plotstd
    
    @property
    def std(self):
        if self._std is None:
            dataset = Dataset(self._calculate_std(self.ifile), 'r')  
            ncvar = dataset.variables[self.variable]
            self._std = ncvar[:].squeeze()
        return self._std    
    
    @property
    def slon(self):
        if self._slon is None:
            self._slon = [self.lon[lo] for (la, lo), value in np.ndenumerate(self.pvalues) if abs(value) < self.alpha and lo % 2 == 0]
        return self._slon
    
    @property
    def slat(self):
        if self._slat is None:
            self._slat = [self.lat[la] for (la, lo), value in np.ndenumerate(self.pvalues) if abs(value) < self.alpha and lo % 2 == 0]
        return self._slat
                
    @property
    def pvalues(self):
        if self._pvalues is None:
            t, self._pvalues = sp.stats.ttest_ind_from_stats(self.plotdata, self.plotstd, self.nobs, self.sister.plotdata, self.sister.plotstd, self.sister.nobs, equal_var=False)
#            t, self._pvalues = sp.stats.ttest_ind(self.full_data, self.sister.full_data, axis=0, equal_var=False)
            print self._pvalues.shape
            print self.pvalues[0]
        return self._pvalues      
        
    @property
    def trendscale(self):
        if self._trendscale is None:
            if self.data_type == 'trends':
                if self.frequency == 'day':
                    self._trendscale = 3650
                if self.frequency == 'mon':
                    self._trendscale = 120
                if self.frequency == 'yr':
                    self._trendscale = 10
            else:
                 self._trendscale = 1
        return self._trendscale
                  
    @property            
    def realm_cat(self):
        if self._realm_cat is not None:
            if self.realm == 'aerosol' or self.realm == 'atmos' or self.realm == 'seaIce':
                self._realm_cat = 'atmos'
            elif self.realm == 'land' or self.realm == 'landIce':
                self._realm_cat = 'land'
            else:
                self._realm_cat = 'ocean'
        return self._realm_cat
            
    def _default_start_date(self):
        if 'piControl' in self.experiment:
            return '2900-01'
        elif 'rcp' in self.experiment:
            return '2070-01'
        return '1900-01'

    def _default_end_date(self):
        if 'piControl' in self.experiment:
            return '3000-01'
        elif 'rcp' in self.experiment:
            return '2100-01'
        return '2005-01'
                        
    def _fill_plot_args(self, **kwargs):
        d_args =  {'fill_continents': False,
                   'draw_parallels': False,
                   'draw_meridians': False,
                   'cblabel': self.units,
                   }
        for key in d_args:
            if key not in kwargs:
                kwargs[key] = d_args[key]
        return kwargs
        
    def _fill_pcolor_args(self):
        if self.data_type == 'trends' or self.divergent:
            d1pca = default_pcolor_args(self.plotdata, anom=True)
        else:
            d1pca = default_pcolor_args(self.plotdata)
       
        if self.sister is not None:
            if self.data_type == 'trends' or self.divergent:
                d2pca = default_pcolor_args(self.sister.plotdata, anom=True)
            else:
                d2pca = default_pcolor_args(self.sister.plotdata) 

            vmin = np.min([d1pca['vmin'], d2pca['vmin']])
            vmax = np.max([d1pca['vmax'], d2pca['vmax']])
            d1pca['vmin'] = vmin
            d1pca['vmax'] = vmax            
        for key in self.pcolor_args:
            d1pca[key] = self.pcolor_args[key]
        self.pcolor_args = d1pca         

    def _default_ax_args(self):
        return {'title': self.variable + ' Model: ' + self.runID}
    
    def _fill_ax_args(self):        
        args = self._default_ax_args()           
        for key in self.ax_args:
            args[key] = self.ax_args[key]       
        self.ax_args = args            
    
    @property
    def plotname(self):
        if self._plotname == None:      
            self._plotname = ('plots/' + self.variable + '_' + self.runID + '_' + self.projection + '_' + 
                              self.data_type + '_' + self.start_date + '-' + self.end_date + 
                              '_' + ''.join(self.seasons))
            try:
                self._plotname.append('_' + str(self.truedepth))
            except:
                pass
        return self._plotname

    def _savefigures(self):
        pdfname = self.plotname + '.pdf'
        pngname = self.plotname + '.png'
        if self.png:
            plt.savefig(pngname, bbox_inches='tight')
        if self.pdf:
            plt.savefig(pdfname, bbox_inches='tight')
    
    @property
    def statstring(self):
         if self._statstring is None:
             vals = self.stats.values()
             ss = [s + ': ' for s in self.stats.keys()]
             val = [s + str(v) for s, v in zip(ss, vals)]
             self._statstring = ' '.join(val)
             
         return self._statstring
    
    @property        
    def stats(self):
        if self._stats is None:
            if self.rmse:
                vals = [str(np.round(self.plotdata.min(), 1)), str(np.round(self.plotdata.max(), 1)), str(np.round(np.sqrt(np.mean(np.square(self.plotdata))), 1))]
                snam = ['min: ', 'max: ', 'rmse: ']
                self._stats = {'rmse': float(vals[2]),
                              'min': float(vals[0]),
                              'max': float(vals[1]),
                              }
            else:
                vals = [str(np.round(self.plotdata.min(), 1)), str(np.round(self.plotdata.max(), 1)), str(np.round(self.plotdata.mean(), 1))]
                snam = ['min: ', 'max: ', 'mean: ']
                self._stats = {'mean': float(vals[2]),
                              'min': float(vals[0]),
                              'max': float(vals[1]),
                              }
        return self._stats
        
    @property
    def fulldepths(self):
        if self._fulldepths is None:
            self._fulldepths = 'surface'
            for dimension in self.ncvar.dimensions:
                try:
                    if self.dataset.variables[dimension].axis == 'Z':
                        self._fulldepths = self.dataset.variables[dimension][:]
                        break
                except:
                     pass # keep looping if the dimension doesn't have an 'axis' attribute
        return self._fulldepths
    
    def plot(self):
        self._fill_pcolor_args()
        self._fill_ax_args()
        self._makeplot()
        self._savefigures()
        
    def _makeplot(self):
        self._draw(**self.plot_args)
    
    @property
    def time_averaged(self):
        if self._time_averaged is not None:
            nc = Dataset(ifile, 'r')
            try: 
                time = nc.variables['time'][:].squeeze()
            except:
                self._time_averaged = True
            else:
                self._time_averaged = time.size == 1
        return self._time_averaged

    @property
    def flatdata(self):
        if self._flatdata is None:
            self._flatdata = self.plotdata.flatten()
        return self._flatdata
    
    @property    
    def standard_deviation(self):
        if self._standard_deviation is None:
            self._standard_deviation = float(self.plotdata.std())
        return self._standard_deviation
    
    @property
    def full_plotdata(self):
        pass
        
    @property
    def full_data(self):
        if self._full_data is None:
            self._full_data = self.full_ncvar[:].squeeze()
        return self._full_data
    
    @property
    def full_ncvar(self):
        if self._full_ncvar is None:
            self._full_ncvar = self.full_dataset.variables[self.variable]
        return self._full_ncvar
    
    @property
    def full_dataset(self):
        if self._full_dataset is None:
            self._full_dataset = Dataset(self.ffile, 'r')
        return self._full_dataset
    
    @property
    def ffile(self):
        if self._ffile is None:
            self._ffile = self._calculate_full(self.ifile)
        return self._ffile   
    
    @property
    def data(self):
        if self._data is None:
            self._data = self.ncvar[:].squeeze()
        return self._data
            
    @property
    def ncvar(self):
        if self._ncvar is None:
            self._ncvar = self.dataset.variables[self.variable]
        return self._ncvar
        
    @property
    def units(self):
        if self._units is None:
            self._units = self.ncvar.units
            if self.data_type == 'trends':
                self._units += ' / decade'                
        return self._units
                
    @property
    def lat(self):
        if self._lat is None:
            self._lat = self.dataset.variables['lat'][:].squeeze() 
        return self._lat
    
    @property
    def lon(self):
        if self._lon is None:
            self._lon = self.dataset.variables['lon'][:].squeeze()
        return self._lon
 
    @property
    def truedepth(self):
        if self._truedepth is None:
            if self.fulldepths == 'surface':
                self._truedepth = 'surface'
            else:
                try:
                    self._truedepth = min(self.fulldepths, key=lambda x: abs(x - self.depthneeded))
                except:
                    try:
                        self._truedepth = min(self.fulldepths, key=lambda x: abs(x - self.idealdepth))
                    except:
                        self._truedepth = self.fulldepths
        return self._truedepth

    @property
    def plotdata(self):
        if self._plotdata is None:
            self._plotdata = (self.unscaleddata + self.shift) * self.scale * self.trendscale          
        return self._plotdata
    
    def _removedepth(self, idata):
        if len(idata.shape) < 3:
            return idata
        else:
            depth_index = np.where(np.round(self.fulldepths) == np.round(self.truedepth))[0][0]
            return idata[depth_index, :, :]  
                 
    @property
    def unscaleddata(self):
        if self._unscaleddata is None:
            self._unscaleddata = self._removedepth(self.data)
        return self._unscaleddata
    
    def _calculate(self, ifile):
        if self.cdostring is not None:
            return self._calculate_cdo(ifile)
        elif self.data_type == 'trends':
            return self._calculate_trends(ifile)
        else:
            return self._calculate_climatology(ifile)
        
    @property
    def ofile(self):
        if self._ofile is None:
            self._ofile = self._calculate(self.ifile)
        return self._ofile
    
    @property
    def dataset(self):
        if self._dataset is None:
            self._dataset = Dataset(self.ofile, 'r')
        return self._dataset
    
    def _calculate_cdo(self, filename):
        return self._cdos(filename)
    
    def _calculate_std(self, filename):
        return self._remap(self._setc(self._time_std(self._mask(self._season(self._sel_var(filename))))))
        
    def _calculate_climatology(self, filename):
        return self._remap(self._setc(self._time_mean(self._mask(self._season(self._sel_var(filename))))))

    def _calculate_trends(self, filename):
        return self._remap(self._trend(self._setc(self._season(self._sel_var(filename)))))

    def _calculate_full(self, filename):
        return self._remap(self._setc(self._sel_date(self._mask(self._season(self._sel_var(filename))))))
    
    def _split(self, name):
        path, filename = os.path.split(name)
        return filename
        
    def _sel_date(self, name):
        if self.time_averaged:
            return name       
        out = 'netcdf/seldate_' + self.start_date + '_' + self.end_date + '_' + self._split(name)
        if not os.path.isfile(out):
            datestring = self.start_date + ',' + self.end_date
            cdo.seldate(datestring, input=name, output=out)
        return out

    def _sel_var(self, name):
        out = 'netcdf/sel_' + self._split(name)
        if not os.path.isfile(out):
            cdo.selvar(self.variable, input=name, output=out)
        return out
        
    def _mask(self, name):
        out = 'netcdf/masked_' + self._split(name)
        if not os.path.isfile(out):
            if self.realm_cat == 'ocean':
                try:
                    cdo.ifthen(input='mask/ocean ' + name, output=out)
                except:
                    with open('logs/log.txt', 'a') as outfile:
                        outfile.write('WARNING: Land data was not masked\n')
                    return name             
            elif self.realm_cat == 'land':
                try:
                    cdo.ifthen(input='mask/land ' + name, output=out) 
                except:
                    with open('logs/log.txt', 'a') as outfile:
                        outfile.write('WARNING: Ocean data was not masked\n')
                    return name 
            else:
                out = name
        return out
    
    def _time_std(self, name):
        if self.time_averaged:
            return name
        out = 'netcdf/std_' + self.start_date + '_' + self.end_date + '_' + self._split(name)
        if not os.path.isfile(out):
            seldatestring = '-seldate,' + self.start_date + ',' + self.end_date
            cdo.timstd(input=seldatestring + ' ' + name, output=out)
        return out
    
    def _time_mean(self, name):
        if self.time_averaged:
            return name    
        out = 'netcdf/climate_' + self.start_date + '_' + self.end_date + '_' + self._split(name)
        if not os.path.isfile(out):
            seldatestring = '-seldate,' + self.start_date + ',' + self.end_date
            cdo.timmean(input=seldatestring + ' ' + name, output=out)
        return out
    
    def _trend(self, name):
        out = 'netcdf/slope_' + self.start_date + '_' + self.end_date + '_' + self._split(name)
        outintercept = 'netcdf/intercept_' + self._split(name)
        if not os.path.isfile(out):
            seldatestring = '-seldate,' + self.start_date + ','  + self.end_date
            cdo.trend(input=seldatestring + ' ' + name, output=outintercept + ' ' + out)
        return out        

    def _detrend(self, name):
        out = 'netcdf/detrend_' + self._split(name)
        if not os.path.isfile(out):
            seldatestring = '-seldate,' + self.start_date + ',' + self.end_date
            cdo.detrend(input=seldatestring + ' ' + name, output=out)
        return out
        
    def _setc(self, name):
        if self.realm_cat == 'atmos':
            return name
        out = 'netcdf/setc_' + self._split(name)
        if not os.path.isfile(out):
            cdo.setctomiss(0, input=name, output=out)
        return out

    def _remap_function(self, remap):
        """ Returns a cdo function from string of the same name.
        """
        def cdoremap(r):
            return {'remapbil': cdo.remapbil,
                    'remapbic': cdo.remapbic,
                    'remapdis': cdo.remapdis,
                    'remapnn': cdo.remapnn,
                    'remapcon': cdo.remapcon,
                    'remapcon2': cdo.remapcon2,
                    'remapplaf': cdo.remaplaf,
                    }[r]  
        return cdoremap(remap) 
                    
    def _remap(self, name):
        out = 'netcdf/' + self.remap + '-' + self.remap_grid + self._split(name)
        if not os.path.isfile(out):
            remap = self._remap_function(self.remap)
            remap(self.remap_grid, input=name, output=out)
        return out

    def _field_mean(self, name):
        out = 'netcdf/fldmean_' + self.start_date + '_' + self.end_date + '_' + self._split(name)
        if not os.path.isfile(out):
            seldatestring = '-seldate,' + self.start_date + ',' + self.end_date
            cdo.fldmean(options='-L', input=seldatestring + ' ' + name, output=out)
        return out
        
    def _zonal_mean(self, name):
        out = 'netcdf/zonmean_' + self._split(name)
        if not os.path.isfile(out):
            cdo.zonmean(input=name, output=out)
        return out

    def _season(self, name):
        if self.seasons == None or self.seasons == ['DJF', 'MAM', 'JJA', 'SON']:
            return name
        seasonstring = ','.join(self.seasons)
        outputstring = ''.join(self.seasons)
        out = 'netcdf/selseason-' + outputstring + '_' + self._split(name)
        if not os.path.isfile(out):
            cdo.selseas(seasonstring, input=name, output=out)
        return out

    def _depthstring(self, depthlist):
        try:
            depthlist = list(depthlist)
        except:
            depthlist = [float(depthlist)]
        depthneed = ["%.2f" % number for number in depthlist]
        for i in xrange(len(depthneed)):
            depthneed[i] = str(depthneed[i])
        return ','.join(depthneed)  
          
    def _intlevel(self, name):
        if self.depthneeded == None or self.depthneeded == 'surface':
            return name
        depth = self._depthstring(self.depthneeded)
        depthname = depth.replace(' ', '')
        if len(depthname) > 100:
            depthname = depthname[:99]
        out = 'netcdf/level-' + str(depthname) + '_' + self._split(name)
        if depth:
            if not os.path.isfile(out):
                try:
                    cdo.intlevelx(str(depth), input=name, output=out)
                except:
                    return name
        else:
            return name
        return out  

    def _cdos(self, name):
        out = 'netcdf/cdo_' + split(name)
        if not os.path.isfile(out):
            s = 'cdo ' + string + ' ' + name + ' ' + out
            os.system(s)
        return out
           
    def _draw(self, latmin=-80, latmax=80, lonmin=0, lonmax=360, lon_0=0,
             fill_continents=False, draw_parallels=False, draw_meridians=False, cblabel='',
             **kwargs):
        if not self.ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))
        else:
            self.fig = self.ax.get_figure()
    
        if self.projection == 'global_map':
            m = Basemap(projection='kav7', llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, lon_0=-180, resolution='c', ax=self.ax)  
            labx, laby = m(10, -88) 
        if self.projection == 'polar_map':
            m = Basemap(projection='npstere', boundinglat=latmin, lon_0=lon_0, resolution='c', ax=self.ax)
            labx, laby = m(-135, 12)
        if self.projection == 'polar_map_south':
            m = Basemap(projection='spstere', boundinglat=latmax, lon_0=lon_0, resolution='c', ax=self.ax)
            labx, laby = m(-45, -12)
        if self.projection == 'mercator':  
            m = Basemap(projection='merc', llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, lat_ts=20, resolution='c', ax=self.ax) 
            labx, laby = m(lonmin + 1, latmin + 1)
    
        lons, lats = np.meshgrid(self.lon, self.lat)
        x, y = m(lons, lats)

        cot = m.pcolor(x, y, self.plotdata, **self.pcolor_args)
        plt.setp(self.ax, **self.ax_args)
        
        if self.alpha is not None:
            a, b = m(self.slon, self.slat)
            m.plot(a, b, '.', markersize=0.2, color='k', zorder=1)

        m.drawcoastlines(linewidth=1.25, ax=self.ax)
        if fill_continents:
            m.fillcontinents(color='0.8', ax=self.ax, zorder=2)
        if draw_parallels:
            m.drawparallels(np.arange(-80, 81, 20), labels=[1, 0, 0, 0], linewidth=0, ax=self.ax)
        if draw_meridians:
            m.drawmeridians(np.arange(0, 360, 90), labels=[0, 0, 0, 1], linewidth=0, yoffset=0.5e6, ax=self.ax)
        m.colorbar(mappable=cot, location='right', label=cblabel)
        self.ax.text(labx, laby, self.statstring, fontsize=8)
        
       
class Section(Plot):
    def __init__(self, **kwargs):
        super(Section, self).__init__(**kwargs) 
        self.set_yscale = kwargs.pop('set_yscale', 'log')
        
    @property
    def unscaleddata(self):
        if self._unscaleddata is None:
            if self.data.ndim == 3:
                self._unscaleddata = self.data.mean(axis=2)
            else:
                self._unscaleddata = self.data
        return self._unscaleddata
    
    def _draw(self, latmin=-80, latmax=80, lonmin=0, lonmax=360, lon_0=0, cblabel='',
              **kwargs):
        if not self.ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(8, 3))
        else:
            self.fig = self.ax.get_figure()
        cot = self.ax.pcolormesh(self.lat, self.fulldepths, self.plotdata, **self.pcolor_args)
        self.ax.contour(self.lat, self.fulldepths, self.plotdata, colors=['k'], 
                        vmin=self.pcolor_args['vmin'], vmax=self.pcolor_args['vmax'])
        self.ax.invert_yaxis()
        self.ax.autoscale(True, axis='both', tight='both')
        self.ax.set_yscale(self.set_yscale)
        if self.ax_args:
            plt.setp(self.ax, **self.ax_args)

        box = self.ax.get_position()
        if self.cbaxis:
            self.fig.colorbar(cot, cax=cbaxis, label=cblabel)
        else:
            tl = self.fig.add_axes([box.x1 + box.width * 0.05, box.y0, 0.02, box.height])
            self.fig.colorbar(cot, cax=tl, label=cblabel)


class Compare(Plot):
    def __init__(self, **kwargs):
        super(Compare, self).__init__(**kwargs)
        self.depthneeded = kwargs.pop('depthneeded', self.sister.truedepth)
    
    def _calculate_climatology(self, filename):
        return self._intlevel(super(self.__class__, self)._calculate_climatology(filename))
    
    def _calculate_trends(self, filename):
        return self._intlevel(super(self.__class__, self)._calculate_climatology(filename))  

    def _default_ax_args(self):
        return {'title': self.variable + ' Model: ' + self.runID}


class CompareSection(Section, Compare):
    def __init__(self, **kwargs):
        super(CompareSection, self).__init__(**kwargs)


class Line(Plot):
    def __init__(self, **kwargs):
        super(Line, self).__init__(**kwargs)
        self._xaxis = kwargs.pop('xaxis', None)
    
    def _gettime(self, dataset):
        nc_time = dataset.variables['time']
        try:
            cal = nc_time.calendar
        except:
            cal = 'standard'
        date = num2date(nc_time[:], nc_time.units, cal)
        x = [datetime.datetime(*item.timetuple()[:6]) for item in date]
        return np.array(x)
    
    def _removedepth(self, idata):
        if len(idata.shape) < 2:
            return idata
        else:
            depth_index = np.where(np.round(self.fulldepths) == np.round(self.truedepth))[0][0]
            if self.projection == 'time_series':
                return idata[:, depth_index]
            else:
                return idata[depth_index, :]             
        
    @property
    def xaxis(self):
        if self._xaxis is None:
            if self.projection == 'time_series':
                self._xaxis = self._gettime(self.dataset)
            else:
                self._xaxis = self.lat
        return self._xaxis
    
    def _calculate_climatology(self, filename):
        if self.projection == 'time_series':
            return self._field_mean(self._setc(self._mask(self._season(self._sel_var(filename)))))
        else:
            return self._zonal_mean(self._remap(self._time_mean(self._setc(self._mask(self._season(self._sel_var(filename)))))))    
    
    def _calculate_trends(self, filename):
        if self.projection == 'time_series':
            return self._field_mean(self._trend(self._setc(self._mask(self._season(self._sel_var(filename))))))
        else:
            return self._zonal_mean(self._remap(self._trend(self._setc(self._mask(self._season(self._sel_var(filename)))))))
    
    def _default_ax_args(self):
        args = {'title': self.variable + ' Model: ' + self.runID}
        args['ylabel'] = str(self.units)
        if self.projection == 'time_series':
            args['xlabel'] = 'Time'
        else:
            args['xlabel'] = 'Latitude'
        return args
        
    def _draw(self, **kwargs):
        if not self.ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))
        else:
            self.fig = self.ax.get_figure()
        
        self.ax.plot(self.xaxis, self.plotdata, label=self.runID, zorder=10)
        plt.setp(self.ax, **self.ax_args)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

class CompareLine(Line, Compare):
    def __init__(self, **kwargs):
        super(CompareLine, self).__init__(**kwargs)
        self.labelled = kwargs.pop('labelled', None)
    
    def _calculate_climatology(self, filename):
        return Compare._calculate_climatology(self, filename)
    
    def _calculate_trends(self, filename):
        return Compare._calculate_trends(self, filename)
    
    def _draw(self, **kwargs):
        if not self.ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))
        else:
            self.fig = self.ax.get_figure()

        
        if self.labelled:
             self.ax.plot(self.xaxis, self.plotdata, label=self.runID, zorder=3)
        else:
             self.ax.plot(self.xaxis, self.plotdata, color='0.5', label='remove', zorder=1)
        
        handles, labels = self.ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        try:
            by_label.pop('remove')
        except: pass
        plt.legend(by_label.values(), by_label.keys())            
        
class Figure(Plot):
    def __init__(self, **kwargs):
        super(Figure, self).__init__(**kwargs) 

    def _calculate_climatology(self, filename):
        return self._remap(self._setc(self._time_mean(self._mask(self._season(self._sel_var(filename))))))

    def _calculate_trends(self, filename):
        return self._setc(self._remap(self._trend(self._season(self._sel_var(filename)))))

class Taylor(Figure, Compare):
    def __init__(self, **kwargs):
        super(Taylor, self).__init__(**kwargs)
        self.plot_args = {}
        self.labelledifiles = kwargs.pop('labelledifiles', {})
        self.unlabelledifiles = kwargs.pop('unlabelledifiles', [])
        self._stdrange = kwargs.pop('stdrange', None)
        
        self._labelledflatdata = kwargs.pop('labelledflatdata', None)
        self._unlabelledflatdata = kwargs.pop('unlabelledflatdata', None)
        self._labelledstandard_deviations = kwargs.pop('labelledstandard_deviations', None)
        self._unlabelledstandard_deviations = kwargs.pop('unlabelledstandard_deviations', None)
        self._labelledcorrelation_coefficients = kwargs.pop('labelledcorrelation_coefficients', None)
        self._unlabelledcorrelation_coefficients = kwargs.pop('unlabelledcorrelation_coefficients', None)
        
        
        self._labelledplotdata = kwargs.pop('labelledplotdata', None)
        self._unlabelledplotdata = kwargs.pop('unlabelledplotdata', None)
        self._unlabelledunscaleddata = kwargs.pop('unlabelledunscaleddata', None)
        self._labelledunscaleddata = kwargs.pop('labelledunscaleddata', None)
        self._labelleddata = kwargs.pop('labelleddata', None)
        self._unlabelleddata = kwargs.pop('unlabelleddata', None)
        self._labelledncvar = kwargs.pop('labelledncvar', None)
        self._unlabelledncvar = kwargs.pop('unlabelledncvar', None)
        self._labelleddatasets = kwargs.pop('labelleddatasets', None)
        self._unlabelleddatasets = kwargs.pop('unlabelleddatasets', None)
        self._labelledofiles = kwargs.pop('labelledofiles', None)
        self._unlabelledofiles = kwargs.pop('unlabelledofiles', None)
        
        self._labelledsamples = kwargs.pop('labelledsamples', None)
        self._unlabelledsamples = kwargs.pop('unlabelledsamples', None)

    def plot(self):
        self._fill_ax_args()
        self._makeplot()
        self._savefigures()
            
    @property
    def stdrange(self):
        if self._stdrange is None:
           combined_stds = list(self.unlabelledstandard_deviations)
           combined_stds.extend(self.labelledstandard_deviations.values())
           self._stdrange = max(combined_stds) * 1.3 / self.sister.standard_deviation
           if self._stdrange <= 1.5:
                self._stdrange = 1.5
        return self._stdrange 
                     
    @property
    def labelledflatdata(self):
        if self._labelledflatdata is None:
            self._labelledflatdata = {label: data.flatten() for label, data in self.labelledplotdata.iteritems()}
        return self._labelledflatdata        
        
    @property
    def unlabelledflatdata(self):
        if self._unlabelledflatdata is None:
            self._unlabelledflatdata = [data.flatten() for data in self.unlabelledplotdata]
        return self._unlabelledflatdata
    
    @property
    def labelledstandard_deviations(self):
        if self._labelledstandard_deviations is None:
            self. _labelledstandard_deviations = {label: float(data.std()) for label, data in self.labelledplotdata.iteritems()}
        return self._labelledstandard_deviations
        
    @property
    def unlabelledstandard_deviations(self):
        if self._unlabelledstandard_deviations is None:
            self._unlabelledstandard_deviations = [float(data.std()) for data in self.unlabelledplotdata]
        return self._unlabelledstandard_deviations
    
    @property
    def labelledcorrelation_coefficients(self):
        if self._labelledcorrelation_coefficients is None:
            self._labelledcorrelation_coefficients = {label: np.ma.corrcoef(self.sister.flatdata, data)[0, 1] for label, data in self.labelledflatdata.iteritems()}
        return self._labelledcorrelation_coefficients
    
    @property
    def unlabelledcorrelation_coefficients(self):
        if self._unlabelledcorrelation_coefficients is None:
            self._unlabelledcorrelation_coefficients = [np.ma.corrcoef(self.sister.flatdata, data)[0, 1] for data in self.unlabelledflatdata]
        return self._unlabelledcorrelation_coefficients
    
    
    @property
    def labelledplotdata(self):
        if self._labelledplotdata is None:
            self._labelledplotdata = {label: (data + self.shift) * self.scale * self.trendscale for label, data in self.labelledunscaleddata.iteritems()}  
        return self._labelledplotdata
    
    @property
    def unlabelledplotdata(self):
        if self._unlabelledplotdata is None:
            self._unlabelledplotdata = [(data + self.shift) * self.scale * self.trendscale for data in self.unlabelledunscaleddata]
        return self._unlabelledplotdata
    
    @property
    def labelledunscaleddata(self):
        if self._labelledunscaleddata is None:
            self._labelledunscaleddata = {label: self._removedepth(data) for label, data in self.labelleddata.iteritems()}
        return self._labelledunscaleddata
    
    @property
    def unlabelledunscaleddata(self):
        if self._unlabelledunscaleddata is None:
            self._unlabelledunscaleddata = [self._removedepth(data) for data in self.unlabelleddata]   
        return self._unlabelledunscaleddata
    
    @property
    def labelleddata(self):
        if self._labelleddata is None:
            self._labelleddata = {label: ncvar[:].squeeze() for label, ncvar in self.labelledncvar.iteritems()}
        return self._labelleddata
    
    @property
    def unlabelleddata(self):
        if self._unlabelleddata is None:
            self._unlabelleddata = [ncvar[:].squeeze() for ncvar in self.unlabelledncvar] 
        return self._unlabelleddata

    @property
    def labelledncvar(self):
        if self._labelledncvar is None:
            self._labelledncvar = {label: dataset.variables[self.variable] for label, dataset in self.labelleddatasets.iteritems()}
        return self._labelledncvar
    
    @property
    def unlabelledncvar(self):
        if self._unlabelledncvar is None:
            self._unlabelledncvar = [dataset.variables[self.variable] for dataset in self.unlabelleddatasets]
        return self._unlabelledncvar
    
    @property
    def labelleddatasets(self):
        if self._labelleddatasets is None:
            self._labelleddatasets = {label: Dataset(ofile, 'r') for label, ofile in self.labelledofiles.iteritems()}
        return self._labelleddatasets
   
    @property
    def unlabelleddatasets(self):
        if self._unlabelleddatasets is None:
            self._unlabelleddatasets = [Dataset(ofile, 'r') for ofile in self.unlabelledofiles]
        return self._unlabelleddatasets
    
    @property
    def labelledofiles(self):
        if self._labelledofiles is None:
            self._labelledofiles = {}
            for label in self.labelledifiles:
                try:
                    self.labelledofiles[label] = self._calculate(self.labelledifiles[label])
                except:
                    print 'could not append: ' + label + ': ' + ifile                  
        return self._labelledofiles
        
    @property
    def unlabelledofiles(self):
        if self._unlabelledofiles is None:
            ofiles = []
            for ifile in self.unlabelledifiles:
                 try:
                     ofiles.append(self._calculate(ifile))
                 except:
                     print 'could not append: ' + ifile
            self._unlabelledofiles = ofiles                
        return self._unlabelledofiles
    
    @property
    def labelledsamples(self):
        if self._labelledsamples is None:
            self._labelledsamples = {label: (self.labelledstandard_deviations[label], self.labelledcorrelation_coefficients[label]) for label in self.labelledstandard_deviations}
        return self._labelledsamples
    
    @property    
    def unlabelledsamples(self):
        if self._unlabelledsamples is None:
            self._unlabelledsamples = zip(self.unlabelledstandard_deviations, self.unlabelledcorrelation_coefficients)            
        return self._unlabelledsamples
    
    def _draw(self, **kwargs):
        if not self.ax:
            self.fig =  plt.figure()
        else:
            self.fig = self.ax.get_figure()
        
        dia = TaylorDiagram(self.sister.standard_deviation, fig=self.fig, label='obs', srange=(0, self.stdrange))
        colors = plt.matplotlib.cm.jet(np.linspace(0,1, len(self.labelledsamples)))
         
        for i, (label, (stddev, corrcoef)) in enumerate(self.labelledsamples.iteritems()):
            dia.add_sample(stddev, corrcoef,
                           marker='.', ms=12, ls='',
                           mfc=colors[i], mec=colors[i],
                           label=label, zorder=2)
        self.fig.legend(dia.samplePoints,
                         [p.get_label() for p in dia.samplePoints],
                         numpoints=1, prop=dict(size='small'), loc='upper right')
        
        for (stddev, corrcoef) in self.unlabelledsamples:
            dia.add_sample(stddev, corrcoef,
                           marker='.', ms=12, ls='',
                           mfc='grey', mec='grey',
                           label=None, zorder=1)
        
        dia.add_grid()
        contours = dia.add_contours(colors='0.5')
        plt.clabel(contours, inline=1, fontsize=10)
        plt.title(self.variable)

class Histogram(Plot):
    def __init__(self, **kwargs):
        super(Histogram, self).__init__(**kwargs)

    def _default_ax_args(self):
        return {'title': self.variable + ' Model: ' + self.runID,
                'ylabel': '# Realizations',
                'xlabel': self.units}
                
    def _removedepth(self, idata):
        if self.fulldepths == 'surface':
            return idata
        else:
            depth_index = np.where(np.round(self.fulldepths) == np.round(self.truedepth))[0][0]
            return idata[depth_index]
       
    def _calculate_climatology(self, filename):
        return self._field_mean(self._remap(self._time_mean(self._setc(self._mask(self._season(self._sel_var(filename)))))))
         
    def _calculate_trends(self, filename):
        return self._field_mean(self._remap(self._trend(self._setc(self._season(self._sel_var(filename))))))
        
class HistogramCompare(Histogram, Taylor):
    def __init__(self, **kwargs):
        super(HistogramCompare, self).__init__(**kwargs)
   
    def _draw(self, **kwargs):
        if not self.ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(8,8))
        else:
            self.fig = self.ax.get_figure()
        
        n, bins, patches = plt.hist(self.unlabelledplotdata, 10, facecolor='grey', alpha=0.75)
        ymax = int(ceil(1.2 * max(n)))
        self.ax.set_ylim(0, ymax)
        
        for key, value in self.labelledplotdata.iteritems():
            plt.axvline(value, label=key, linewidth=4,
                        color=next(self.ax._get_lines.color_cycle))
        plt.setp(self.ax, **self.ax_args)
        plt.title(self.variable)
        
        self.ax.legend(loc='best')
            
            
if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    s = Plot(ifile='/raid/rc40/data/ncs/historical-cvu/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-1_historical-cvu_r1i1p1_185001-200512.nc',
             runID='cvu',
             experiment='historical',
             variable='ta',
             depth=100000,
             data_type='climatology',
             projection='mercator',
             realm='ocean',                         
             )
    q = Compare(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID='edr',
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='mercator',
             realm='atmos',
             sister=s
             )
    s.sister = q
    r = Plot(plotdata = s.plotdata - q.plotdata,
             runID='comp',
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='mercator',
             realm='atmos',
             alpha=0.05,
             pvalues = s.pvalues,
             divergent=True,
             dataset = s.dataset
             )
    r.plot()  
    """    
    s = Plot(ifile='/raid/rc40/data/ncs/historical-cvu/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-1_historical-cvu_r1i1p1_185001-200512.nc',
             runID='cvu',
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='time_series',

             realm='ocean',                         
             )
   
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    q = Line(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID='edr',
             experiment='historical',
             variable='ta',
             depth=100000,
             data_type='climatology',
             projection='zonal_mean',
             realm='atmos',
             ax=ax,
             )
    s = CompareLine(ifile='/raid/rc40/data/ncs/historical-cvu/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-1_historical-cvu_r1i1p1_185001-200512.nc',
             runID='cvu',
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='zonal_mean',
             realm='ocean',
             sister=q,
             ax=ax,
             labelled=True                      
             )
       
    r = Plot(runID='comparison',
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='global_map',
             realm='atmos',
             plotdata=s.plotdata - q.plotdata,
             dataset=s.dataset,
             rmse=True,
             divergent=True                  
             )
    
    s.plot()
    q.plot()
    r.plot()
    
    q = Plot(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID='edr',
             experiment='historical',
             variable='ta',
             depth=100000,
             data_type='climatology',
             projection='taylor',
             realm='atmos',)
    r = Taylor(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID= 'edr',
             unlabelledifiles = ['/raid/rc40/data/ncs/historical-cvu/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-1_historical-cvu_r1i1p1_185001-200512.nc'],
             experiment='historical',
             variable='ta',
             data_type='climatology',
             projection='taylor',
             realm='atmos',
             sister=q)
    r.plot()

    q = Histogram(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID='edr',
             experiment='historical',
             variable='ta',
             depth=100000,
             data_type='climatology',
             projection='taylor',
             realm='atmos',)
    r = HistogramCompare(ifile='/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc',
             runID= 'edr',
             labelledifiles = {'edr':'/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc'},
             unlabelledifiles = ['/raid/rc40/data/ncs/historical-cvu/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-1_historical-cvu_r1i1p1_185001-200512.nc','/raid/rc40/data/ncs/historical-edr/mon/atmos/ta/r1i1p1/ta_Amon_DevAM4-2_historical-edr_r1i1p1_185001-200012.nc'],
             experiment='historical',
             variable='ta',
             data_type='trends',
             projection='taylor',
             realm='atmos',
             sister=q)
    
    r.plot()
    """    
        
