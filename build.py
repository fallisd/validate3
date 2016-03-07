import classes

PLOT_OBJECTS = []
MAPS = ['global_map', 'mercator', 'polar_map', 'polar_map_south', 'zonal_mean']

def buildPlots(plots):
    for plot in plots:
        if depths in plot:
            for depth in plot['depths']:
                if plot['plot_projection'] == 'taylor':
                    buildTaylor(plot, depth)
                if plot['plot_projection'] in MAPS:
                    buildMap(plot, depth)
                if plot['plot_projection'] == 'section'
                    buildSection(plot, depth)
            
            if plot['compare_climatology']:
                buildCompare(plot, 'climatology')
            if plot['compare_trends']:
                buildCompare(plot, 'trends')

def buildTaylor():
    pass
    

def buildMap(plot, depth=None):
    plotlist = []
    plotlist.append(classes.Plot(**plot, depth=depth,
                                 pcolor_args=plot['data1_args']['pcolor_args'],
                                 ax_args=plot['data1_args']['ax_args']))
    
    for model in plot['comp_models']:
        plotobject = classes.Plot(**plot, depth=depth,
                                  pcolor_args=plot['data2_args']['pcolor_args'],
                                  ax_args=plot['data2_args']['ax_args']))
        plotobject.ifile = plot['model_file'][model]
        plotobject.runID = model
        plotlist.append(plotobject)
    
    for runid in plot['comp_ids']:
        plotobject = classes.Plot(**plot, depth=depth,
                                  pcolor_args=plot['data2_args']['pcolor_args'],
                                  ax_args=plot['data2_args']['ax_args']))
        plotobject.ifile = plot['id_file'][runid]
        plotobject.runID = model
        plotlist.append(plotobject)
    
    for cmip in plot['comp_cmips']:
        plotobject = classes.Plot(**plot, depth=depth,
                                  pcolor_args=plot['data2_args']['pcolor_args'],
                                  ax_args=plot['data2_args']['ax_args']))
        plotobject.ifile = plot['cmip_file'][cmip]
        plotobject.runID = cmip
        plotlist.append(plotobject)
    
    for obs in plot['comp_observations']:
        plotobject = classes.Plot(**plot, depth=depth,
                                  pcolor_args=plot['data2_args']['pcolor_args'],
                                  ax_args=plot['data2_args']['ax_args']))
        plotobject.ifile = plot['obs_file'][obs]
        plotobject.runID = obs
        plotlist.append(plotobject)        
    
    PLOT_OBJECTS.append(plotlist)       
