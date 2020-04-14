from matplotlib import colors            
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe
import numpy as np
from obspy import UTCDateTime



def hourmap(data,
            bans=None,
            ax=None,
            scale = 1e9):
    width = data.index[1]-data.index[0]
    width = np.pi * 2 / 24 / 60 /60 * width.seconds 
    theta = np.asarray([(d.hour/24+d.minute/60/24)*np.pi*2-width/2 for d in data.index])
    radii = np.asarray([int(d.to_julian_date()+0.5) for d in data.index])
    radii = radii-min(radii)
    norm = colors.Normalize(vmin=np.nanpercentile(data,1),
                            vmax=np.nanpercentile(data,95))
    c_m = plt.cm.viridis
    s_m = plt.cm.ScalarMappable(cmap=c_m, 
                                norm=norm)
    s_m.set_array([])
    valid = np.where(np.isfinite(data))[0][::-1]
    
    if ax is None:
        ax=plt.figure(figsize=(7,9)).add_subplot(111, projection='polar')
    ax.grid(color='w',
            #path_effects=[pe.withStroke(linewidth=2,foreground='w')]
            )
    ax.set_xticks(np.linspace(0,np.pi*2*23/24,24))
    ax.set_xticklabels(['%dh'%h for h in range(24)])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rmax(max(radii))
    if bans is not None:
        rticks = [(UTCDateTime(ban).datetime - data.index.min().to_pydatetime()).days for iban,ban in enumerate(bans.keys())]
        xticks = [(UTCDateTime(ban).datetime.hour/24+UTCDateTime(ban).datetime.minute/60/24)*np.pi*2 for iban,ban in enumerate(bans.keys())]
        labels = [bans[iban] for iban in bans.keys()]
        xticks = [xticks[i] for i,d in enumerate(rticks) if d>0]
        labels = [labels[i] for i,d in enumerate(rticks) if d>0]
        rticks = [d for d in rticks if d>0]
        ax.set_rticks(rticks)
        for x,r,l,c in zip(xticks,
                           rticks,
                           labels,
                           range(len(labels))):
            ax.plot(x,r,'o',
                    label=l,
                    color='C%d'%c,
                    path_effects=[pe.withStroke(linewidth=5,
                                                foreground='w'),
                                  pe.withStroke(linewidth=3,
                                                foreground='k')])
    ax.set_yticklabels([])
    ax.set_rorigin(max(radii[valid])/-2)
    ax.text(np.pi,max(radii[valid])/-2,
            data.index[0].strftime("%Y-%m-%d"),
            ha='center',va='center')    
    ax.set_xlabel(data.index[-1].strftime("%Y-%m-%d"))
    plt.legend(loc='lower left',
               bbox_to_anchor= (0.0, -0.2), 
               ncol=2,
               borderaxespad=0, 
               frameon=False)
    cb=plt.colorbar(s_m,orientation='horizontal')#,pad=0.07)
    ticks = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x*scale))
    cb.ax.xaxis.set_major_formatter(ticks)
    cb.ax.set_xlabel("Displacement (nm)")    
    
    ax.bar(theta[valid], radii[valid], 
           color=s_m.to_rgba(data[valid]),
           bottom=radii[valid]-1,
           width=width)
    
    return ax

