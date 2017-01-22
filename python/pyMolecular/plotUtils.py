
import numpy as np
import matplotlib.pyplot as plt
from elements import ELEMENTS
from matplotlib import collections  as mc

def plotAtoms( es, xs, ys, scale=0.9, edge=True, ec='k', color='w' ):
    '''
    sizes   = [ ELEMENTS[ int(ei) ][7]*scale for ei in es ]
    colors  = [ '#%02x%02x%02x' %(ELEMENTS[ int(ei) ][8])  for ei in es ]
    print("sizes ",sizes)
    print("colors", colors)
    plt.scatter(xs,ys,s=sizes,c=colors)
    '''
    
    plt.fig = plt.gcf()
    ax = plt.fig.gca()
    plt.scatter(xs,ys)
    for i in range( len(es) ):
        element  = ELEMENTS[ int(es[i]) ]
        atomSize = element[7]*0.4
        fc = '#%02x%02x%02x' %element[8]
        if not edge:
            ec=fc
        circle=plt.Circle( ( xs[i], ys[i] ), atomSize, fc=fc, ec=ec  )
        ax.add_artist(circle)
    

def plotBonds( bonds, xs, ys, color='k', width=1 ):
    lines = [ ((xs[i],ys[i]),(xs[j],ys[j])) for (i,j) in bonds ]
    lc = mc.LineCollection( lines, colors=color, linewidths=width )
    plt.fig = plt.gcf()
    ax = plt.fig.gca()
    ax.add_collection(lc)
    '''
    for b in bonds:
        i=b[0]; j=b[1]
        lines = [ (  )])
        #plt.arrow( xs[i], ys[i], xs[j]-xs[i], xs[j]-ys[i], head_width=0.0, head_length=0.0,  fc='k', ec='k', lw= 1.0,ls='solid' )
        plt.draw.line(((xs[i],ys[i]),(xs[j],ys[j])), fill=color, width=width )
    '''


