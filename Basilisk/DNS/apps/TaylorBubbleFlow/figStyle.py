import numpy as np

markEv = 100
markedLineTypes = ['-o', '-s', '-d', '-^', '-v']
markedLineMarkers = ['o', 's', 'd', '^', 'v']

colors = ['dodgerblue', 'aqua', 'gold', 'limegreen', 'orangered', 'firebrick', 'black']

#colors4 = ['dodgerblue', 'limegreen', 'gold', 'orangered']
colors4 = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']

def prep(p, size):

    p.rc('font', family='serif', size=9, serif='STIXGeneral')
    p.rc('mathtext', fontset='stix')
    p.rc('axes', labelsize=9, titlesize=9)
    p.rc('lines', dash_capstyle='round', linewidth=1, color='black', markersize=3)

    p.rc('xtick.major', size=4, width=0.25, pad=2)
    p.rc('xtick.minor', size=2, width=0.25, pad=2)
    p.rc('ytick.major', size=4, width=0.25, pad=2)
    p.rc('ytick.minor', size=2, width=0.25, pad=2)

    p.rc('legend', fontsize=6, handlelength=2, numpoints=1, labelspacing=0.2, borderpad=0.2, fancybox=False, )

    p.rc('savefig', format='pdf')

    p.rc('figure', max_open_warning=0)

    p.autoscale(tight=True)

    if size == 'square':
        p.rc('figure', figsize=(2.75, 2.65))

    elif size == 'column':
        p.rc('figure', figsize=(1.5, 6))

    elif size == 'Iguchi':
        p.rc('figure', figsize=(2.5, 4))


def post(f, size, l=False):

    ax = f.gca()

    from matplotlib.ticker import ScalarFormatter

    for pos in ['bottom', 'top', 'right', 'left']:
        ax.spines[pos].set_linewidth(0.5)

    if l:
        l.get_frame().set_linewidth(0.25)
        l.get_frame().set_edgecolor('black')

    figSize = f.get_size_inches()

    if size == 'square':
        ax.set_position([0.25, 0.15, 0.72, 0.72*figSize[0]/figSize[1]])
    elif size == 'column':
        ax.set_position([0.32, 0.07, 0.6, 0.6*6*figSize[0]/figSize[1]])
    elif size == 'Iguchi':
        ax.set_position([0.23, 0.1, 0.71, 0.89])

    if ax.get_yscale() != 'log':
        fmt = ScalarFormatter()
        fmt.set_powerlimits((-2,2))
        ax.yaxis.set_major_formatter(fmt)

    return

def markevery(i, j):

    return (int(np.mod(i*markEv/j,markEv)), markEv)

def noLines(cnt):
    for c in cnt.collections:
        c.set_edgecolor("face")

