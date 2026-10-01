"""PPT 용 공통 스타일."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def setup():
    plt.rcParams.update({
        'font.size': 13, 'axes.labelsize': 15, 'axes.titlesize': 16,
        'xtick.labelsize': 13, 'ytick.labelsize': 13, 'legend.fontsize': 12,
        'axes.linewidth': 1.4, 'xtick.major.width': 1.4, 'ytick.major.width': 1.4,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True,
        'figure.dpi': 110, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
        'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
    })

GREY = '#5a5a5a'
BLUE = '#1f5fd0'
