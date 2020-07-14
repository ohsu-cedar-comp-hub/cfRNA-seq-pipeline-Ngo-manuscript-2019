import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib import rc
import numpy as np
import sys
from six import string_types
import statsmodels.nonparametric.api as smnp

rc('font', size = 5)
rc('font', family = 'arial')
rc('axes', labelsize = 8)


def _kde_support(data, bw, gridsize, cut, clip):
    """Establish support for a kernel density estimate."""
    support_min = max(data.min() - bw * cut, clip[0])
    support_max = min(data.max() + bw * cut, clip[1])
    return np.linspace(support_min, support_max, gridsize)


def _statsmodels_bivariate_kde(x, y, bw='scott', gridsize=100, cut=3, clip = [(-np.inf, np.inf), (-np.inf, np.inf)]):
    """Compute a bivariate kde using statsmodels."""
    if isinstance(bw, string_types):
        bw_func = getattr(smnp.bandwidths, "bw_" + bw)
        x_bw = bw_func(x)
        y_bw = bw_func(y)
        bw = [x_bw, y_bw]
    elif np.isscalar(bw):
        bw = [bw, bw]

    if isinstance(x, pd.Series):
        x = x.values
    if isinstance(y, pd.Series):
        y = y.values

    kde = smnp.KDEMultivariate([x, y], "cc", bw)
    x_support = _kde_support(x, kde.bw[0], gridsize, cut, clip[0])
    y_support = _kde_support(y, kde.bw[1], gridsize, cut, clip[1])
    xx, yy = np.meshgrid(x_support, y_support)
    z = kde.pdf([xx.ravel(), yy.ravel()]).reshape(xx.shape)
    return xx, yy, z

def plot_3d_scatter(lda_file='LDA_table_for_3D_plot.txt', elev=21, azim=56):
    """ Generate 3D scatter plot using LDA loadings

    Args:
        lda_file: precomputed LDA loading with columns : LD1, LD2, LD3, Groups, Colour

    """
    lda = pd.read_csv(lda_file,sep='\t',index_col=0)

    x = lda.LD1
    y = lda.LD2
    z = lda.LD3
    c = lda.Colour

    fig = plt.figure(figsize = (3.2, 3.2))
    ax1 = fig.add_subplot(111, projection = '3d')
    plt.subplots_adjust(left = 0, right = 0.90, top = 1, bottom = 0, wspace = 0.22)
    ax1.scatter(x,y,z,c=c,marker='o',**{'edgecolors':'white','alpha':1.0,'s':60})
    ax1.set_xlabel('\n LD1', size=8)
    ax1.set_ylabel('\n LD2', size=8)

    ax1.xaxis.set_rotate_label(False)
    ax1.yaxis.set_rotate_label(False)
    ax1.zaxis.set_rotate_label(False)

    ax1.view_init(elev = 10, azim = 135)
    # ax1.grid(False)
    ax1.xaxis.pane.set_edgecolor('black')
    ax1.yaxis.pane.set_edgecolor('black')
    ax1.zaxis.pane.set_edgecolor('black')

    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.yaxis.set_major_locator(MultipleLocator(5))

    [t.set_va('center') for t in ax1.get_yticklabels()]
    [t.set_ha('right') for t in ax1.get_yticklabels()]

    [t.set_va('center') for t in ax1.get_xticklabels()]
    [t.set_ha('left') for t in ax1.get_xticklabels()]

    [t.set_va('center') for t in ax1.get_zticklabels()]
    [t.set_ha('right') for t in ax1.get_zticklabels()]

    ax1.xaxis._axinfo['tick']['inward_factor'] = 0
    ax1.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax1.yaxis._axinfo['tick']['inward_factor'] = 0
    ax1.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax1.zaxis._axinfo['tick']['inward_factor'] = 0
    ax1.zaxis._axinfo['tick']['outward_factor'] = 0.4
    ax1.zaxis._axinfo['tick']['outward_factor'] = 0.4

    color_dict  = {'#000000':'Greys','#9000FF':'Purples','#00B700':'Greens','#FF010A':'Reds'}
    lda['cmap']=lda.Colour.map(color_dict)

    groups = lda.groupby('Groups')
    yoffset = ax1.get_w_lims()[2]
    xoffset = ax1.get_w_lims()[4]

    for group,frame in groups:
        cmap = frame.cmap.unique()[0]
        xx,yy,zz = _statsmodels_bivariate_kde(frame.LD1, frame.LD2)
        ax1.contour(xx, yy, zz, zdir = 'z', cmap = cmap, offset=-4, linewidths=0.5)
        #ax1.contour(xx, yy, zz, zdir = 'x', offset = -5.973452472012009, cmap = cmap, linewidths=0.5)
        #ax1.contour(xx, yy, zz, zdir = 'y', offset = -6.5, cmap = cmap, linewidths=0.5)
    if elev and azim != None:
        ax1.view_init(elev=elev, azim=azim)

    #plt.show()

    plt.savefig('test_3d_small_size.png',dpi=500)


if __name__ == '__main__':
    plot_3d_scatter()
