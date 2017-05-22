import matplotlib.pyplot as plt
import math
import shapely.geometry as sg
import descartes

def vennplot(AintB,cardA,RatioBoverA,radB,alpha_value):
    if AintB > cardA:
        AintB = cardA

    radA = math.sqrt(1/RatioBoverA)*radB
    xpntA= .5+radB+(1-2*(AintB/cardA))*radA
    B = sg.Point(.5, .5).buffer(radB)
    A = sg.Point(xpntA, .5).buffer(radA)

    left = A.difference(B)
    right = B.difference(A)
    middle = A.intersection(B)

    ax = plt.gca()
    ax.add_patch(descartes.PolygonPatch(left, fc='b', ec='k', alpha=alpha_value))
    ax.add_patch(descartes.PolygonPatch(right, fc='r', ec='k', alpha=alpha_value))
    ax.add_patch(descartes.PolygonPatch(middle, fc='g', ec='k', alpha=alpha_value))

    font = {'family': 'serif', 'color':  'black', 'weight': 'normal', 'size': 18}

    ax.text(0.5,.5+radB , r'$B$', horizontalalignment='left', verticalalignment='bottom', fontdict=font, usetex=True)
    ax.text(xpntA+radA, 0.5, r'$A$', horizontalalignment='left', verticalalignment='bottom', fontdict=font, usetex=True)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_aspect('equal')
    plt.xticks([])
    plt.yticks([])
    pltfigpath = '../Paper/Figs/Venn-AintB=%s-Ratio=%s.png' % (AintB,int(RatioBoverA))
    plt.savefig(pltfigpath, transparent=True)
    plt.close()
    #plt.show()
