import numpy as np

def RatioBoverANumofHashes(cardB, cont, conf, p, delta):
    m = np.arange(1,201,5)
    cardA = np.zeros(m.size)
    cardAIntB = np.zeros(m.size)
    cardAUB = np.zeros(m.size)
    cardA = cardB/m  # computing card of A
    cardAIntB = c * cardA  # computing card of intersection
    cardAUB = cardA + cardB - cardAIntB  # computing card of union
    kj = (-3 * np.log(conf/2) * cardAUB) / (np.sqrt(delta)*cardAIntB)
    kc = (-3*(cardAIntB/cardA+p)*np.log(conf/2))/np.sqrt(cardAIntB/cardA*delta)
    k_est = (-3 * (cont+p) * np.log(conf/2) * np.sqrt(cardAUB+(1+delta)*cardAIntB)) / (np.sqrt(cont*delta*cardAUB))
    return m, kj, kc, k_est


import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import os

# initial values
cardB = 1000  # Cardinal of B
cardRatio = [10, 200]  # Ratio of cardinal of B over cardinal A
cont = [.1, .9]  # Containment of A in B
conf = .01  # Confident
p = 0.01  # False positive rate
deltaStep = 0.001  # Delta step
deltabeginpnt = 0.0001
alpha_value = .6
k = 1000  # fix number of hash functions which we use
delta = 0.1
i=1
for c in cont:

    font_prop = {'family': 'serif', 'weight': 'normal', 'size': 14}

    # Relative error (delta)- Number of Hashes (k) Plot
    m, kj, kc, k_est = RatioBoverANumofHashes(cardB, c, conf, p, delta)
    fig = plt.figure()
    plt.semilogy(m, kj, 'r--', label='Classic Min Hash')
    #plt.semilogy(m, kc, 'b--', label='Containment Min Hash',alpha= alpha_value)
    #plt.semilogy(m,k_est, 'g--', label ='Jaccard Estimation by Containment Method',alpha=alpha_value)
    plt.semilogy(m, k_est, 'b--', label='Containment Min Hash', alpha=alpha_value)

    plt.xlabel(r"$|B|/|A|$",**font_prop)
    plt.ylabel('Number of hashes',**font_prop)
    plt.xlim([1, 200])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, prop= {'family': 'serif', 'weight': 'normal', 'size': 11}, ncol=2, mode="expand",
                   borderaxespad=0)

    # Save the final delta-k plot
    pltfigpath = '../Paper/Figs/increasingRatioWithCont%s.png' % (i)
    plt.savefig(pltfigpath)
    i=i+1
    #plt.show()
