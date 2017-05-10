import numpy as np
import matplotlib.pyplot as plt


def deltaHashNumPlot(cardAIntB, cardAUB, cont, conf, p, m, deltaStep, deltabeginpnt):
    delta= np.arange(deltabeginpnt, 1., deltaStep)
    kj= (-3 * np.log(conf/2) * cardAUB) / (np.sqrt(delta)*cardAIntB)
    kc = (-3 * (cont+p) * np.log(conf/2) * np.sqrt(cardAUB+(1+delta)*cardAIntB)) / (np.sqrt(cont*delta*cardAUB))
    return delta,kj,kc


cardB = 1000000  # Cardinal of B
cardRatio = [10, 200]  # Ratio of cardinal of B over cardinal A
cont = [.1, .9]  # Containment of A in B
conf = .98  # Confident
p = 0.01  # False positive rate
deltaStep = 0.001  # Delta step
deltabeginpnt = 0.01

for m in cardRatio:
    for c in cont:
        cardA = int(cardB / m)  # computing card of A
        cardAIntB = int(c * cardA)  # computing card of intersection
        cardAUB = cardA + cardB - cardAIntB  # computing card of union

        delta, kj, kc = deltaHashNumPlot(cardAIntB, cardAUB, c, conf, p, m, deltaStep,deltabeginpnt)
        fig = plt.figure()
        plotTitle = 'Ratio = %d, Containment = %s' % (m, c)
        plt.plot(delta, kj, 'r--', label='Classic Min Hash')
        plt.plot(delta, kc, 'b-', label='Containment Min Hash')
        fig.suptitle(plotTitle, fontsize=20)
        plt.xlabel('Delta')
        plt.ylabel('Number of hash functions')
        plt.xlim([deltabeginpnt, 1])
        plt.legend()
        pltfigpath = '../Paper/Figs/%s%d.png' % (m,c*100)
        plt.savefig(pltfigpath)
        #plt.show()

print(kj[1])
print(kc[1])