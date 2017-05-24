import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from VennPlot import vennplot
from DeltaNumOfHashes import deltaHashNumPlot
from DeltaConfident import deltaconfident
import os

# initial values
cardB = 100000  # Cardinal of B
cardRatio = [10, 200]  # Ratio of cardinal of B over cardinal A
cont = [.1, .9]  # Containment of A in B
conf = .98  # Confident
p = 0.01  # False positive rate
deltaStep = 0.001  # Delta step
deltabeginpnt = 0.0001
alpha_value = .6
k = 1000  # fix number of hash functions which we use

for m in cardRatio:
    for c in cont:

        cardA = int(cardB / m)  # computing card of A
        cardAIntB = c * cardA  # computing card of intersection
        cardAUB = cardA + cardB - cardAIntB  # computing card of union

        vennplot(cardAIntB, cardA, cardB / cardA, .2, alpha_value)  # make a Venn plot
        venn_file = '../Paper/Figs/Venn-AintB=%s-Ratio=%s.png' % (cardAIntB, int(cardB / cardA))  # save the Venn plot

        venn_image_file = cbook.get_sample_data(os.path.abspath(venn_file))  # load the Venn plot file
        venn_plot = plt.imread(venn_image_file)

        font_prop = {'family': 'serif', 'weight': 'normal', 'size': 14}

        # Relative error (delta)- Number of Hashes (k) Plot
        delta, kj, kc = deltaHashNumPlot(cardAIntB, cardAUB, c, conf, p, m, deltaStep, deltabeginpnt)
        fig = plt.figure()
        plt.semilogy(delta, kj, 'r--', label='Classic Min Hash')
        plt.semilogy(delta, kc, 'b-', label='Containment Min Hash')

        plt.xlabel('Relative error ($\delta$)',**font_prop)
        plt.ylabel('Number of hashes',**font_prop)
        plt.xlim([deltabeginpnt, 1])
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, prop= font_prop, ncol=2, mode="expand",
                   borderaxespad=0)
        # Add the Venn Graph to delta-k plot
        ax = plt.axes([.24, .26, .9, .9])
        ax.imshow(venn_plot)
        ax.axis('off')
        # Save the final delta-k plot
        pltfigpath = '../Paper/Figs/deltaK-%s%d.png' % (m,c*100)
        plt.savefig(pltfigpath)
        plt.close()

        # Relative error (delta)- Confident Plot
        delta, tj, tc, t_est = deltaconfident(cardAIntB, cardAUB, c, k, p, deltaStep, deltabeginpnt)
        fig = plt.figure()
        plt.plot(delta, tj, 'r--', label='Classic Min Hash')
        plt.plot(delta, tc, 'b--', label='Containment Min Hash')
        plt.plot(delta, t_est, 'g--', label='Jaccard Estimation by Containment Method')

        plt.axis([0, 1, 0, 1])
        plt.xlabel('Relative error ($\delta$)', **font_prop)
        plt.ylabel('Confidence', **font_prop)
        plt.xlim([deltabeginpnt, 1])
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, prop={'family': 'serif', 'weight': 'normal', 'size': 11}, ncol=2, mode="expand",
                   borderaxespad=0)
        # Add the Venn Graph to delta-k plot
        ax = plt.axes([.24, .26, .9, .9])
        #ax = plt.axes([.28, .42, .9, .65])
        ax.imshow(venn_plot)
        ax.axis('off')
        # Save the final delta-k plot
        pltfigpath = '../Paper/Figs/deltaConfident-%s%d.png' % (m, c * 100)
        plt.savefig(pltfigpath)
        plt.close()
