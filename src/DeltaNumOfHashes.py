import numpy as np

def deltaHashNumPlot(cardAIntB, cardAUB, cont, conf, p, m, deltaStep, deltabeginpnt):
    delta = np.arange(deltabeginpnt, 1., deltaStep)
    kj = (-3 * np.log(conf/2) * cardAUB) / (np.sqrt(delta)*cardAIntB)
    kc = (-3*(cont+p)*np.log(conf/2))/np.sqrt(cont*delta)
    k_est = (-3 * (cont+p) * np.log(conf/2) * np.sqrt(cardAUB+(1+delta)*cardAIntB)) / (np.sqrt(cont*delta*cardAUB))
    return delta, kj, kc, k_est