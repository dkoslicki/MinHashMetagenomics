import numpy as np

def deltaconfident(cardAIntB, cardAUB, cont, k, p, deltaStep, deltabeginpnt):
    delta= np.arange(deltabeginpnt, 1., deltaStep)
    tj = np.exp((np.power(delta,2) * k*cardAIntB)/(-3*cardAUB))  # confident of Jaccard approach
    tc = np.exp((np.power(cont,2)*np.power(delta,2)*k)/(-3*(cont+p))) # confident of containment approach
    t_est = np.exp((np.power(cont,2)*np.power(delta,2)*np.power(cardAUB,2)*k)
                   /(-3*(cont+p)*np.power(cardAUB+(1+delta)*cardAIntB,2)))  # confident of estimated Jaccard with containment method
    return delta, tj, tc, t_est