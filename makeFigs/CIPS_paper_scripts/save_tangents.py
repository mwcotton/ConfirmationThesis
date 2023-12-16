import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import better_mob as bm
from datetime import datetime

savefile = datetime.now().strftime('%Y%m%d%H%M%S') + '_tangent.csv'
# savefile = '../tangents.csv'
print('Saving to file: ' + savefile)
enz_star=0.15
del_mu=8.0
del_e=-5.0
k_spo=1.0
k_cat=1.0
Dpe=4.0
Dse=1.0
v_rat=20.0
#make float with .0 incase ran in python2 (due to division)

dilute, dense = {}, {}
f = open(savefile, 'a')

#only if making a new file
f.write('del_mu,del_e,k_spo,k_cat,Dpe,Dse,v_rat,phiE1low,phiE1hi,phiE2low,phiE2hi\n')

# vars = np.arange(0.5, 2, 0.05)
# for k_cat in vars:
#     print(k_cat, end=' -- ')
#     tan1, tan2 = bm.find_binodals(1,0.05,1e-7,1e-7,del_mu,del_e,k_spo,k_cat,Dpe,Dse,v_rat)
#     dilute[k_cat] = (tan1[0]+tan1[1])/2
#     dense[k_cat] = (tan2[0]+tan2[1])/2
#     f.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(del_mu,del_e,k_spo,k_cat,Dpe,Dse,v_rat,tan1[0],tan1[1],tan2[0],tan2[1]) )
# f.close()

vars = np.arange(2, 10, 0.08)
for del_mu in vars:
    print(del_mu, end=' -- ')
    tan1, tan2 = bm.find_binodals(1,0.05,1e-7,1e-7,del_mu,del_e,k_spo,k_cat,Dpe,Dse,v_rat)
    dilute[del_mu] = (tan1[0]+tan1[1])/2
    dense[del_mu] = (tan2[0]+tan2[1])/2
    f.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(del_mu,del_e,k_spo,k_cat,Dpe,Dse,v_rat,tan1[0],tan1[1],tan2[0],tan2[1]) )
f.close()

print()

plt.scatter([dilute[elem] for elem in vars], vars, c='b', alpha=0.5)
plt.scatter([dense[elem] for elem in vars],  vars, c='k', alpha=0.5)
plt.ylabel(bm.labels['del_mu'])
plt.xlabel(bm.labels['enz_star'])
plt.draw()
plt.show()