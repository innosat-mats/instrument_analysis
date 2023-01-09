import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('/home/waves/projects/instrument_analysis/UtilityFunctions/src/retrieval/logs/time.txt')
a = a[a[:, 0].argsort()]

plt.scatter(a[:,0],a[:,1]/60)
plt.grid()
plt.xlabel('number of along track pts')
plt.ylabel('minutes')
plt.show()