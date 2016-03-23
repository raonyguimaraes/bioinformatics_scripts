#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! 
#fig = plt.figure()

#ax = Axes3D(fig) #<-- Note the difference from your original code...

#X, Y, Z = axes3d.get_test_data(0.05)
#cset = ax.contour(X, Y, Z, 16, extend3d=True)
#ax.clabel(cset, fontsize=9, inline=1)
#plt.show()

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
plt.show()
