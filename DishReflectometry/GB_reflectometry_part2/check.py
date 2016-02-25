#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
	
	
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

a = np.cos(2*np.pi*np.linspace(0, 1, 60.))

ax1.plot(range(60), a)
ax2.plot(range(100), np.ones(100)) # Create a dummy plot
ax2.cla()
plt.show()