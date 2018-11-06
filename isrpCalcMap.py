
from isrp import *

sensors = isrpLoadSensorParameters()

demFilename = 'dem/goms/goms.asc'
fig = plt.figure(num=1,figsize=(12, 8))

ax1 = fig.add_subplot(121)
xdem, ydem, zdem = isrpLoadDem(demFilename,sensors,1)
plt.contourf(xdem,ydem, zdem, cmap="gray",levels=list(range(0, 5000, 60)))
plt.title("Elevation Contours Goms (CH) area")
cbar = plt.colorbar(orientation="horizontal")
plt.gca().set_aspect('equal', adjustable='box')
lineMap, = plt.plot(np.nan,np.nan, linestyle='dashed', linewidth=1, color='b')

ax2 = fig.add_subplot(122)
lineProfile, = plt.plot(0, 0, linewidth=1, color='b')
lineSound, = plt.plot(0, 0,marker='o', markersize=6, linestyle='dashed', linewidth=1, color='r')
#ax2.set_xlim([0,20000])
ax2.set_ylim([1000,5000])
plt.grid(color='k', linestyle='-', linewidth=.1)
plt.draw()
plt.pause(0.1)

#dT=isrpDemTravelDt(demFilename,xdem,ydem,zdem,sensors,340,50)
loadT=np.load(demFilename+'dT.npZ')
dT=loadT['dT']
T=loadT['T']
overSampling = 1
S,dS=isrpDemTravelDs(-dT,T,50,overSampling)
corrM=np.array([[1, 1, 1, 1, 0, 0, 0, 0],
               [1, 1, 1, 1, 0, 0, 0, 0],
               [1, 1, 1, 1, 0, 0, 0, 0],
               [1, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 1, 1, 1],
               [0, 0, 0, 0, 1, 1, 1, 1],
               [0, 0, 0, 0, 1, 1, 1, 1],
               [0, 0, 0, 0, 1, 1, 1, 1]])

dMap=isrpArrange(demFilename,dS,S,sensors,25)

fig = plt.figure(num=2,figsize=(12, 8))
aa=dT[0,1,:,:].T
plt.subplot(2,2,1)
plt.imshow(aa)
aa=dT[0,2,:,:].T
plt.subplot(2,2,2)
plt.imshow(aa)
aa=dT[0,3,:,:].T
plt.subplot(2,2,3)
plt.imshow(aa)
aa=dT[1,2,:,:].T
plt.subplot(2,2,4)
plt.imshow(aa)
plt.colorbar()
plt.draw()
plt.show()
