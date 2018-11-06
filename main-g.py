from correlogram import *
#import pandas as pd
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.colors import LightSource
from matplotlib.animation import FFMpegWriter

#plt.rcParams['animation.ffmpeg_path'] = '/Users/giacomo/Downloads/ffmpeg-20181018-f72b990-macos64-static/bin/ffmpeg'

def nanargmax(a):
    idx = np.argmax(a, axis=None)
    multi_idx = np.unravel_index(idx, a.shape)
    if np.isnan(a[multi_idx]):
        nan_count = np.sum(np.isnan(a))
        # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
        idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
        multi_idx = np.unravel_index(idx, a.shape)
    return multi_idx

smp=50
rsmp=1

asensors=[0, 1, 2, 3, 4, 5, 6, 7]
corrM=np.array([
               [0, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 1, 1, 1],
               [0, 0, 0, 0, 0, 0, 1, 1],
               [0, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0]])

# NETWORK PARAMETERS
demFilename = 'dem/goms/goms.asc'
sensors = isrpLoadSensorParameters()

# DEM DATA
xdem, ydem, zdem = isrpLoadDem(demFilename,sensors,1)
zdem[zdem == 0]=np.nan

# DT MATRIX from 3D-DEM
loadT=np.load(demFilename+'dMap.npZ')
dMap=loadT['dMap']
sMap=loadT['sMap']
minCorr=np.asscalar(loadT['minCorr'])
maxCorr=np.asscalar(loadT['maxCorr'])

#TEST DATA
# datafile = 'goms-20170114-060000-cav.mat' #750-800
# xObs = [437059, 437185, 437364, 436724, 437655, 437526, 437207, 436841]
# yObs = [5144416, 5145264, 5145114, 5145658, 5145842, 5146059, 5146222, 5144561]
# #start = 1620
# #stop = 1750
# fb = [1, 10]
#     # cotrolled 18 sec from Blitzigen
# start = 1480
# stop = 1880

datafile = 'goms-20170114-070000-cav.mat' # GazEx Reckigen and size2-3 cotrolled avalaches
xObs = [439815, 441130]
yObs = [5147174, 5148367]
start = 0
stop = 180
fb = [.5, 20]

# datafile = 'goms-20170305-050000-cav.mat' # WT Blitzingen and size2 cotrolled avalaches
# xObs = [437059, 437185, 437364, 436724, 437655, 437526, 437207, 436841]
# yObs = [5144416, 5145264, 5145114, 5145658, 5145842, 5146059, 5146222, 5144561]
# start = 550
# stop = 750
#fb = [1, 10]

# datafile = 'goms-20170305-050000-cav.mat' # WT Blitzingen and size2 cotrolled avalaches
# xObs = [437059, 437185, 437364, 436724, 437655, 437526, 437207, 436841]
# yObs = [5144416, 5145264, 5145114, 5145658, 5145842, 5146059, 5146222, 5144561]
# start = 550
# stop = 750
#fb = [1, 10]

# datafile = 'goms-20170117-130000-nav.mat' # GazEx Reckigen and size2-3 cotrolled avalaches
# xObs = []
# yObs = []
# start = 2800
# stop = 3050
# fb = [1, 5]

# datafile = 'goms-20170207-190000-nav.mat' # ??? treni d'onda
# xObs = []
# yObs = []
# start = 3200
# stop = 4000
# fb = [1, 5]

# datafile = 'goms-20170304-200000-nav.mat'
# xObs = []
# yObs = []
# start = 3200
# stop = 3400
# fb = [1, 10]

# datafile = 'goms-20170309-194800-nav-bigone.mat' #1060-1200
# xObs = 443310
# yObs = 5143823
# start = 1060
# stop = 1200
#fb = [1, 10]

# datafile = 'goms-20170306-032600-nav-wet.mat' #640-760 valanga verificata versante sud accanto a BLZ - 7 sensori (RCK#2 no) - elab ok "ma solo per coerenze alte"
# xObs = 439464
# yObs = 5141393
# start = 640
# stop = 760
#fb = [1, 10]#

# datafile = 'goms-20170309-194800-nav.mat' #175-280
# xObs = 0
# yObs = 0
#fb = [1, 10]

# datafile = 'goms-20170309-121800-nav.mat' #200-320
# xObs = 0
# yObs = 0
#fb = [1, 10]

#datafile = 'goms-20170309-100#800-plain.mat' #500-600

# datafile = 'goms-20170304-184800-nav.mat' #175-275
# xObs = np.nan
# yObs = np.nan
# start = 175
# stop = 275
#fb = [1, 10]

# datafile = 'goms-20170304-205400-nav.mat' #550-650
# xObs = np.nan
# yObs = np.nan
# start = 550
# stop = 650
#fb = [1, 10]

# da SE ... azimuth paralleli! segnale forte - isrp sembrao 2 aerei !!!
# datafile = 'goms-20170304-160000-nav.mat'
# xObs = []
# yObs = []
# start = 3400
# stop = 3750
# fb = [.8, 4]

# datafile = 'goms-20170304-225100-nav.mat' #335-460
# xObs = []
# yObs = []
# start = 310
# stop = 460
#fb = [1, 10]

# valanga verificata versante sud quasi nel mezzo - 7 sensori (RCK#2 no) - elab ok "ma solo per coerenze alte"
# datafile = 'goms-20170301-233500-nav.mat'
# xObs = 441798
# yObs = 5144189
# start = 40
# stop = 120
#fb = [1, 10]

# valanga verificata versante nord vicino BLZ - 8 sensori - elab ok
# datafile = 'goms-20170321-121000-nav-wet.mat'
# xObs = [437988]
# yObs = [5145034]
# start = 200
# stop = 380
#fb = [1, 10]

# datafile = 'goms-20170321-120000-nav-wet.mat' #220-340 valanga verificata versante nord vicino BLZ - 8 sensori - elab ok
# xObs = [437988]
# yObs = [5145034]
# start = 200
# stop = 380
#fb = [1, 10]

#metadata = dict(title='goms-20170305-050000-cav-reckigenGazEx', artist='G', comment='isrp!!')
#writer = FFMpegWriter(fps=10, metadata=metadata)

#datafile = 'goms-20170114-071500-cav.mat'
sgn=isrpGetMatLocalData('data/'+datafile)
# EVENT DATA SELECTION
start=start*smp*rsmp
stop=stop*smp*rsmp

sgn=sgn[:,start:stop]

i=0
j=3

defCorr=20*rsmp
l=np.arange(minCorr,maxCorr)
minCorrIdx=int(np.where(l>-defCorr)[0][0])
maxCorrIdx=int(np.where(l<defCorr)[0][-1])

sgnFlt = butter_bandpass_filter(sgn, fb[0], fb[1], smp, order=4)
time = np.linspace(0, sgnFlt.shape[1]/(smp*rsmp), sgnFlt.shape[1], endpoint=False)
sgnFlt, time =signal.resample(sgnFlt,sgnFlt.shape[1]*rsmp, time, axis=1)
time = np.linspace(0, sgnFlt.shape[1]/(smp*rsmp), sgnFlt.shape[1], endpoint=False)

shift=int(.5*smp*rsmp)
wnd=5*smp*rsmp
w=np.zeros([8,8,len(xdem),len(ydem)])
wc=np.zeros([8,len(xdem),len(ydem)])
kk=0
nI=np.sum(np.triu(corrM,1),axis=(0,1))
nD=len(np.where(np.sum(corrM, axis=(1)))[0])



# FIGURE
fig = plt.figure(num=1, figsize=(12, 8))
fig.suptitle(datafile, fontsize=16)
#fig.set_size_inches(10,8)
#ax1 = plt.subplot(2,1,1)
ax1 = fig.add_subplot(311)
ax1.set_position([.05, .3, .32, .6])



#sgnFlt[:,t:t+wnd]
#smax = np.nanmax([0, np.nanmax(sgnFlt[:, start:stop])])
smax = np.nanmax([0, np.nanmax(sgnFlt[:,500:-1])])
for i in range(0,8):
    #sigi = sgnFlt[i, start:stop]
    sigi = sgnFlt[i, 200:-1]
    sigi = sigi / np.max(sigi)
    if i in asensors:
        plt.plot(time[200:-1], sigi + 2 * i, linewidth=0.15, color='red')
        #plt.plot(time[start:stop], sigi + 2 * i, linewidth=0.15, color='red')
    else:
        plt.plot(time[200:-1], sigi + 2 * i, linewidth=0.1, color='black')
        #plt.plot(time[start:stop], sigi + 2 * i, linewidth=0.1, color='black')
#ax1.set_xlim([start/(smp*rsmp),stop/(smp*rsmp)])
ax1.set_xlim([0,len(time)/(smp*rsmp)])
slabels = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8']
plt.yticks(np.arange(0, 16, step=2), labels = slabels, size=8)

plt.ylabel('Sensors (NA)')
plt.xlabel('time (s)')



#trcsRck = plt.plot(time[start:stop], sgnFlt[4,start:stop], linewidth=0.2, color='blue')
startw, = plt.plot(0,0,linewidth=1, color='black')
ax1.grid(True)

#ax2 = plt.subplot(2,1,2)
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                      linewidth=0, antialiased=False)
xdem2, ydem2 = np.meshgrid(xdem, ydem)

#ax2 = fig.add_subplot(212, projection='3d')
ax2 = fig.add_subplot(312)
ax2.set_position([0.38, 0.3, .6, .6])

#rstride=1, cstride=1, linewidth=0
#ax2.plot_surface(xdem2, ydem2, zdem, cmap=cm.gray, rstride=2, cstride=2, linewidth=0, alpha=.5)
ax2.contour(xdem,ydem,zdem, cmap='gray', levels = list(range(1200, 4000, 25)), alpha=0.3)
ax2.plot(sensors['X'][3], sensors['Y'][3], marker='*', color='r', markersize=12)
ax2.annotate('IDA BLZ', (sensors['X'][3]-200, sensors['Y'][3]+250), color='red', fontsize=8,
             bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 2},
             verticalalignment='center')
ax2.plot(sensors['X'][7], sensors['Y'][7], marker='*', color='r', markersize=12)
ax2.annotate('IDA RCK', (sensors['X'][7], sensors['Y'][7]+250), color='red', fontsize=8,
             bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 2},
             verticalalignment='center')


cThr = 0
cm1 = zdem / np.nanmax(zdem)
cm1[cm1 < cThr] = np.nan

plt.plot(xObs,yObs,marker='<', color='green', markersize=8, ls='' ,alpha=.5)
# for i in range(0,len(xObs)):
#     ax2.annotate('Obs', (xObs[i]+200, yObs[i]), color='green', fontsize=8,
#                  bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 2},
#                  verticalalignment='center')


pSou, = ax2.plot(sensors['X'][7], sensors['Y'][7], marker='o', markersize=3, color='b', ls='')

cMap = ax2.contourf(xdem2,ydem2, cm1)

ax2.set_aspect("equal")
#ax2.set_yticklabels([])
#ax2.set_xticklabels([])
#ax2.grid(True)
#ax2.set_axis_off()
ax2.patch.set_facecolor('None')
hwide = 8000
ax2.set_xlim(440111-hwide, 440111+hwide)
ax2.set_ylim(5144500-hwide, 5144500+hwide)
ax2.annotate("", xy=(440111, 5144500-hwide+100), xytext=(440111+1000, 5144500-hwide+100),
             arrowprops=dict(arrowstyle="<->"), fontsize=14)
ax2.annotate('1 km', (440111+500, 5144500-hwide+300), color='black', fontsize=10,
             bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 1},
             verticalalignment='center', horizontalalignment='center')
#ax2.view_init(90,-90)
#ax2.dist = 6

axlogo = fig.add_subplot(111)
axlogo.set_position([0.04, 0.93, .17, .071])
img = plt.imread('imgs/GeCo.png')
axlogo.imshow(img, aspect='equal')
axlogo.set_axis_off()

ax3 = fig.add_subplot(313)
ax3.set_position([.1, .06, .8, .12])
pSoutim, = plt.plot([0, 0], [0, 0], marker='o', markersize=3, color='b', ls='')
plt.ylabel('Elevation (m)')
plt.xlabel('time (s)')
ax3.patch.set_facecolor('None')
#ax3.set_xlim(start/(smp*rsmp), stop/(smp*rsmp))
ax3.set_xlim([0,len(time)/(smp*rsmp)])
ax3.set_ylim(cThr, 1)
ax3.grid(True)

plt.draw()
plt.ion()
plt.show()


# REAL_TIME SYMULATION LOOP
print('REAL_TIME SYMULATION')
xsou = []
ysou = []
zsou = []
csou = []
tsou = []
nSensor=corrM.shape[0]
tSpan=60
tMap=np.zeros([tSpan,nSensor,len(xdem),len(ydem)],'float16')
tSpan2=int(tSpan/2)

CutDist=4500
iR = np.where((ydem > (sensors['Y'][3] - CutDist)) & (ydem < (sensors['Y'][7] + CutDist)))
iC = np.where((xdem > (sensors['X'][3] - CutDist)) & (xdem < (sensors['X'][7] + CutDist)))
xdem2R, ydem2R = np.meshgrid(xdem[iC], ydem[iR])
iCi, iRi = np.meshgrid(iR, iC)
zdemR = zdem[iRi, iCi]


#with writer.saving(fig, "goms-20170305-050000-cav-reckigenGazEx.mp4", 100):
for t in range(0, len(time)-shift, shift):
    tIns=int(t/shift) % tSpan
    for i in range(0,8):
        k=0

        for j in range(i+1,8):
            if(corrM[i,j]>0):
                k += 1
                w[i,j,:,:] ,m =isrpCorrelogram(sgnFlt[:,t:t+wnd], i, j, zdem, dMap,minCorrIdx,maxCorrIdx,defCorr)

        tMap[tIns, i, :, :] = np.sum(w, axis=(1))[i,:,:] / k
        tRead=tIns-tSpan2

    #wc[0,:,:]=tMap[tRead,0,:,:]
    for j in range(0, 8):
        wc[j, :, :] = isrpCorrelogram2(tMap, tIns, j, zdem, sMap, shift)

    # print(" tIns "+str(tIns)+" tRead "+str(tRead))
    #
    # print("tr "+str(tMap[tRead,0,100,100]))
    # print("ti "+str(tMap[tIns, 0, 100, 100]))
    wwS=np.nansum(wc,axis=(0))/nD
    print(int(np.nanmax(wwS) * 100))
    #wwe = wwS / np.nansum(wwS, axis=(0, 1))
    #print(int(np.nanmax(wwe) * 100))


    # Ritaglio

    wwc = wwS[iRi, iCi]
    print(np.max(wwc))
    wweR = wwc / np.nansum(wwc, axis=(0, 1))

    wwgood = wweR
    wwgood[np.where(wwc < 0.7)] = np.nan
    print(np.any(np.isfinite(wwgood)))

    startw.set_xdata([time[t], time[t]])
    startw.set_ydata([0, 8 * 2])

    if t/(smp*rsmp) > 35:
        if np.any(np.isfinite(wwgood)):
            print("LOCATION")

            for coll in cMap.collections:
                coll.remove()
            we = wwgood * zdemR.size / 2
            cMap = ax2.contourf(xdem2R, ydem2R, we.T*100, cmap='jet', levels = list(range(int(cThr*100), 100, 5)), alpha=0.9)

            pdmax = np.nanmax([0, np.nanmax(wwgood)])
            cmax= np.nanmax([0, np.nanmax(we)])

            wenodo = wwgood

            uu = np.where(we >= (np.nanmax(we) * .9))
            print(uu)
            xx = np.nanmedian(xdem[iC][uu[0]])
            yy = np.nanmedian(ydem[iR][uu[1]])
            zz = np.nanmedian(zdemR[uu])

            xsou = np.append(xsou, xx)
            ysou = np.append(ysou, yy)
            zsou = np.append(zsou, zz)
            csou = np.append(csou, cmax)
            tsou = np.append(tsou, t/(smp*rsmp))
            #
            pSoutim.set_xdata(tsou)
            pSoutim.set_ydata(csou)
            #
            pSou.set_xdata(xsou)
            pSou.set_ydata(ysou)


            plt.draw()


        plt.pause(0.001)

        #writer.grab_frame()


scipy.io.savemat(datafile, {'xsou': xsou, 'ysou': ysou, 'zsou': zsou, 'csou': csou, 'tsou': tsou})
# fig2 = plt.figure(num=2)
# ax2 = Axes3D(fig2)
# ax2.scatter(xsou, ysou, zsou, color='r')
#


