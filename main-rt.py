from correlogram import *
from isrp import *
#import pandas as pd
from scipy import interpolate
import time

    # STATION, SENSOR
station_id = 'T01'
sensors = isrpLoadSensorParameters(station_id)

    # STATION ELAB PARAMETRERS
elabPrms = isrpLoadElabDtsConf(station_id)

    # DEM DATA
#TODO define criteria for load dem data
xdem, ydem, zdem = isrpLoadDem(('dem/' + elabPrms[0]['dtmFilename']),sensors,1)
zdem[zdem == 0]=np.nan
# TDT matrix
#TODO define criteria for load dt matrix
loadT=np.load(('tdt/' + elabPrms[0]['dtFilename']+'dMap.npZ'))

    # DTS
dtsType = np.dtype({'names': ('time', 'ampAv', 'bkz', 'bkzSd', 'aVel', 'aVelSd',
                                          'cohAv', 'fqPk', 'fqMn', 'cnsAv', 'acSns', 'snr',
                                          'latSou', 'lonSou', 'eleSou', 'pdSou', 'station_id'),
                                'formats': ('float', 'float', 'float', 'float', 'float', 'float',
                                          'float', 'float', 'float', 'float', 'float', 'float',
                                            'float', 'float', 'float', 'float', 'str')})
#dts = np.zeros(0, dtype=dtsType)

    # REAL TIME LOOP PARAMETERS
pTimeLogFile = 'log/pTimeLogFile.txt'
if os.path.exists(pTimeLogFile):
    file = open(pTimeLogFile, 'r')
    pT1 = file.read()
    pT1 = datetime.strptime(pT1, '%Y-%m-%d %H:%M:%S')
else:
    pT1 = dt.datetime.utcnow()
    pT1 = pT1.replace(second=0)
    pT1 = pT1.replace(microsecond=0)
    file = open(pTimeLogFile, 'w')
    file.write(str(pT1))
    file.close()
#print('Last Processing Time: ' + str(pT1) + 'UTC')
offLineLag = 0 #86400 * 20
serverResposeTimeLag = 20 # [s] processig latency due to real time data trasmission and sharing latency
pLatency = 5 # algorith refresh time
pWindow = 30 # window analysis (seconds)
cTime = dt.datetime.utcnow()
#cTime = cTime.replace(second=0)
#cTime = cTime.replace(microsecond=0)
pT1 = cTime - dt.timedelta(seconds=+offLineLag)
pStartTime = cTime + dt.timedelta(seconds=+serverResposeTimeLag)

    # LOOP
while 1:
    cTime = dt.datetime.utcnow()
    if cTime > pStartTime:
        #tic = time.clock()
        # reading RealTime data from server
        pT0 = pT1 - dt.timedelta(seconds=+pWindow)
        print('... reading data: ' + str(pT0) + ' - ' + str(pT1))
        data, timeStamp = isrpGetWyData(sensors,pT0, pT1)

        if timeStamp[0] == 0:
            ss=0
            #print('!!! WARNING !!! No data from the server')
        else:
            # isrp Processing
            #print('... starting isrp processing')
            lastProcessedTime2 = datetime.utcfromtimestamp(timeStamp[-1] / 1000)
            lastProcessedTime1 = datetime.utcfromtimestamp(timeStamp[0] / 1000)
            print('... current time: ' + str(cTime))
            print('... readed data: ' + str(lastProcessedTime1) + ' - ' + str(lastProcessedTime2))

            # TODO  ... isrp signal to detection processing
            # dts = isrpS2D(elabPrms, sensors, xdem, ydem, zdem, loadT, timeStamp, data, dtsType)

            # result = isrpInsertElabDts(station_id, dts)



        pT1 = pT1 + dt.timedelta(seconds=+pLatency)

        alertLatency = dt.datetime.utcnow() - lastProcessedTime2 - dt.timedelta(seconds=+offLineLag)
        print('Alert letency: ' + str(alertLatency.seconds) + ' secs')

        pStartTime = pStartTime + dt.timedelta(seconds=+serverResposeTimeLag)
        #print(str('Start Processing Time: ' + str(pStartTime)))
       # toc = time.clock()
        #etime = int((toc - tic) * 1000)
        #print('#1 Loop cycle: ' + str(etime) + ' msecs')

    print('... waiting')
    time.sleep(1)














