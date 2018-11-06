# Earthworm lib
from obspy import Stream, UTCDateTime
from obspy.clients.earthworm import Client
import numpy as np


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def isrpGetWsData(StzPrms,Ti,Te):
    # ti_str = "20180225_11:00:00"
    # tf_str = "20180225_12:00:00"

    ti = UTCDateTime(Ti)
    tf = UTCDateTime(Te)

    client = Client(StzPrms['WssClient'][0], StzPrms['WssPort'][0], timeout=20)

    Wfrm = Stream()
    for i in range(0,StzPrms['WssChannel'][0]):
        response = client.get_availability(StzPrms['WssNetwork'][0], StzPrms['StationName'][0],
                                           StzPrms['WssLocation'][0], 'CH'+str(i+1))
        print(response)
        out = client.get_waveforms(StzPrms['WssNetwork'][0], StzPrms['StationName'][0], StzPrms['WssLocation'][0], 'CH'+str(i+1), ti, tf)
        Wfrm += out

    return Wfrm


StzPrmsType = np.dtype({'names': ('WssClient', 'WssPort', 'WssNetwork', 'WssLocation', 'StationName', 'WssChannel'),
                          'formats': ('U16', 'uint16', 'U2', 'U2', 'U3', 'uint8')})
StzPrms = np.zeros(1, dtype=StzPrmsType)

StzPrms['WssClient'] = '148.251.122.130'
StzPrms['WssPort'] = 8081
StzPrms['WssNetwork'] = 'FI'
StzPrms['WssLocation'] = '00'
StzPrms['WssChannel'] = 5
StzPrms['StationName'] = 'NO1'

Ti = "20181025_22:00:00"
Te = "20181025_22:01:00"

Wfrm = isrpGetWsData(StzPrms,Ti,Te)

npts = Wfrm.traces[0].stats.npts
samprate = Wfrm.traces[0].stats.sampling_rate
tstrt = Wfrm.traces[0].stats.starttime.timestamp  # >> NB: UTCDateTime(tstrt)=tr.stats.starttime
UnixTimeStamp = np.arange(tstrt, npts / samprate + tstrt, 1 / samprate)  # >> unix timestamp ()


fig = plt.figure(num=1)

for i in range(0,StzPrms['WssChannel'][0]):

    Sig = Wfrm.traces[i].data
    plt.subplot(StzPrms['WssChannel'][0],1,i+1)
    plt.plot(UnixTimeStamp,Sig)





