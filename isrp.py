from correlogram import *
import os
import json
import numpy as np
from numpy import hstack
from numpy import vstack
from scipy import signal
from scipy.signal import butter, lfilter
import scipy.io as sio
import scipy.ndimage
#import scilab2py

# Earthworm lib
from obspy import Stream, UTCDateTime
#from obspy.clients.earthworm import Client

import urllib.request
from urllib.parse import unquote
import datetime as dt
import time
from datetime import datetime

import math
from math import cos, pi, sqrt, floor

import utm

import imp

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import fnmatch

import requests
import wget

from collections import namedtuple

from joblib import Parallel, delayed
import multiprocessing

import pymysql.cursors



sensorsType = np.dtype({'names': ('name', 'configFilename', 'X', 'Y', 'Z'),
                          'formats': ('U10', 'U60', 'f8', 'f8', 'f8')})

#""" Convert python datetime object to matlab serial number """
def datetime2matlab(t_datetime):
    mdn = t_datetime + dt.timedelta(days = 366)
    frac_seconds = (t_datetime - dt.datetime(t_datetime.year,t_datetime.month,t_datetime.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    frac_microseconds = t_datetime.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    t_matlab = mdn.toordinal() + frac_seconds + frac_microseconds
    return t_matlab


#------------------------ LOAD INPUTS from configfile
#print '... loading input parameters'
def isrpLoadSensorParameters():
    modname = 'stationparameters'
    network = 'goms'
    path_config = os.curdir+'/configfiles/'+network+'/'
    #file_config = glob.glob(path_config+'*.txt')
    listOfFiles = os.listdir(path_config)
    pattern = "*.txt"

    i = 0
    file_config = {}
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            file_config[i] = path_config + entry
            i += 1

    nSensors=len(file_config)
    sensors=np.zeros(nSensors, dtype=sensorsType)
    i = 0
    for entry in file_config:
        print(file_config[entry])
        I = imp.load_source(modname, file_config[entry])
        sensors['name'][i] = I.id
        sensors['configFilename'] = file_config[entry]
        sensors['X'][i] = I.x
        sensors['Y'][i] = I.y
        sensors['Z'][i] = I.z
        i += 1
    return sensors


def isrpLoadDem(demFilename,sensors,undersamplingFactor):
    f = open(demFilename, 'r')
    whl = f.read()
    whl = whl.split('\n')
    # header
    hdr = whl[0:6]
    deminfo = np.dtype({'names': ('ncols', 'nrows', 'xllcorner', 'yllcorner', 'dx', 'dy'),
                        'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
    nInfo = len(hdr)
    info = np.zeros(1, dtype=deminfo)
    s = hdr[0].split()
    info['ncols'] = float(s[1])
    s = hdr[1].split()
    info['nrows'] = float(s[1])
    s = hdr[2].split()
    info['xllcorner'] = float(s[1])
    s = hdr[3].split()
    info['yllcorner'] = float(s[1])
    s = hdr[4].split()
    info['dx'] = float(s[1])
    s = hdr[5].split()
    info['dy'] = float(s[1])

    xllcu, yllcu, zonen, zonel = utm.from_latlon(info['yllcorner'], info['xllcorner'])
    print('utm zone: ' + str(zonen) + zonel)

    # Latitude:  1 deg = 110.54 km
    # Longitude: 1 deg = 111.320*cos(latitude) km
    dy = info['dy'] * 110.54 * 1000
    dx = info['dx'] * 111.32 * cos(info['yllcorner'] * pi / 180) * 1000

    xu = np.arange(xllcu, xllcu + info['ncols'] * dx, dx)
    #print(info['ncols'])
    yu = np.arange(yllcu + info['nrows'] * dy, yllcu, -dy)
    #print(info['nrows'])

    # elevation matrix
    bb = whl[7:len(whl)]
    print(len(bb))
    idx = 0
    for i in range(0, len(yu)):
        bw = bb[i].split(' ')
        bw.remove('')
        aa = np.asarray(bw, dtype=np.float32)
        if idx == 0:
            z = aa
        else:
            z = vstack((z, aa))
        idx += 1

    z[z < 0] = 0

    # TODO insert resampling factor as fuction input
    z=z[0::undersamplingFactor, 0::undersamplingFactor]
    xu=xu[0::undersamplingFactor]
    yu=yu[0::undersamplingFactor]
    #f.close()

    print('ncols: ', str(len(xu)))
    print('nrows: ', str(len(yu)))
    print('zmatrix: ' + z.shape[0].__str__() + ' x ' + z.shape[1].__str__())

    # int(np.amax(z))
    #plt.contourf(xu, yu, z, cmap="gray",
     #            levels=list(range(0, 5000, 50)))
    #plt.title("Elevation Contours Goms (CH) area")
    #cbar = plt.colorbar()
    #plt.gca().set_aspect('equal', adjustable='box')
    # lat=5146098.00 m N
    # lon=442227.00 m E
    print(len(sensors))
    #for i in range(0,len(sensors)):
    #    print(i)
    #    plt.plot(sensors['X'][i], sensors['Y'][i], marker='o', markersize=2, color='r', ls='')

    #plt.show()
    return xu, yu, z

#TODO a serius travel time function
def isrpTravelTimesComput(x0, y0, xS, yS, x1, y1, xdem ,ydem ,zdem, c, demres):
    # x0 = 13  node xposition index
    # x1 = 255 sensor xposition index
    # y0 = 135 node yposition index
    # y1 = 222 sensor yposition index
    # demres resolution
    #fig2 = plt.figure(num=1,figsize=(12, 8))

    zdem[np.isnan(zdem)]=0


    dPlane = (sqrt((xS - xdem[x0]) ** 2 + (yS - ydem[y0]) ** 2))

    if(dPlane>100):
        num = floor(dPlane/demres)

        xp, yp = np.linspace(x0, x1, num), np.linspace(y0, y1, num)
        zp = scipy.ndimage.map_coordinates(zdem, np.vstack((yp, xp)))
        dr = dPlane/(len(zp)-1)
        # (len(zp) - 1)
        cr = dr * np.arange(0, len(zp), dtype=np.int)
        nd = len(cr)

        mm = (zp[-1] - zp[0]) / cr[-1]  # diff quota / distanza (coeff ang o elevation)
        r = mm * cr + zp[0]  # retta con coeff. ang pari a elevation tra array e nodo
        r.transpose()
        h = np.fix((zp - r) / 10) * 10  # 10 = z resolution (m)
        h[np.nonzero(h < 0)] = 0
        h.transpose()

        k = 0
        io = 0
        ii=[]
        xx=[]
        yy=[]

        while io < nd - 1:
            for i in range(io, nd):

                m1 = (h[i] - h[io]) / (cr[i] - cr[io])
                n1 = (h[i] * cr[io] - h[io] * cr[i]) / (cr[io] - cr[i])
                r1 = m1 * cr[io:i] + n1
                h1 = np.fix((h[io:i] - r1) / 10) * 10
                ih = h1 > 0
                ih = ih.astype(np.int)
                if np.sum(ih) == 0:
                    ik = i

                ii.append(ik)
                xx.append(cr[io])
                yy.append(zp[io])
                k += 1

            io = ik

        xx = np.append(xx, cr[-1])
        yy = np.append(yy, zp[-1])

        dTopo = 0
        for i in range(1, k + 1):
            dTopo = dTopo + sqrt((xx[i] - xx[i - 1]) ** 2 + (yy[i] - yy[i - 1]) ** 2)
        tT = dTopo/c

        # lineMap.set_xdata([xu[x0] , xS])
        # lineMap.set_ydata([yu[y0] , yS])
        #
        # lineProfile.set_xdata(cr)
        # lineProfile.set_ydata(zp)
        #
        # lineSound.set_xdata(xx)
        # lineSound.set_ydata(yy)
        # plt.sca(ax2)
        # plt.gca().relim()
        # plt.gca().autoscale_view()
        # plt.pause(0.05)

    else:
        tT=np.nan
        dTopo=np.nan
    # print("dPlane"+str(dPlane))
    #print("dtopo "+str(dTopo))
    return tT


def isrpParallelDemTravelDt(i,xdem, ydem, zdem, sensors:sensorsType, c, demRes):
    tT = np.zeros(( len(xdem), len(ydem)))

    dx = xdem[1] - xdem[0]
    dy = ydem[1] - ydem[0]

    x1 = np.where(xdem > sensors['X'][i])
    y1 = np.where(ydem < sensors['Y'][i])

    x1 = x1[0][0]
    y1 = y1[0][0]
    xS = sensors['X'][i]
    yS = sensors['Y'][i]

    dxS = (xdem[x1] - xS) / dx
    dyS = (ydem[y1] - yS) / dy

    x1 = x1 + dxS
    y1 = y1 + dyS
    for xi in range(0,len(xdem)):
        print("s " +str(i) +"xi "+str(xi))
        for yi in range(0,len(ydem)):
            #if not np.isnan(zdem[yi, xi]):
            if (zdem[yi, xi]>0):
                tT[xi, yi] = isrpTravelTimesComput(xi, yi,xS, yS, x1, y1,xdem,ydem,zdem, c,demRes)
            else :
                tT[xi, yi] = np.nan

   # T[i,:,:]+=tT[:,:]
    #T.put(i,tT)
    return tT


def isrpDemTravelDt(demFilename, xdem, ydem, zdem, sensors:sensorsType,c,demRes):

    #xdem, ydem, zdem, = isrpLoadDem(demFilename,sensors,1)
    #TODO esclude sensor matrix
    dT=np.zeros((len(sensors), len(sensors), len(xdem), len(ydem)))
    T = np.zeros((len(sensors), len(xdem), len(ydem)))

    # for i in range(0,len(sensors)):
    #     print( "isrpDemTravelDt dElaborating sernsor "+str(i))
    #
    #     dx = xdem[1] - xdem[0]
    #     dy = ydem[1] - ydem[0]
    #
    #     x1 = np.where(xdem > sensors['X'][i])
    #     y1 = np.where(ydem < sensors['Y'][i])
    #
    #     x1 = x1[0][0]
    #     y1 = y1[0][0]
    #     xS = sensors['X'][i]
    #     yS = sensors['Y'][i]
    #
    #     dxS = (xdem[x1] - xS) / dx
    #     dyS = (ydem[y1] - yS) / dy
    #
    #     x1 = x1 + dxS
    #     y1 = y1 + dyS
    #
    #     for xi in range(0,len(xdem)):
    #         for yi in range(0,len(ydem)):
    #             print("s " +str(i) +"xi "+str(xi)+"  yi "+str(yi))
    #             #if not np.isnan(zdem[yi, xi]):
    #             if (zdem[yi, xi]>0):
    #                 T[i, xi, yi] = isrpTravelTimesComput(xi, yi,xS, yS, x1, y1,xdem,ydem,zdem, c,demRes)
    #             else :
    #                 T[i, xi, yi] = np.nan
    #           #  print(T[i, xi, yi])

    num_cores = multiprocessing.cpu_count()

    T[:]=Parallel(n_jobs=num_cores)(delayed(isrpParallelDemTravelDt)(i,xdem, ydem, zdem, sensors, c, demRes) for i in range(0,len(sensors)))
    print(T)
    for i in range(0,len(sensors)):
        for ii in range(i+1, len(sensors)):
            print( "isrpDemTravelDt dElaborating sernsor couple "+str(i)+" "+str(ii))
            for xi in range(0,len(xdem)):
                for yi in range(0,len(ydem)):
                    dT[i, ii, xi, yi] = T[ii,xi,yi]-T[i,xi,yi]
                   # print(dT[i,ii, xi, yi])


                    #dT[ii, i, xi, yi] = -dT[i, ii, xi, yi]
    np.savez(demFilename+'dT',dT=dT,T=T)
    return dT

def isrpDemTravelDs(dT:np.array, T, smp,overSampling):
    dSample = dT*smp*overSampling
    Sample = T * smp * overSampling
    #dSample = dSample.astype(int)
    return Sample, dSample

def isrpArrange (demFileName,dS:np.array,S,sensors,sShift):
    dMap=[]
    sMap=[]
    dsMin = np.nanmin(dS)
    dsMax = np.nanmax(dS)#np.where(corrM>0)
    for i in range(0,len(sensors)):
        dii=[]
        for ii in range(0, len(sensors)):
            #print("isprArrange dElaborating sernsor couple " + str(i) + " " + str(ii))
            dj=[]
            for j in range(int(dsMin-1),int(dsMax+1)):
              #  print("isprArrange Elaborating T " + str(j))
                dj.append(np.where((dS[i,ii,:,:] >= j) & (dS[i,ii,:,:]< j+1)))
            dii.append(dj)
        dMap.append(dii)

    for i in range(0, len(sensors)):
        dj = []
        for j in range(0, int(np.nanmax(S)),sShift):
            print("isprArrange Elaborating T " + str(j))
            dj.append(np.where((S[i, :, :] >= j) & (S[i, :, :] < j + sShift)))
        sMap.append(dj)


    np.savez(demFileName + 'dMap', dMap=dMap,sMap=sMap,minCorr=dsMin,maxCorr=dsMax)
    return dMap, sMap

def isrpGetWsData(sensors,ti_str,tf_str):
    # ti_str = "20180225_11:00:00"
    # tf_str = "20180225_12:00:00"

    ti = UTCDateTime(ti_str)
    tf = UTCDateTime(tf_str)

    # response = client.get_availability('FI', 'NO1')
    # print(response)

    st = Stream()
    idx=0
    for i in sensors:
        try:
            I = imp.load_source('stationparameters', i['configFilename'])
            client = Client(I.wsclient, I.wsport, timeout=20)
            a = client.get_waveforms(I.wsnetwork, I.stationame, I.wslocation, I.wschannel, ti, tf)
            # ,metadata=True)
            # NB: last character of channel can be a wildcard ('?' or '*') to fetch 'Z', 'N', and 'E' component.
        except:
            print("    > channel  failed")
            continue
        st += a
        if idx == 0:
            data = a
        else:
            data = vstack((data, a.traces[0].data))
        idx +=1

        if len(st) == 0:
            print('... WARNING: NO DATA retrieved from WWS >>> quitting !')
            quit()

    npts = st.traces[0].stats.npts
    samprate = st.traces[0].stats.sampling_rate
    tstrt = st.traces[0].stats.starttime.timestamp  # >> NB: UTCDateTime(tstrt)=tr.stats.starttime
    timestamps = np.arange(tstrt, npts / samprate + tstrt, 1 / samprate)  # >> unix timestamp ()

    print(data.shape)
    return data, timestamps

def isrpGetWyData(sensors,ti_str,tf_str):
    tmin = ti_str.strftime('%Y-%m-%d%%20%H:%M:%S.%F')
    tmax = tf_str.strftime('%Y-%m-%d%%20%H:%M:%S.%F')
    idx = 0
    data = np.empty([0])
    for i in range(0,len(sensors)):
        try:
            domain = sensors[i]['serverAdrs']
            location = sensors[i]['serverApi']
            key = sensors[i]['serverKey']
            station = sensors[i]['id']
            args = "key={key}&id={id}&json=true&tmin={tmin}&tmax={tmax}".format(id=station, tmin=tmin, tmax=tmax,
                                                                                key=key)
            req = "{d}/{p}?{args}".format(d=domain, p=location, args=args)
            #print(req)
            #TODO check server running and send ServerStatusFlag
            # tic = time.clock()
            with urllib.request.urlopen(req) as url:
                data = json.loads(url.read().decode())
            # toc = time.clock()
            # etime = int((toc - tic) * 1000)
            # print('Server Response Time ID ' + str(sensors[i]['id']) + ': ' + str(etime) + ' msecs')
            out = np.array(data["values"])
            a = out[:, 1]

            if out.size == 0:
                idx += 1

        except:
            print("    > ID " + str(sensors[i]['id']) + " failed")

            continue

        if idx == 0:
            data = a
        else:
            data = vstack((data, a))

        #out.tt = datenum(1970, 1, 1) + time'/86400;
    #[ref_time + datetime.timedelta(seconds=i) for i in secs]
    print(data)
    if data.size == 0:
        timestamp = 0
        print('!!! W A R N I N G !!! ... Empty Data')
    else:

        secs = out[:, 0]
        #timestamp = [dt.datetime(1970, 1, 1) + datetime.timedelta(seconds=i) for i in secs]
        timestamp = secs

    return data, timestamp


def isrpGetMatData(sensors,ti_str,tf_str):
    #tmin.strftime('%Y-%m-%d%%20%H:%M:00'), tmax.strftime('%Y-%m-%d%%20%H:%M:%S')
    idx=0
    for i in sensors:
        try:
            I = imp.load_source('stationparameters', i['configFilename'])
            filename = I.stationame+'_'+ti_str[0:]+'.mat'
            url = I.matclient+'/'+I.matlocation+'/'+I.stationame+'/'+ti_str[0:4]+'/'+I.stationame+'_'+ti_str[0:8]+'/'+filename
            url = url.replace(':', '')
            url = 'http://'+url
            #http://148.251.122.130/matfiles/GMS/2017/GMS_20170309/GMS_20170309_010000.mat
            # Load the data
            filedwnl = wget.download(url,'temp.mat')
            dat = sio.loadmat(filedwnl)
            out = dat.get('data',['values'])
            fields = out.dtype.names
            a = out[I.channel]
            a1 = a[0][0].transpose()
            #a1 = a1.transpose()

            os.remove('temp.mat')

        except:
            print("    > channel " + i + " failed")
            continue

        if idx == 0:
            data = a1
        else:
            data = vstack((data, a1))
        idx += 1

    print(data.shape)
    timestamp = out['tt']

    return data, timestamp

def isrpGetMatLocalData(fileName):
    dat = sio.loadmat(fileName)
    out = dat.get('gomsdata',['values'])
   # fields = out.dtype.names
 #  a = out[I.channel]
   # a1 = a[0][0].transpose()

    return out

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def isrpLoadSensorParameters(station_id):
    # Connect to the database
    conn = pymysql.connect(host='85.10.202.61',
                           user='geco',
                           password='geco-company-2018',
                           db='GeCoMonit',
                           cursorclass=pymysql.cursors.DictCursor)
    try:

        with conn.cursor() as cursor:
            # Read a single record
            sql = "SELECT * FROM `sensors` WHERE `station_id` = " + "'" + station_id + "'"
            out = cursor.execute(sql)
            out = cursor.fetchall()
            sensorsType = np.dtype({'names': ('smp', 'id', 'gain', 'type', 'lon', 'lat', 'ele', 'statio_id'),
                                    'formats': ('uint8', 'uint8', 'f16', 'U32', 'f8', 'f8', 'uint8', 'U3')})
            nSensors = len(out)
            sensors = np.zeros(nSensors, dtype=sensorsType)
            sensors = out
    finally:
        conn.close()

    return sensors

def isrpLoadElabDtsConf(station_id):
    # Connect to the database
    conn = pymysql.connect(host='85.10.202.61',
                           user='geco',
                           password='geco-company-2018',
                           db='GeCoMonit',
                           cursorclass=pymysql.cursors.DictCursor)
    try:

        with conn.cursor() as cursor:
            # Read a single record
            sql = "SELECT * FROM `elabDtsConf` WHERE `station_id` = " + "'" + station_id + "'"
            out = cursor.execute(sql)
            out = cursor.fetchall()
            elabConfType = np.dtype({'names': ('version', 'algorithmName', 'smp', 'rsmp', 'window', 'overlap', 'fqMin', 'fqMax', 'filterType', 'filterOrder', 'xcorrMaxLag', 'dtmFilename', 'dtFileame', 'corrM', 'station_id'),
                                    'formats': ('U32', 'U32', 'uint8', 'uint8', 'f8', 'f8', 'f8', 'f8', 'U16', 'uint8', 'f8', 'U32', 'U32', 'uint8', 'U3')})
            nSensors = len(out)
            elabConf = np.zeros(nSensors, dtype=elabConfType)
            elabConf = out
    finally:
        conn.close()

    return elabConf

def isrpInsertElabDts(station_id):

    pymysql.converters.encoders[np.float64] = pymysql.converters.escape_float
    pymysql.converters.conversions = pymysql.converters.encoders.copy()
    pymysql.converters.conversions.update(pymysql.converters.decoders)
    # Connect to the database
    conn = pymysql.connect(host='85.10.202.61',
                           user='geco',
                           password='geco-company-2018',
                           db='GeCoMonit',
                           cursorclass=pymysql.cursors.DictCursor)
    try:

        with conn.cursor() as cursor:
            #TODO dts parameter list & right formatting
            dtsType = np.dtype({'names': ('time', 'ampAv', 'bkz', 'bkzSd', 'aVel', 'aVelSd',
                                          'cohAv', 'fqPk', 'fqMn', 'cnsAv', 'acSns', 'snr',
                                          'latSou', 'lonSou', 'eleSou', 'pdSou', 'station_id'),
                                'formats': ('float', 'float', 'float', 'float', 'float', 'float',
                                          'float', 'float', 'float', 'float', 'float', 'float',
                                            'float', 'float', 'float', 'float', 'str')})
            data = np.zeros(3, dtype=dtsType)
            ndata = data.shape[0]
            stzId = np.repeat(station_id, ndata)
            data['station_id'] = stzId
            sql = "INSERT INTO recordsDts VALUES ('%s', '%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')"
            number_of_rows = cursor.executemany(sql,map(tuple,data))
            conn.commit()

    finally:
        conn.close()

    return str(number_of_rows) + " records written on MySQL dB"







def isrpS2D(elabPrms, sensors, xdem, ydem, zdem, loadT, time, sgn, dtsType):

    dMap = loadT['dMap']
    sMap = loadT['sMap']
    minCorr = np.asscalar(loadT['minCorr'])
    maxCorr = np.asscalar(loadT['maxCorr'])

    # Elab Parameters
    shift = int(elabPrms[0]['overlap'] * elabPrms[0]['smp'] * elabPrms[0]['rsmp'])
    wnd = elabPrms[0]['window'] * elabPrms[0]['smp'] * elabPrms[0]['rsmp']

    w = np.zeros([8, 8, len(xdem), len(ydem)])
    wc = np.zeros([8, len(xdem), len(ydem)])

    corrM = elabPrms[0]['corrM']
    nI = np.sum(np.triu(corrM, 1), axis=(0, 1))
    nD = len(np.where(np.sum(corrM, axis=(1)))[0])

    defCorr = 20 * elabPrms[0]['rsmp']
    l = np.arange(minCorr, maxCorr)
    minCorrIdx = int(np.where(l > -defCorr)[0][0])
    maxCorrIdx = int(np.where(l < defCorr)[0][-1])

    # filtering & resamples
    sgnFlt = butter_bandpass_filter(sgn, elabPrms[0]['fqMin'], elabPrms[0]['fqMin'], elabPrms[0]['smp'], order=elabPrms[0]['filterOrder'])
    time = np.linspace(0, sgnFlt.shape[1] / (elabPrms[0]['smp'] * elabPrms[0]['rsmp']), sgnFlt.shape[1], endpoint=False)
    sgnFlt, time = signal.resample(sgnFlt, sgnFlt.shape[1] * elabPrms[0]['rsmp'], time, axis=1)
    time = np.linspace(0, sgnFlt.shape[1] / (elabPrms[0]['smp'] * elabPrms[0]['rsmp']), sgnFlt.shape[1], endpoint=False)

    # REAL_TIME SYMULATION LOOP
    print('REAL_TIME SYMULATION')
    nSensor = corrM.shape[0]
    tSpan = 60
    tMap = np.zeros([tSpan, nSensor, len(xdem), len(ydem)], 'float16')
    tSpan2 = int(tSpan / 2)

    # CutDist = 4500
    # iR = np.where((ydem > (sensors['Y'][3] - CutDist)) & (ydem < (sensors['Y'][7] + CutDist)))
    # iC = np.where((xdem > (sensors['X'][3] - CutDist)) & (xdem < (sensors['X'][7] + CutDist)))
    # xdem2R, ydem2R = np.meshgrid(xdem[iC], ydem[iR])
    # iCi, iRi = np.meshgrid(iR, iC)
    # zdemR = zdem[iRi, iCi]

    # with writer.saving(fig, "goms-20170305-050000-cav-reckigenGazEx.mp4", 100):
    for t in range(0, len(time) - shift, shift):
        tIns = int(t / shift) % tSpan
        for i in range(0, 8):
            k = 0

            for j in range(i + 1, 8):
                if (corrM[i, j] > 0):
                    k += 1
                    w[i, j, :, :], m = isrpCorrelogram(sgnFlt[:, t:t + wnd], i, j, zdem, dMap, minCorrIdx, maxCorrIdx,
                                                       defCorr)

            tMap[tIns, i, :, :] = np.sum(w, axis=(1))[i, :, :] / k
            tRead = tIns - tSpan2

        # wc[0,:,:]=tMap[tRead,0,:,:]
        for j in range(0, 8):
            wc[j, :, :] = isrpCorrelogram2(tMap, tIns, j, zdem, sMap, shift)

        # print(" tIns "+str(tIns)+" tRead "+str(tRead))
        #
        # print("tr "+str(tMap[tRead,0,100,100]))
        # print("ti "+str(tMap[tIns, 0, 100, 100]))
        wwS = np.nansum(wc, axis=(0)) / nD
        print(int(np.nanmax(wwS) * 100))
        # wwe = wwS / np.nansum(wwS, axis=(0, 1))
        # print(int(np.nanmax(wwe) * 100))

        # Ritaglioì
        wwc = wwS
        print(np.max(wwc))
        wweR = wwc / np.nansum(wwc, axis=(0, 1))

        wwgood = wweR
        wwgood[np.where(wwc < 0.7)] = np.nan

        if np.any(np.isfinite(wwgood)):
            print("LOCATION")
            # 'time', 'ampAv', 'bkz', 'bkzSd', 'aVel', 'aVelSd',
            # 'cohAv', 'fqPk', 'fqMn', 'cnsAv', 'acSns', 'snr',
            # 'latSou', 'lonSou', 'eleSou', 'pdSou', 'station_id')

            # for coll in cMap.collections:
            #     coll.remove()
            we = wwgood * zdem.size / 2
            # cMap = ax2.contourf(xdem2R, ydem2R, we.T * 100, cmap='jet', levels=list(range(int(cThr * 100), 100, 5)),
            #                     alpha=0.9)
            pdmax = np.nanmax([0, np.nanmax(wwgood)])
            cmax = np.nanmax([0, np.nanmax(we)])
            wenodo = wwgood
            uu = np.where(we >= (np.nanmax(we) * .9))
            xx = np.nanmedian(xdem[uu[0]])
            yy = np.nanmedian(ydem[uu[1]])
            zz = np.nanmedian(zdem[uu])

            xsou = np.append(xsou, xx)
            ysou = np.append(ysou, yy)
            zsou = np.append(zsou, zz)
            csou = np.append(csou, cmax)
            tsou = np.append(tsou, t / (elabPrms[0]['smp'] * elabPrms[0]['rsmp']))

            dts = np.zeros(1,dtype=dtsType)
            dts[0][0] = tsou
            dts[0][6] = csou
            dts[0][12] = xsou
            dts[0][13] = ysou
            dts[0][14] = zsou



        #plt.pause(0.001)

            # writer.grab_frame()
    return dts
    #scipy.io.savemat(datafile, {'xsou': xsou, 'ysou': ysou, 'zsou': zsou, 'csou': csou, 'tsou': tsou})
    # fig2 = plt.figure(num=2)
    # ax2 = Axes3D(fig2)
    # ax2.scatter(xsou, ysou, zsou, color='r')
    #


