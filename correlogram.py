
from isrp import *



def gcc_phat(sig, refsig, fs=1, max_tau=None, interp=16):
    '''
    This function computes the offset between the signal sig and the reference signal refsig
    using the Generalized Cross Correlation - Phase Transform (GCC-PHAT)method.
    '''

    # make sure the length for the FFT is larger or equal than len(sig) + len(refsig)
    n = sig.shape[0] + refsig.shape[0]

    # Generalized Cross Correlation Phase Transform
    SIG = np.fft.rfft(sig, n=n)
    REFSIG = np.fft.rfft(refsig, n=n)
    R = SIG * np.conj(REFSIG)
    l = np.arange(0, n)
    cc = np.fft.irfft(R / np.abs(R), n=n)#(interp * n))


   # max_shift = int(interp * n / 2)
   # if max_tau:
   # #      max_shift = np.minimum(int(interp * fs * max_tau), max_shift)
   # #
   # #  cc = np.concatenate((cc[-max_shift:], cc[:max_shift + 1]))
   # #
   # #  # find max cross correlation index
   # #  shift = np.argmax(np.abs(cc)) - max_shift
   #
   #  tau = shift / float(interp * fs)

    return cc , l

def xcorr(x, y, scale='coeff'):
    # Pad shorter array if signals are different lengths
    if x.size > y.size:
        pad_amount = x.size - y.size
        y = np.append(y, np.repeat(0, pad_amount))
    elif y.size > x.size:
        pad_amount = y.size - x.size
        x = np.append(x, np.repeat(0, pad_amount))

    corr = np.correlate(x, y, mode='full')  # scale = 'none'
    lags = np.arange(-(x.size - 1), x.size)

    if scale == 'biased':
        corr = corr / x.size
    elif scale == 'unbiased':
        corr /= (x.size - abs(lags))
    elif scale == 'coeff':
        corr /= np.sqrt(np.dot(x, x) * np.dot(y, y))

    return corr, lags

def isrpCorrelogram(sgn, i, j, zdem, dMap,minCorrIdx,maxCorrIdx,defCorr,norm=False):
    w=np.zeros(zdem.T.shape)
    c, l=xcorr(sgn[i,:],sgn[j,:],'coeff')
    #c,l =gcc_phat(sgn[i,:],sgn[j,:])
    l=int(len(l)/2)
    for k in range(minCorrIdx,maxCorrIdx):
        ll=l-(minCorrIdx-k)-defCorr
        w[dMap[i][j][k][0],dMap[i][j][k][1]]=c[ll]

    #TODO normalizzation
    w=(w+1)/2
    maxW = 0
    if(norm):
        maxW=np.nanmax(w)
        tot=np.nansum(w,axis=(0,1))
        w=w/tot
    return w, maxW

def isrpCorrelogram2(tMap,tRead, i, zdem, sMap,sft):
     w=np.zeros(zdem.T.shape)

     _sft=sft
     maxT=tMap.shape[0]
     for k in range(0,maxT):
         _tRead=tRead-k
         _k=maxT-k
         w[sMap[i][_k][0],sMap[i][_k][1]]=tMap[_tRead,i,sMap[i][_k][0],sMap[i][_k][1]]

     #w=np.nansum(w,axis=(0))
     #TODO normalizzation

     #tot=np.nansum(w,axis=(0,1))
     #w=zdem.size*w/tot
     return w


