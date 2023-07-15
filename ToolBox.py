from scipy import signal
from scipy.fft import fftshift
import numpy as np
import scipy.io as sio
import copy
from matplotlib import pyplot as plt
from lmfit.models import GaussianModel, ConstantModel
from pathlib import Path

import mat73
from pathlib import Path

from .DPI import hl2adb as hldb
from .DPI import changeDriver


def hl_specgram(shotNumber, channelName=['CPS07', 'CPS08'], ts=100, te=200, fs=2e6,sys='CPS', *args):
    t_I, y0, U = hldb(shotNumber, channelName[0], ts/1e3, te/1e3, fs, sys)
    t_Q, y1, U = hldb(shotNumber, channelName[1], ts/1e3, te/1e3, fs, sys)
    A_I=signal.detrend(y0)
    A_Q=signal.detrend(y1)
    signalWave = A_I+A_Q*1J

    f,t,s=signal.spectrogram(signalWave,fs,'hanning',noverlap=512,nfft=1024,detrend='linear',mode='psd')
    # fig, ax = plt.subplots()
    plt.pcolormesh(t*1e3+t_I[0],fftshift(f),10*np.log10(fftshift(s, axes=0)),cmap='jet')
    plt.colorbar()
    return t*1e3+t_I[0], fftshift(f), U   # s->ms

def smooth(x,num=5):
    if num // 2 == 0: # 偶数转奇数
        num -= 1
    length = len(x)
    y = np.zeros(length)
    N = (num - 1) / 2
    for i in range(0, length):
        cont_0 = i
        cont_end = length - i - 1
        if cont_0 in range(0,int(N)) or cont_end in range(0,int(N)):
            cont = min(cont_0,cont_end)
            y[i] = np.mean(x[i - cont : i + cont + 1])
        else:
            y[i] = np.mean(x[i - int(N) : i + int(N) + 1])
    return y

def decorate(**options):
    """Decorate the current axes.

    Call decorate with keyword arguments like

    decorate(title='Title',
             xlabel='x',
             ylabel='y')

    The keyword arguments can be any of the axis properties

    https://matplotlib.org/api/axes_api.html

    In addition, you can use `legend=False` to suppress the legend.

    And you can use `loc` to indicate the location of the legend
    (the default value is 'best')
    """
    loc = options.pop("loc", "best")
    if options.pop("legend", True):
        legend(loc=loc)

    plt.gca().set(**options)
    plt.tight_layout()


def legend(**options):
    """Draws a legend only if there is at least one labeled item.

    options are passed to plt.legend()
    https://matplotlib.org/api/_as_gen/matplotlib.plt.legend.html

    """
    underride(options, loc="best", frameon=False)

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, **options)

def underride(d, **options):
    """Add key-value pairs to d only if key is not in d.

    If d is None, create a new dictionary.

    d: dictionary
    options: keyword args to add to d
    """
    if d is None:
        d = {}

    for key, val in options.items():
        d.setdefault(key, val)

    return d

class Signal():
    
    pass
    
class complexSignal():

    """Represents a discrete-time waveform.

    """

    def __init__(self, shotNum, channelList, ts_s=100, ts_e=1500, step=1, Device = 'HL-2M'):
        """Initializes the wave.segment = wave4.segment(start=1.15, duration=0.85)

        ys: wave array
        ts: array of times
        framerate: samples per second
        """
        # 0=driver mode   0=hl2a,1=local, 2=exl50,3=east,4=hl2m
        if Device == 'HL-2A':
            changeDriver(0)
        elif Device == 'HL-2M':
            changeDriver(4)
        else:
            changeDriver(1)
        self.Device = Device
        self.shotNum = shotNum
        self.chl_I = channelList[0]
        self.chl_Q = channelList[1]
        ts1, y1, U = hldb(shotNum,channelList[0], ts_s/1e3, ts_e/1e3,step)
        ts2, y2, U = hldb(shotNum,channelList[1], ts_s/1e3, ts_e/1e3,step)

        self.ys = np.asanyarray(y1+y2*1J)
        self.framerate_k =  round(1/(ts1[1]-ts1[0]))
        self.framerate = round(1e3/(ts1[1]-ts1[0]))
        self.ts = np.asanyarray(ts1)
        self.fftPara = None
        self.spect_Sxx = None
        self.fd = None

    def copy(self):
        """Makes a copy.

        Returns: new Wave
        """
        return copy.deepcopy(self)

    def __len__(self):
        return len(self.ys)

    @property
    def start(self):
        return self.ts[0]

    @property
    def end(self):
        return self.ts[-1]

    @property
    def duration(self):
        """Duration (property).

        returns: float duration in microseconds
        """
        return len(self.ys) / self.framerate_k
    
    def find_index(self, t):
        """Find the index corresponding to a given time."""
        n = len(self)
        start = self.start
        end = self.end
        i = round((n - 1) * (t - start) / (end - start))
        return int(i)

    def segment(self, start=None, end=None):
        """Extracts a segment.

        start: float start time in microseconds
        end: float duration in microseconds

        returns: Wave
        """
        if start is None:
            start = self.ts[0]
            i = 0
        else:
            i = self.find_index(start)

        j = None if end is None else self.find_index(end)
        return self.slice(i, j)

    def slice(self, i, j):
        """Makes a slice from a Wave.

        i: first slice index
        j: second slice index
        """
        Signal_slice = self.copy()
        Signal_slice.ys = self.ys[i:j].copy()
        Signal_slice.ts = self.ts[i:j].copy()

        if self.spect_Sxx is not None:
            ts0 = self.spect_t[0]
            te0 = self.spect_t[-1]
            start = Signal_slice.ts[0]
            end = Signal_slice.ts[-1]
            n = len(self.spect_t)
            i = round((n - 1) * (start - ts0) / (te0 - ts0))
            j = round((n - 1) * (end - ts0) / (te0 - ts0))
            Signal_slice.spect_Sxx = self.spect_Sxx[...,i:j].copy()
            Signal_slice.spect_t = self.spect_t[i:j].copy()

            if self.fd is not None:
                Signal_slice.fd = self.fd[i:j].copy()
            else:
                Signal_slice.fd = None

        else:
            Signal_slice.spect_Sxx = None
        return Signal_slice
    def get_xfactor(self, options):
        try:
            xfactor = options["xfactor"]
            options.pop("xfactor")
        except KeyError:
            xfactor = 1
        return xfactor
    
    def plot(self, **options):
        """Plots the wave.

        If the ys are complex, plots the real part.

        """
        xfactor = self.get_xfactor(options)
        # TODO:
        # if isinstance(self.ys,complex):
        plt.plot(self.ts * xfactor, np.real(self.ys), label=self.chl_I, **options)
        plt.plot(self.ts * xfactor, np.imag(self.ys), label=self.chl_Q, **options)

        decorate(title=self.Device+':#'+str(self.shotNum),xlabel='Time (ms)',ylabel='Amplitude',legend=True)
        # else:
            # plt.plot(self.ts * xfactor, np.real(self.ys), **options)
    
    def spectrogram(self,para=[512*2,350*2,1024*2],isPlot=False,isDecorate=True,ylim=(-500,500),window='hamming',cmap='jet',mode='psd',scaling='density',**options):
        """
        scaling : { 'density', 'spectrum' }, optional
            Selects between computing the power spectral density ('density')
            where `Sxx` has units of V**2/Hz and computing the power
            spectrum ('spectrum') where `Sxx` has units of V**2, if `x`
            is measured in V and `fs` is measured in Hz. Defaults to
            'density'.
        mode : str, optional
            Defines what kind of return values are expected. Options are
            ['psd', 'complex', 'magnitude', 'angle', 'phase']. 'complex' is
            equivalent to the output of `stft` with no padding or boundary
            extension. 'magnitude' returns the absolute magnitude of the
            STFT. 'angle' and 'phase' return the complex angle of the STFT,
            with and without unwrapping, respectively.
        """
        if self.fftPara != para:
            self.fftPara = para
            f,t,Sxx=signal.spectrogram(self.ys,fs=self.framerate,window=window,nperseg=para[0],noverlap=para[1],nfft=para[2],detrend='linear',mode=mode,scaling=scaling,return_onesided=False)
            self.spect_t = t*1e3 + self.start
            self.spect_f = fftshift(f/1e3)
            self.spect_Sxx = fftshift(Sxx, axes=0)
        else:
            pass
        if isPlot is True:
            im = plt.pcolormesh(self.spect_t, self.spect_f, 10*np.log10(self.spect_Sxx), cmap=cmap, shading='auto',**options)
            plt.ylim(ylim)
            return im
            if isDecorate is True:
                decorate(title=self.Device+':#'+str(self.shotNum),xlabel='Time (ms)',ylabel='Freq (kHz)',legend=True)
                plt.colorbar()
        # signal.periodogram()
        # plt.specgram(self.ys,NFFT=para[0],noverlap=para[1],Fs = self.framerate,detrend='linear',cmap=cmap, scale_by_freq=True,mode=mode,sides='default',**options)
        else:
            pass
        

    def plot_power(self, start=None, end=None,xlim=(-800,800),bias=0,isPlotFlip=False,isPlotDelta=False,label=None,smoothSize=None):
        """Extracts a segment.

        start: float start time in microseconds
        end: float duration in microseconds

        returns: Wave
        """
        n = len(self.spect_t)
        ts0 = self.ts[0]
        te0 = self.ts[-1]

        if start is None:
            start = ts0
            i = 0
        else:
            i = round((n - 1) * (start - ts0) / (te0 - ts0))

        if end is None:
            end = te0
            j = None
        else:
            j = round((n - 1) * (end - ts0) / (te0 - ts0))

        
        f = self.spect_f
        Pxx0 = self.spect_Sxx[...,i:j]
        Pxx = np.mean(Pxx0,axis=1)
        self.power_Pxx = Pxx
        self.power_Pxx_log = 10*np.log10(Pxx)
        if label is None:
            label=str(round(start,1))+'-'+str(round(end,1))+' ms'

        if smoothSize is not None:
            p = plt.plot(smooth(f,smoothSize), smooth(10*np.log10(Pxx)+bias,smoothSize),label=label)
        else:
            p = plt.plot(f, 10*np.log10(Pxx)+bias,label=label)
        if isPlotFlip:
            plt.plot(np.flip(f), 10*np.log10(Pxx),label='flip:'+label)
        plt.xlim(xlim)
        if self.fd is not None:
            fd = np.mean(self.fd)
            plt.axvline(x=fd,ls='--',c=p[0].get_color())


        if isPlotDelta:
            p = plt.plot(f,10*np.log10(np.abs(Pxx-np.flip(Pxx))))
            # p = plt.plot(f,10*np.log10(np.mean(np.abs(Pxx0-np.flip(Pxx0,axis=0)),axis=1)))
        # plt.axvline(x=0,ls='--',c='black')
        decorate(title=self.Device+':#'+str(self.shotNum),xlabel='Freq (kHz)',ylabel='Power(V^2)',legend=True)
        # return self.slice(i, j)

    def gaussian3(self,data=None,para=[2e-8,-200,1,6e-8,200,1],isLog=True):
        if data is None:
            x = self.spect_f
            y = self.power_Pxx
        else:
            x = data[0]
            y = data[1]
        gmodel = GaussianModel(prefix='p1_') + GaussianModel(prefix='p2_')+GaussianModel(prefix='p3_')#+GaussianModel(prefix='p4_')
        # make parameters with starting values:
        # params = gmodel.make_params(p1_amplitude=2e-8, p1_center=-200, p1_sigma=0.2,
        #                             p2_amplitude=6e-8, p2_center=10, p2_sigma=0.2)
        params = gmodel.make_params(p1_amplitude=para[0], p1_center=para[1], p1_sigma=para[2],
                                    p2_amplitude=para[3], p2_center=para[4], p2_sigma=para[5],
                                    #p4_amplitude=para[6], p4_center=para[7], p4_sigma=para[8],
                                    p3_amplitude=np.max(y), p3_center=0, p3_sigma=2)
        # it's not really needed for this data, but you can put bounds on
        # # parameters like this (or set .vary=False to fix a parameter)
        # params['peak_sigma'].min = 0         # sigma  > 0
        # params['peak_amplitude'].min = 1e-8     # amplitude < 0
        # params['peak_center'].min =-300
        # params['peak_center'].max = 200
        params['p3_center'].set(vary=False)
        params['p1_amplitude'].min = 0
        params['p2_amplitude'].min = 0
        params['p3_amplitude'].min = 0
        # params['p4_amplitude'].min = 0

        # run fit
        result = gmodel.fit(y, params, x=x)
        self.gausssianFit = result
        # amp1 = result.best_values["p1_amplitude"]
        cen1 = result.best_values["p1_center"]
        # wid1 = result.best_values["p1_sigma"]
        # amp2 = result.best_values["p2_amplitude"]
        cen2 = result.best_values["p2_center"]
        # wid2 = result.best_values["p2_sigma"]
        # amp3 = result.best_values["p3_amplitude"]
        # cen3 = result.best_values["p3_center"]
        # wid3 = result.best_values["p3_sigma"]
        # g1 = (amp1/(np.sqrt(2*np.pi)*wid1)) * np.exp(-np.power(x - cen1, 2.) / (2 * np.power(wid1, 2.)))
        # g2 = (amp2/(np.sqrt(2*np.pi)*wid2)) * np.exp(-np.power(x - cen2, 2.) / (2 * np.power(wid2, 2.)))
        # g3 = (amp3/(np.sqrt(2*np.pi)*wid3)) * np.exp(-np.power(x - cen3, 2.) / (2 * np.power(wid3, 2.)))
        comps = result.eval_components(x=x)
        g1 = comps['p1_']
        g2 = comps['p2_']
        g3 = comps['p3_']

        if isLog:
            plt.plot(x, 10*np.log10(result.best_fit), 'r-',label='Sum')
            l1=plt.plot(x, 10*np.log10(g1),label='peak1:'+str(round(cen1,1))+'KHz')
            l2=plt.plot(x, 10*np.log10(g2),label='peak2:'+str(round(cen2,1))+'KHz')
            l3=plt.plot(x, 10*np.log10(g3),label='Symmetric')
        else:
            plt.plot(x, result.best_fit, 'r-',label='Sum')
            l1=plt.plot(x, g1,label='peak1:'+str(round(cen1,1))+'KHz')
            l2=plt.plot(x, g2,label='peak2:'+str(round(cen2,1))+'KHz')
            l3=plt.plot(x, g3,label='Symmetric')

        plt.legend()

    def periodogram(self,nfft=512,xlim=(-500,500),window='hamming',scaling='spectrum',**options):
        """
        window : str or tuple or array_like, optional
            Desired window to use. If `window` is a string or tuple, it is
            passed to `get_window` to generate the window values, which are
            DFT-even by default. See `get_window` for a list of windows and
            required parameters. If `window` is array_like it will be used
            directly as the window and its length must be equal to the length
            of the axis over which the periodogram is computed. Defaults
            to 'boxcar'.
        nfft : int, optional
            Length of the FFT used. If `None` the length of `x` will be
            used.
        return_onesided : bool, optional
            If `True`, return a one-sided spectrum for real data. If
            `False` return a two-sided spectrum. Defaults to `True`, but for
            complex data, a two-sided spectrum is always returned.
        scaling : { 'density', 'spectrum' }, optional
            Selects between computing the power spectral density ('density')
            where `Pxx` has units of V**2/Hz and computing the power
            spectrum ('spectrum') where `Pxx` has units of V**2, if `x`
            is measured in V and `fs` is measured in Hz. Defaults to
            'density'
        """
        f, Pxx = signal.periodogram(self.ys,fs=self.framerate,window=window,nfft=nfft,detrend='constant',return_onesided=False,scaling=scaling,**options)
        f = fftshift(f/1e3)
        Pxx = fftshift(Pxx, axes=0)

        plt.plot(f, 10*np.log10(Pxx),label=str(round(self.start,1))+'-'+str(round(self.end,1))+' ms')
        plt.xlim(xlim)
        decorate(title=self.Device+':#'+str(self.shotNum),xlabel='Freq (kHz)',ylabel='Power(V^2)',legend=True)
        # plt.xlabel('frequency [Hz]')
        # plt.ylabel('PSD [V**2]')
        # plt.show()
        return f,Pxx
    def extractfd(self,para=None,isPlot=True,smoothSize=5,ignore=[0,1000] ,**options):
        if para is not None and self.fftPara != para:
            self.spectrogram(para)
        elif para is None and self.fftPara is None:
            print('Input Para!')
        # elif para is not None and self.fftPara == para:
        #     pass
        # elif para is None and self.fftPara != para:
        #     pass
        else:
            pass
        # f = self.spect_f[1:] # 使f对称
        # index1 = (f<0)
        # index2 = (f>0)
        # f1 = np.flip(f[index1])
        # f2 = f[index2]
        # Sxx = self.spect_Sxx[1:] # 使f对称
        # Sxx1 = np.flip(Sxx[index1,...])
        # Sxx2 = Sxx[index2,...]
        # self.fd = (f1.dot(Sxx1)+f2.dot(Sxx2))/np.sum(np.abs(Sxx1-Sxx2),axis=0)

        f = self.spect_f[1:] # 使f对称
        # index1 =np.logical_and(f<-1*ignore,f>-1*(self.framerate_k/2-ignore))
        # index2 = np.logical_and(f>ignore,f<self.framerate_k/2-ignore)
        index1 =np.logical_and(f<-1*ignore[0],f>-1*ignore[1])
        index2 = np.logical_and(f>ignore[0],f<ignore[1])
        f1 = np.flip(f[index1])
        f2 = f[index2]
        Sxx = self.spect_Sxx[1:,...] # 使f对称
        Sxx1 = np.flip(Sxx[index1,...],axis=0)
        Sxx2 = Sxx[index2,...]

        # Sxx1 = smooth(Sxx1,smoothSize)
        # Sxx2 = smooth(Sxx2,smoothSize)
        self.fd = (f1.dot(Sxx1)+f2.dot(Sxx2))/np.sum(np.abs(Sxx1-Sxx2),axis=0)
        # fd1=sum(ff1.*pp1+flip(ff2).*flip(pp2),1)./sum(abs(pp1-flip(pp2)),1);
        # fd2=sum(F1/1000.*y0)./sum(y0,1);
        if isPlot is True:
            plt.plot(smooth(self.spect_t,smoothSize),smooth(self.fd,smoothSize),**options)
            # plt.plot(smooth(self.spect_t,smoothSize),smooth(self.fd2,smoothSize),**options)
            # plt.plot(smooth(self.spect_t,smoothSize),smooth(self.fd2,smoothSize),**options)



