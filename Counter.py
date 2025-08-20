import nptdms as td
from scipy.signal import find_peaks as sp_find_peaks
import numpy as np
import json

channels = ['Neut_Mon_1', 'Neut_Mon_2', 'Neut_Mon_3']
paramsFile = open('Params.json', encoding='ascii')
ERROR = 0
try:
    params = json.load(paramsFile)
except UnicodeDecodeError:
    ERROR = 1
    print('КИРИЛЛИЦА НА ПУТИ')
paramsFile.close()


class Signaler:
    def __init__(self, file, filetype, channel, inverted):
        self.inverted = inverted
        self.tdms_file = td.TdmsFile.read(file)
        self.all_groups = self.tdms_file.groups()
        self.group = self.tdms_file[params['counterData']['types'][filetype]]
        self.channel = self.group[channel]
        self.offset = self.channel.properties['Offset']
        if filetype == 'MIV':
            self.channel_data = -self.channel.read_data(scaled=True)
        else:
            self.channel_data = self.channel.read_data(scaled=True)
        self.Gain = self.channel.properties['GAIN']

    @property
    def true_data(self):
        val = (self.channel_data * self.Gain + self.offset) * self.inverted
        return val

    def peak(self):
        return round(np.max(self.true_data), 1)

class Calculator:
    def __init__(self, file, noise_zero, noise_peak, channel, inverted, diff=False):
        self.diff = diff
        self.inverted = inverted
        self.noise_level = noise_zero
        self.threshold_noise = noise_peak
        self.tdms_file = td.TdmsFile.read(file)
        all_channels = []
        for group in self.tdms_file.groups():
            for chan in group.channels():
                all_channels.append(chan.name)
                if chan.name == channel:
                    self.group = group
        if channel not in all_channels:
            return
        self.channel = self.group[channel]
        #self.gain = self.channel.properties['GAIN']
        #self.offset = self.channel.properties['Offset']
        self.rate = self.channel.properties['RATE']
        self.time = np.zeros(len(self.true_data()))
        for i in range(0, len(self.true_data())):
            self.time[i] = i/self.rate

    def true_data(self):
        self.base = 0
        channel_data = self.channel.read_data()
        val = channel_data * self.inverted
        #val = (channel_data * self.gain + self.offset) * self.inverted
        print(4.07)
        for i in range(0, 10):
            self.base += val[i]
        self.base = self.base/10
        val = val - self.base
        return val

    def integrate(self, time_constant):
        data = self.true_data()
        val = np.zeros(len(data))
        signal = 0
        for i in range(1, len(data)):
            if data[i] > self.noise_level:
                cur = data[i]
            else:
                cur = 0
            signal = signal / (2.7182818 ** ((self.time[i]-self.time[i-1])/time_constant)) + cur
            val[i] = signal
        return val

    def differentiate(self):

        data = self.true_data()
        val = np.zeros(len(data))
        for i in range(1, len(data)):
            if data[i] - data[i-1] > 0:
                val[i] = data[i] - data[i-1]
        return val

    def find_peaks(self):
        if self.diff:
            data = self.differentiate()
        else:
            data = self.true_data()
        self.peaks = sp_find_peaks(data, height=0, threshold=self.threshold_noise)[0]
        self.data_peaks = data[self.peaks]
        peaks_2 = np.argwhere(data > self.noise_level)
        self.peaks = np.intersect1d(self.peaks, peaks_2)
        self.data_peaks = self.data_peaks[self.data_peaks > self.noise_level]
        self.amount = len(self.data_peaks)
        return self.peaks, self.data_peaks, self.amount

class MultiCalcultor():
    def __init__(self, files, noise_zero, noise_peak):
        self.Calc = []
        for file in files:
            for channel in channels:
                temp = Calculator(file, noise_zero, noise_peak, channel)
                self.Calc.append(temp)

    def find_peaks(self):
        self.locamount = []
        for calc in self.Calc:
            _, _, locamount = calc.find_peaks()
            self.locamount.append(locamount)
        return self.locamount


"""
DIAGNOSTICS

print(f'true_max = {val_max}, true_min = {val_min}')
print(f'max*gain = {max(channel_data)*Gain}, min*gain = {min(channel_data)*Gain}')
print(f'Sdvigs = {Sdvig}')
print(f'Sdvig/Offset = {Sdvig/Offset}, Sdvig*Offset = {Sdvig*Offset}, Sdvig/Gain = {Sdvig/Gain}, Sdvig*Gain = {Sdvig*Gain}')
print(f'1 - 9 {channel_data[0:9]}')
print(f'third = {true_data[2]}')
print(f'min = {min(true_data)}')
print(f'Мои {(val_max-val_min)}, их {max(channel_data)*Gain-min(channel_data)*Gain}')
"""
