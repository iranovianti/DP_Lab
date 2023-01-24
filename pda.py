import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks

class PDA:
	def __init__(self, file_path, name=''):
		self.name = name

		self.data = self.get_file(file_path)
		self.t_step = self.data.index.to_numpy()
		self.wl_step = self.data.columns.to_numpy()

		self.t_unit = 'min'

		self.peaks = {}

		plt.rc('font', size=15)
		plt.rc('legend', fontsize=14)

	def get_file(self, file_path):
		#function to read output PDA file from HPLC JASCO

		#1. Reading the header to get measurement information
		header = open(file_path, 'r').read().split('\n')[1:9]
		header = [text.split('\t') for text in header]
		units = {s[0]: get_value_unit(s[1]) for s in header}

		#2. Get time information
		t_start = units['START_TIME']['value']
		t_points = int(units['TIME_NPOINTS']['value'])
		t_int = units['TIME_INTERVAL']['value'] / 60000 #time interval in minutes (min = 60000 ms)

		#3. Get wavelength information
		wl_start = units['START_WL']['value']
		wl_points = int(units['WL_NPOINTS']['value'])
		wl_int = units['WL_INTERVAL']['value']

		#4. Read intensity values
		pda_df = pd.read_csv(file_path, sep='\t', skiprows=12, header=None)
		pda_df = pda_df.dropna(axis='columns')

		#5. Assign wavelength and time to each measured point
		pda_df.index = steps(t_start, t_points, t_int)
		pda_df.columns = steps(wl_start, wl_points, wl_int)
		data = pda_df

		return data

	def info(self):
		print(f'PDA name: {self.name}')
		print(f'Wavelength range: {self.wl_step[0]} - {self.wl_step[-1]} nm ({len(self.wl_step)} points)')
		print(f'Monitored time: {self.t_step[-1]:.4} {self.t_unit}')

	def _plot_func(self, xs, ys, xlabel, ylabel,
		legend=False, **attr):
		"""
		kwargs:
		xlim: limit of x axis, default -> min and max value of xs
		ylim: limit of y axis, default -> not set
		title: title
		colors: list of colors to use for each data, must be the same length as the data
		labels: list of labels
		"""
		assert len(xs) == len(ys)

		xlim = attr.get('xlim', (np.min(xs), np.max(xs)))
		ylim = attr.get('ylim')
		title = attr.get('title')
		colors = attr.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color'])
		labels = attr.get('labels', np.full(len(xs), ' '))
		figsize= attr.get('figsize', (6,6))

		plt.figure(figsize=figsize)
		for i in range(len(xs)):
			plt.plot(xs[i], ys[i], color=colors[i], label=labels[i], lw=1.5)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		if labels and legend:
			plt.legend(frameon=False)
		plt.xlim(xlim)
		if ylim:
			plt.ylim(ylim)
		if title:
			plt.title(title)
		plt.gca().spines['right'].set_color('none')
		plt.gca().spines['top'].set_color('none')
		plt.show()

	def plot_at_rt(self, *rt, **attr):
		rts = list(rt)
		adj_rts = adjust_values(rts, self.t_step)

		xs = []
		ys = []
		labels = []

		for r_t in adj_rts:
			xs.append(self.wl_step)
			ys.append(self.data.loc[r_t])
			labels.append(f'at {float(r_t):.4} {self.t_unit}')

		xlabel = 'Wavelength (nm)'
		ylabel = 'Intensity'

		if len(adj_rts) == 1:
			self._plot_func(xs, ys,
							xlabel=xlabel, ylabel=ylabel,
							title=f'{self.name}\nspectra at {labels[0]}', **attr)
		else:
			self._plot_func(xs, ys, labels=labels,
							xlabel=xlabel, ylabel=ylabel,
							title=self.name, legend=True, **attr)

	def plot_at_wl(self, *wl, **attr):
		wls = list(wl)
		adj_wls = adjust_values(wls, self.wl_step)

		xs = []
		ys = []
		labels = []

		for w_l in adj_wls:
			xs.append(self.t_step)
			ys.append(self.data.loc[:][w_l])
			labels.append(f'{int(w_l)} nm')

		xlabel = f'retention time ({self.t_unit})'
		ylabel = 'Intensity'

		figsize = attr.get('figsize', (10,4))

		if len(adj_wls) == 1:
			self._plot_func(xs, ys, figsize=figsize,
							xlabel=xlabel, ylabel=ylabel,
							title=f'{self.name}\nmonitored at {labels[0]}', **attr)
		else:
			self._plot_func(xs, ys, labels=labels, figsize=figsize,
							xlabel=xlabel, ylabel=ylabel,
							title=self.name, legend=True, **attr)

	def print_peak(self):
		print(f'{len(self.peaks)} peaks registered')
		if len(self.peaks) > 0:
			for k,v in self.peaks.items():
				print(f'{k} at {v:.4} {self.t_unit}')

	def show_peak(self, wl=200):
		if len(self.peaks) > 0:
			data = np.array(self.data.loc[:][wl])
			plt.figure(figsize=(10,4))
			plt.plot(self.t_step, data)
			for k,v in self.peaks.items():
				x = v
				y = self.data.loc[v][wl]
				plt.plot(x, y, 'o', c='red', ms=10)
				plt.annotate(k, (x,y), size=15)
			plt.xlim(np.min(self.t_step), np.max(self.t_step))
			plt.xlabel(f'retention time ({self.t_unit})')
			plt.ylabel('Intensity')
			plt.show()

	def detect_peak(self, wl=200, show=True):
		wl_list = np.array(self.wl_step)
		if wl in wl_list:
			wl_ref = wl
		else:
			print(f'{wl} is not listed')
			adj_wl = find_nearest(wl_list, wl)
			wl_ref = adj_wl
			
		data = np.array(self.data.loc[:][wl_ref])

		print(f'detecting peaks from chromatogram at {wl_ref} nm')
		
		peak_prom = np.std(data) * 2
		peak_ids, _ = find_peaks(data, prominence=(peak_prom, None))

		if len(peak_ids) == 0:
			print('no peaks detected')

		else:
			peaks = {f'peak_{i+1}': self.t_step[peak] for i,peak in enumerate(peak_ids)}
			self.peaks = peaks
			self.print_peak()
			if show:
				self.show_peak(wl=wl_ref)

	def plot_peaks(self, *peaks, **attr):
		if len(self.peaks) > 0:
			rts = [self.peaks[f'peak_{p}'] for p in list(peaks)]
			self.plot_at_rt(*rts, **attr)

#Helper functions

def get_value_unit(string):
	unit = string.strip('0123456789.')
	value = float(string.replace(unit, ''))
	return {'value': value, 'unit': unit}

def steps(start, npoints, interval):
	return np.arange(start, start+(npoints*interval), interval)

def find_nearest(a, num):
	return a[np.abs(a - num).argmin()]

def adjust_values(val_list, a):
	adjusted = []
	a = np.array(a)
	for val in val_list:
		if val in a:
			adjusted.append(val)
		else:
			adj_val = find_nearest(a, val)
			adjusted.append(adj_val)
			if val != round(adj_val, 2):
				print(f'{val} is not listed, adjusted to {adj_val}')
	return adjusted
