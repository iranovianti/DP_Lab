import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks

class PDA:
	def __init__(self, file_path, name=''):
		self.name = name
		header = open(file_path, 'r').read().split('\n')[1:9]
		header = [text.split('\t') for text in header]
		self.units = {s[0]: get_value_unit(s[1]) for s in header}

		t_start = self.units['START_TIME']['value']
		t_points = int(self.units['TIME_NPOINTS']['value'])
		t_int = self.units['TIME_INTERVAL']['value'] / 60000 #time interval in minutes (min = 60000)
		self.t_step = steps(t_start, t_points, t_int)

		wl_start = self.units['START_WL']['value']
		wl_points = int(self.units['WL_NPOINTS']['value'])
		wl_int = self.units['WL_INTERVAL']['value']
		self.wl_step = steps(wl_start, wl_points, wl_int)

		pda_df = pd.read_csv(file_path, sep='\t', skiprows=12, header=None)
		pda_df = pda_df.dropna(axis='columns')
		pda_df.index = self.t_step
		pda_df.columns = self.wl_step
		self.data = pda_df

		self.peaks = {}

		plt.rc('font', size=15)
		plt.rc('legend', fontsize=14)

	def info(self):
		print(f'PDA name: {self.name}')
		for k,v in self.units.items():
			print(f"{k}: {v['value']} {v['unit']}")

	def _plot_func(self, xs, ys, xlabel, ylabel, figsize=(6,6),
					xlim=None, ylim=None, legend=False,
					title='', labels=None, colors=None):
		assert len(xs) == len(ys)
		plt.figure(figsize=figsize)

		if colors:
			cols = colors
		else:
			cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

		if labels:
			lbls = labels
		else:
			lbls = np.full(len(xs), ' ')

		for i in range(len(xs)):
			plt.plot(xs[i], ys[i], color=cols[i], label=lbls[i], lw=1.5)

		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		if labels and legend:
			plt.legend(frameon=False)
		if xlim:
			plt.xlim(xlim)
		if ylim:
			plt.ylim(ylim)
		plt.gca().spines['right'].set_color('none')
		plt.gca().spines['top'].set_color('none')
		plt.title(title)
		plt.show()

	def plot_at_rt(self, *rt):
		rts = list(rt)
		adj_rts = adjust_values(rts, self.t_step)

		xs = []
		ys = []
		labels = []

		for r_t in adj_rts:
			xs.append(self.wl_step)
			ys.append(self.data.loc[r_t])
			labels.append(f'at {float(r_t):.4} min')

		xlabel = 'Wavelength (nm)'
		ylabel = 'Intensity'
		xlim = (np.min(self.wl_step), np.max(self.wl_step))
		
		if len(adj_rts) == 1:
			self._plot_func(xs, ys, xlim=xlim,
							xlabel=xlabel, ylabel=ylabel,
							title=f'{self.name}\nspectra at {labels[0]}')
		else:
			self._plot_func(xs, ys, xlim=xlim, labels=labels,
							xlabel=xlabel, ylabel=ylabel,
							title=self.name, legend=True)

	def plot_at_wl(self, *wl):
		wls = list(wl)
		adj_wls = adjust_values(wls, self.wl_step)

		xs = []
		ys = []
		labels = []

		for w_l in adj_wls:
			xs.append(self.t_step)
			ys.append(self.data.loc[:][w_l])
			labels.append(f'{w_l} nm')

		xlabel = 'retention time (min)'
		ylabel = 'Intensity'
		figsize = (10,4)
		xlim = (np.min(self.t_step), np.max(self.t_step))

		if len(adj_wls) == 1:
			self._plot_func(xs, ys, xlim=xlim, figsize=figsize,
							xlabel=xlabel, ylabel=ylabel,
							title=f'{self.name}\nmonitored at {labels[0]}')
		else:
			self._plot_func(xs, ys, xlim=xlim, labels=labels, figsize=figsize,
							xlabel=xlabel, ylabel=ylabel,
							title=self.name, legend=True)

	def print_peak(self):
		print(f'{len(self.peaks)} peaks registered')
		if len(self.peaks) > 0:
			for k,v in self.peaks.items():
				print(f'{k} at {v:.4} min')

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
			plt.xlabel('retention time (min)')
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
				self.show_peak()

	def plot_peaks(self, *peaks):
		if len(self.peaks) > 0:
			rts = [self.peaks[f'peak_{p}'] for p in list(peaks)]
			self.plot_at_rt(*rts)

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