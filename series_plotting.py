import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob

from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.lines as mlines

def make_colormap(seq):
	"""source: https://stackoverflow.com/a/16836182
	Return a LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	"""
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = seq[i - 1]
			r2, g2, b2 = seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)

class SpectraSeries:

	def __init__(self, mode='Absorbance', series_measure='time', series_unit='mins', series_start=0, series_interval=1,
		font_size=18, legend_font_size=15):
		self.data_series = []
		self.mode = mode
		self.x_label = 'Wavelength (nm)'
		self.series_measure = series_measure
		self.series_unit = series_unit
		self.series_start = series_start
		self.series_interval = series_interval

		self.font_size = font_size
		self.legend_font_size = legend_font_size

		plt.rc('font', size=self.font_size)
		plt.rc('legend', fontsize=self.legend_font_size)

	def open_abs(self, file_path):
		_abs = pd.read_csv(file_path, sep='\t', skiprows=1)
		x = _abs['Wavelength nm.']
		y = _abs['Abs.']
		return x,y

	def add_file(self, file_path, file_type, preview=False):
		assert file_type in ['single', 'singles', 'multi']

		print('adding...')

		if file_type == 'single':
			f = self.open_abs(file_path)
			self.data_series.append(f)
		
		elif file_type == 'singles':
			filenames = sorted(glob(file_path+'/*.txt'))
			abs_list = [self.open_abs(fn) for fn in filenames]
			for f in abs_list:
				self.data_series.append(f)
		
		elif file_type == 'multi':
			file_stack = pd.read_csv(file_path, sep='\t', header=None, skiprows=1)
			abs_list = [(file_stack[0], file_stack[i]) for i in file_stack.columns[1:]]
			for f in abs_list:
				self.data_series.append(f)
		
		print(f'total number of data: {len(self.data_series)}')
		
		if preview:
			self.spectra(title='Preview')
  
	def remove_data(self, pos):
		self.data_series.pop(pos-1)
		print(f'total number of data: {len(self.data_series)}')
  
	def step(self):
		return np.arange(self.series_start,
			self.series_start+(len(self.data_series)*self.series_interval),
			self.series_interval)
  
	def spectra(self, figure_size=(6,6), color='black', xlim='auto', ylim='auto', title='title'):
		plt.figure(figsize=figure_size)
		for s in self.data_series:
			x,y = s
			plt.plot(x, y, color=color)
		
		if xlim == 'auto':
			plt.xlim(np.min(self.data_series[0][0]), np.max(self.data_series[0][0]))
		else:
			plt.xlim(xlim)

		if ylim == 'auto':
			y_max = np.max([a[1] for a in self.data_series])
			y_min = np.min([a[1] for a in self.data_series])
			gap = 0.15*(y_max-y_min)
			plt.ylim(y_min-(gap/2), y_max+gap)
		else:
			plt.ylim(ylim)
		
		plt.gca().spines['right'].set_color('none')
		plt.gca().spines['top'].set_color('none')
		plt.xlabel(self.x_label)
		plt.ylabel(self.mode)
		plt.title(title)
		plt.show()
  
	def check_data(self, pos, figure_size=(6,6)):
		assert pos in range(1, len(self.data_series)+1)
		
		print(f'Showing data at position: {pos}')
		
		plt.figure(figsize=figure_size)
		x,y = self.data_series[pos-1]
		plt.plot(x,y)
		plt.show()

	def spectra_gradient(self, figure_size=(6,6), color1='blue', color2='red', xlim='auto', 
		title='title', comment='', ylim='auto', cbar=True):

		data = [np.column_stack(s) for s in self.data_series]

		ccv = mcolors.ColorConverter().to_rgb
		colormap = make_colormap([ccv(color1), ccv(color2)])

		fig, ax = plt.subplots(figsize=figure_size)

		line_segments = LineCollection(data, linestyles='solid', cmap=colormap)
		line_segments.set_array(self.step())
		ax.add_collection(line_segments)

		if cbar:
			cbaxes = inset_axes(ax, width="70%", height="3%", loc=1) 
			axcb = fig.colorbar(line_segments, cax=cbaxes, orientation='horizontal')
			axcb.ax.tick_params(labelsize=15)
			axcb.set_label(f'{self.series_unit}\n{comment}', fontsize=15)

		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.mode)
		ax.set_title(title)

		if xlim == 'auto':
			ax.set_xlim(np.min(self.data_series[0][0]), np.max(self.data_series[0][0]))
		else:
			ax.set_xlim(xlim)

		if ylim == 'auto':
			y_max = np.max([a[1] for a in self.data_series])
			y_min = np.min([a[1] for a in self.data_series])
			gap = 0.15*(y_max-y_min)
			ax.set_ylim(y_min-(gap/2), y_max+gap)
		else:
			ax.set_ylim(ylim)
		
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		plt.show()
	
	def spectra_disc(self, labels=None, color_palette='viridis', figure_size=(6,6), legend=False, title='title', custom_colors=None, xlim='auto',
			 ylim=None, legend_title='', savefig=None):
		if labels:
			assert len(labels) == len(self.data_series)
			label_list = labels
		else:
			label_list = [f'data_{i+1}' for i in range(len(self.data_series))]
		
		if custom_colors:
			assert len(custom_colors) == len(self.data_series)
			colors = custom_colors
		else:
			colors = plt.get_cmap(color_palette)(np.linspace(0,1,len(self.data_series)))
		
		plt.figure(figsize=figure_size)
		for i in range(len(self.data_series)):
			x,y = self.data_series[i]
			plt.plot(x, y, color=colors[i], label=label_list[i])
		plt.xlabel(self.x_label)
		plt.ylabel(self.mode)
		plt.title(title)
		
		if xlim == 'auto':
			plt.xlim(np.min(x), np.max(x))
		else:
			plt.xlim(xlim)
		
		if ylim:
			plt.ylim(ylim)
		if legend:
			plt.legend(title=legend_title, frameon=False)
		plt.gca().spines['right'].set_color('none')
		plt.gca().spines['top'].set_color('none')
		if savefig:
			plt.savefig(savefig, dpi=300, bbox_inches='tight')
		plt.show()
  
	def plot_at(self, wavelength=350, figure_size=(6,6), title='', xlim='auto', ylim='auto', color='black', gradient=False, color1='blue', color2='red', lw=3):
		data = self.data_series

		wl_index = data[0][0].tolist().index(wavelength)
		abs_wl = [a[1][wl_index] for a in data]
		x = self.step()

		plt.figure(figsize=figure_size)
		if gradient:
			points = np.array([x, abs_wl]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			ccv = mcolors.ColorConverter().to_rgb
			colormap = make_colormap([ccv(color1), ccv(color2)])
			
			lc = LineCollection(segments, cmap=colormap, norm=plt.Normalize(0.0, 1.0), alpha=1.0)
			lc.set_array(np.linspace(0.0, 1.0, len(x)))
			lc.set_linewidth(lw)
			plt.gca().add_collection(lc)
		else:
			plt.plot(x, abs_wl, color=color, lw=lw)

		plt.xlabel(f'{self.series_measure} ({self.series_unit})')
		plt.ylabel(f'{self.mode} at {wavelength} nm')
		plt.title(title)

		if xlim == 'auto':
			plt.xlim(np.min(x), np.max(x))
		else:
			plt.xlim(xlim)
		
		if ylim == 'auto':
			y_max = np.max(abs_wl)
			gap = np.abs(y_max)*0.15
			plt.ylim(-gap, y_max+gap)
		else:
			plt.ylim(ylim)

		plt.gca().spines['right'].set_color('none')
		plt.gca().spines['top'].set_color('none')
		plt.show()

	def find_peak_wavelength(self, range_wl=(450, 600), pos=1):
		"""
		Find the wavelength of maximum absorbance within a given range.
		
		Parameters
		----------
		range_wl : tuple
			Wavelength range (min, max) to search within
		pos : int
			Position in data_series (1-indexed)
		
		Returns
		-------
		float
			Wavelength at maximum absorbance
		"""
		_wavelength = (self.data_series[pos-1][0]).tolist()
		_abs = (self.data_series[pos-1][1]).tolist()

		start = _wavelength.index(range_wl[0])
		end = _wavelength.index(range_wl[1])

		max_abs = np.max(_abs[start:end])
		max_idx = _abs.index(max_abs)
		max_wl = _wavelength[max_idx]

		return max_wl

	def plot_at_peak(self, range_wl=(450,600), pos=1, figure_size=(6,6), title='', xlim='auto', color='black', lw=3):
		"""Find peak wavelength and plot absorbance at that wavelength over time."""
		max_wl = self.find_peak_wavelength(range_wl=range_wl, pos=pos)
		print(f'Peak found at {max_wl} nm')
		self.plot_at(max_wl, figure_size=figure_size, title=title, xlim=xlim, color=color, lw=lw)
