
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
# import cartopy.crs as ccrs
import subprocess
import os
from ftplib import FTP_TLS
import platform
import cdflib
import glob
from pathlib import Path

import datetime
from numpy.lib.arraypad import pad

from tqdm import tqdm

#https://notebook.community/daniestevez/jupyter_notebooks/IONEX
#https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html

def ionex_filename_to_date(filename):
    if os.sep in filename: 
        filename = filename.split(os.sep)[-1]
    a = re.findall(r'\d+', filename)
    day_number = a[0][:-1]
    year = '20'+a[1]
    # print(day_number,year,a)
    return day_number_to_date(day_number,year)

def day_number_to_date(day_number,year):
    day_number = str(day_number)
    year = str(year)
    # print(day_number,year)
    res = datetime.datetime.strptime(year + "-" + day_number, "%Y-%j")
    return res

def date_to_year_day_number(year,month,day):
	today = datetime.date(year, month, day)
	return int(today.strftime('%j'))

def time_to_index_day_number(hours,minutes,step_in_minutes):
	h_step = int(60 / step_in_minutes)
	h = hours * h_step
	m = int(minutes / step_in_minutes)
	return int(h+m)

def index_day_number_to_time(index,step_in_minutes):
	h_step = 60 / step_in_minutes
	h = int(index / h_step) % 24
	m = int((index % h_step) * step_in_minutes)

	return h,m

class IONEX(object):

	def __init__(self, save_directory , centre = 'esa'):
		self.centre = centre
		self.directory = save_directory
	
	def parse_map(self,tecmap, exponent = -1):
		tecmap = re.split('.*END OF TEC MAP', tecmap)[0]
		return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n',tecmap)[1:]])*10**exponent
	
	def parse_rms(self,tecmap, exponent = -1):
		tecmap = re.split('.*END OF RMS MAP', tecmap)[0]
		return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n',tecmap)[1:]])*10**exponent

	def get_numpy_tecmaps(self,years_list,days_list):
		if not type(days_list) == list: days_list = [days_list]
		if not type(years_list) == list: years_list = [years_list]
		result = None
		file_names = []
		# for year in tqdm(years_list,desc="Year",position=0, leave=True):
		for year in years_list:
			for day in tqdm(days_list,desc="Year : {}, Days : ".format(year),position=0, leave=True):
				try:
					file_name = self.ionex_local_path(year,day)
					np_tmap = np.array(self.get_tecmaps(file_name))
					file_names.append(file_name)
				except Exception as e:
					print(e,file_name)
					continue
				if np_tmap.shape[0] > 13:
					np_tmap = np_tmap[list(range(0,26,2)),:,:]
				d,h,w = np_tmap.shape
				np_tmap = np_tmap.reshape(1,d,h,w)

				#print(np_tmap.shape)
				if result is None:
					result = np_tmap.copy()
				else:
					result = np.vstack((result,np_tmap))
				# print(result.shape)
		return result,file_names

	def get_numpy_rmsmaps(self,years_list,days_list):
		if not type(days_list) == list: days_list = [days_list]
		if not type(years_list) == list: years_list = [years_list]
		result = None
		file_names = []
		# for year in tqdm(years_list,desc="Year",position=0, leave=True):
		for year in years_list:
			for day in tqdm(days_list,desc="Year : {}, Days : ".format(year),position=0, leave=True):
				try:
					file_name = self.ionex_local_path(year,day)
					np_tmap = np.array(self.get_rmsmaps(file_name))
					file_names.append(file_name)
				except Exception as e:
					print(e,file_name)
					continue
				if np_tmap.shape[0] > 13:
					np_tmap = np_tmap[list(range(0,26,2)),:,:]
				# print(file_name)
				d,h,w = np_tmap.shape
				np_tmap = np_tmap.reshape(1,d,h,w)

				#print(np_tmap.shape)
				if result is None:
					result = np_tmap.copy()
				else:
					result = np.vstack((result,np_tmap))
				# print(result.shape)
		return result,file_names

	def get_tecmaps(self,filename):
		with open(filename) as f:
			ionex = f.read()
			return [self.parse_map(t) for t in ionex.split('START OF TEC MAP')[1:]]

	def get_rmsmaps(self,filename):
		with open(filename) as f:
			ionex = f.read()
			return [self.parse_rms(t) for t in ionex.split('START OF RMS MAP')[1:]]

	def get_tec(self,tecmap, lat, lon):
		i = round((87.5 - lat)*(tecmap.shape[0]-1)/(2*87.5))
		j = round((180 + lon)*(tecmap.shape[1]-1)/360)
		return tecmap[i,j]

	def ionex_filename(self,year, day, zipped = True):
		return '{}g{:03d}0.{:02d}i{}'.format(self.centre, day, year % 100, '.Z' if zipped else '')

	def ionex_ftp_path(self,year, day):
		"""
		'gps/products/ionex/2010/001/esag0010.10i.Z'
		"""
		return '/gps/products/ionex/{:04d}/{:03d}/{}'.format(year, day, self.ionex_filename(year, day))
		# return 'ftp://cddis.gsfc.nasa.gov/gps/products/ionex/{:04d}/{:03d}/{}'.format(year, day, self.ionex_filename(year, day, self.centre))

	def ionex_local_path(self,year, day, zipped = False):
		return os.path.join(self.directory,str(year),self.ionex_filename(year, day, zipped))

	def create_dir_path(self,filename):
		if not os.path.exists(os.path.dirname(filename)):
			try:
				os.makedirs(os.path.dirname(filename))
			except OSError as exc: # Guard against race condition
				pass
		
	def download_ionex(self, year, day, debug=False):
		if not type(year) == list:
			year = [year]
		if not type(day) == list:
			day = [day]

		ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
		ftps.login(user='anonymous', passwd='1234@gmail.com')
		ftps.prot_p()

		for y in year:
			for d in day:
				ftp_path = self.ionex_ftp_path(y, d)

				filename_zip = os.path.join(self.directory,str(y),'zip',self.ionex_filename(y, d))
				filename = os.path.join(self.directory,str(y),self.ionex_filename(y, d))[:-2]

				self.create_dir_path(filename_zip)
				self.create_dir_path(filename)

				if not os.path.isfile(filename_zip) or os.path.getsize(filename_zip) < 1:
					try:
						os.remove(filename_zip)
					except: pass
					if debug : print('Downloading... : ',filename_zip)
					try:
						ftps.retrbinary("RETR " + ftp_path, open(filename_zip, 'wb').write)
					except OSError as ose:
						print('OSError : {}, Another attempt!'.format(ose))
						ftps.retrbinary("RETR " + ftp_path, open(filename_zip, 'wb').write)
						pass
				else:
					if debug : print(filename_zip,' exist!')

				if not os.path.isfile(filename):
					if debug : print('Extracting... : ',filename)
					subprocess.call(['7z', 'e', filename_zip,'-o'+os.path.dirname(filename)])
					
					# if platform.system() == 'Windows':
					# 	subprocess.call(['7z', 'e', filename_zip,'-o'+os.path.dirname(filename)])
					# else:
					# 	subprocess.call(['gzip', '-d', filename])
				else:
					if debug : print(filename,' exist!')

		ftps.close()

		if debug : print('DONE!')
		
	# def plot_tec_map(self,tecmap):
	# 	proj = ccrs.PlateCarree()
	# 	f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=proj))
	# 	ax.coastlines()
	# 	h = plt.imshow(tecmap, cmap='viridis', vmin=0, vmax=100, extent = (-180, 180, -87.5, 87.5), transform=proj)
	# 	plt.title('VTEC map')
	# 	divider = make_axes_locatable(ax)
	# 	ax_cb = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)

	# 	f.add_axes(ax_cb)
	# 	cb = plt.colorbar(h, cax=ax_cb)
	# 	plt.rc('text', usetex=True)
	# 	cb.set_label('TECU ($10^{16} \\mathrm{el}/\\mathrm{m}^2$)')	
		

class IONEX_CDF(object):
    
	def __init__(self,save_directory):
    	
		self.directory = save_directory

	def get_numpy_tecmaps(self,years_list):
    		
		days_array = []
		
		for year in years_list:
    	
			fpath = os.path.join(self.directory,str(year),'*.cdf')
			cdf_files_list = glob.glob(fpath)

			
			cdf_files_list.sort()
			for file in cdf_files_list:#,desc="CDF {} : ".format(fpath)):
				# print(file)
				cdf = cdflib.CDF(file)
				# print(cdf.varget('lat'),cdf.varget('lon'))


				days_array.append(cdf.varget('tecUQR'))
				
			
		days_array = np.array(days_array)
		return days_array



def map_matrix2string_ionex(matrix,inonex_map_columns=16):
    n=inonex_map_columns
    test_map_str = [["   "+'   '.join(map(str, arr[i:i + n])) for i in range(0, len(arr), n)] for arr in matrix]
    test_map_str = ['\n'.join(map(str, l)) for l in test_map_str]
    return test_map_str


def dmd_ionex(c1p_file_path,dmd_predicted_maps,_replace='c1p',_replace_with='dmd',debug=False):
    data = ''
    file_path = c1p_file_path
    maps = dmd_predicted_maps[0]
    with open(file_path,'r') as f:
        line  = f.readline()
        map_count = 0
        lat_count = 0
        test_map_str = map_matrix2string_ionex(maps[map_count])
        while line:
            
            data += line
            if "LAT/LON1/LON2/DLON/H" in line:
                try:
                    for _ in range(5):
                        line  = f.readline()
                    data += test_map_str[lat_count]+"\n"
                    lat_count+=1
                except IndexError as ie:
                    if debug: print('ERROR : ',map_count,lat_count,line)
                    pass
            if "END OF TEC MAP" in line:
                map_count+=1
                lat_count=0
                try:
                    test_map_str = map_matrix2string_ionex(maps[map_count])
                except:
                     if debug: print('ERROR : ',map_count,maps.shape)

            line  = f.readline()
            if "START OF RMS MAP" in line:
                break
        while line:
            data += line
            line  = f.readline()

    dmd_file = file_path.replace(_replace,_replace_with)

    if debug : print(dmd_file,os.path.dirname(dmd_file))
    Path(os.path.dirname(dmd_file)).mkdir(parents=True, exist_ok=True)
    with open(dmd_file,'w') as f:
        f.write(data)
    
	


if __name__ == '__main__':

	# from utilities.IONEX import IONEX,IONEX_15MIN
	YEAR = 2017
	DAYS = list(range(1,366))
	# ionex_gps = IONEX(save_directory=os.path.join(os.path.abspath('.'),'ionex'),centre='cod')

	# ionex_gps.download_ionex(YEAR,DAYS)

	ionex_igs15 = IONEX_CDF(save_directory=os.path.join(os.path.abspath('.'),'igs_15min'))

	# igs15 = ionex_igs15.get_numpy_tecmaps([2013,2014])

	day = date_to_year_day_number(2013,12,31)+date_to_year_day_number(2014,6,10)
	time = time_to_index_day_number(11,00,15)
	print(day,time)
	# print(date_to_year_day_number(2013,12,31)+date_to_year_day_number(2014,6,10))
	# print(time_to_index_day_number(16,00,15)-time_to_index_day_number(11,00,15))
	igs15 = ionex_igs15.get_numpy_tecmaps([2017])
	print(igs15.shape)


	fig,axs = plt.subplots(4,5,figsize=(8,8))
	plt.subplots_adjust(bottom=12)
	
	for r,raxs in enumerate(axs):
		for c,ax in enumerate(raxs):
			# ax.imshow(np.eye(256))
			# print(index_day_number_to_time(time,15))
			h,m = index_day_number_to_time(time,15)
			ax.set_title("{}:{}".format(format(h, '02d'),format(m, '02d')))
			im = ax.imshow(igs15[day,time])
			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			plt.colorbar(im, cax=cax)
			time+=1
			
	# ax1.imshow(igs15[day,time-1])
	# ax2.imshow(igs15[day,time])
	# ax3.imshow(igs15[day,time+1])
	# ax1.imshow(np.eye(256))
	plt.show()