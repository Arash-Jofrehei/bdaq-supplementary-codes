# Author: Arash Jofrehei - University of Zurich
# Pixel-wise crosstalk scan to be used for BDAQ53

from tqdm import tqdm
import numpy as np
import math
import os

from bdaq53.scan_base import ScanBase
from bdaq53.analysis import analysis
from bdaq53.analysis import plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# The parameters below will be written over the scan parameters and the chip parameters. Don't forget the mask file!

chip = 'irr-HLL-Bern'  # Arbitrary chip name. To be used in the name of the saved plots and numpy maps.
PA_BIAS = [350,380,410,450,500]
mainThrVcalRange = [500,680,2]      # ['VCAL_HIGH_start' , 'VCAL_HIGH_stop' , 'VCAL_HIGH_step']
secondThrVcalRange = [680,2780,20]  # ['VCAL_HIGH_start' , 'VCAL_HIGH_stop' , 'VCAL_HIGH_step']
vthLin = 372
krumenacher = 20
thrMapsAvailable = True
showPlots = False
maxXtalk = 20 # percent
binsPerPercent = 5 # for the correlation plot

injectionPatterns = ['firstPeak','secondPeakEven','secondPeakOdd'] # [first threshold , second threshold when injecting in even columns , second threshold when injecting in odd columns]

scan_parameters = {
    'start_column': 128,
    'stop_column': 264,
    'start_row': 0,
    'stop_row': 192,
    'maskfile': 'output_data/meta_tuning_HLL_irr_Bern_linear_380_thr1400_noise_occupancy_conservative.masks.h5',
    'use_good_pixels_diff': False,
    'VCAL_MED': 500,
    'VCAL_HIGH_start': 500,
    'VCAL_HIGH_stop': 1500,
    'VCAL_HIGH_step': 10
}


class ThresholdScan(ScanBase):
    scan_id = 'threshold_scan'
    VTH = 400
    PA = 410
    krum = 28
    vcalRange = range(500,1000,10)
    injectionPattern = 'firstPeak'
    def configure(self, start_column=0, stop_column=400, start_row=0, stop_row=192, TDAC=None, **_):
        self.chip.masks['enable'][start_column:stop_column, start_row:stop_row] = True
        
        start = 0
        step = 1
        if (self.injectionPattern is injectionPatterns[1]):
            start = 0
            step = 2
        if (self.injectionPattern is injectionPatterns[2]):
            start = 1
            step = 2
        
        self.chip.masks['injection'][start_column+start:stop_column:step, start_row:stop_row] = True
        self.chip.masks.apply_disable_mask()


        if TDAC is not None:
            if type(TDAC) == int:
                self.chip.masks['tdac'][start_column:stop_column, start_row:stop_row] = TDAC
            elif type(TDAC) == np.ndarray:
                self.chip.masks['tdac'][start_column:stop_column, start_row:stop_row] = TDAC[start_column:stop_column, start_row:stop_row]

        self.chip.masks.update(force=True)

    def scan(self, n_injections=100, VCAL_MED=500, VCAL_HIGH_start=1000, VCAL_HIGH_stop=4000, VCAL_HIGH_step=100, **_):

        vcal_high_range = range(VCAL_HIGH_start, VCAL_HIGH_stop, VCAL_HIGH_step)
        
        self.log.info('Starting scan...')
        pbar = tqdm(total=self.chip.masks.get_mask_steps() * len(vcal_high_range), unit=' Mask steps')
        for scan_param_id, vcal_high in enumerate(vcal_high_range):
            self.chip.setup_analog_injection(vcal_high=vcal_high, vcal_med=VCAL_MED)
            self.chip.registers['Vthreshold_LIN'].write(value=self.VTH)
            self.chip.registers['PA_IN_BIAS_LIN'].write(value=self.PA)
            self.chip.registers['KRUM_CURR_LIN'].write(value=self.krum)
            with self.readout(scan_param_id=scan_param_id):
                for fe in self.chip.masks.shift(masks=['enable', 'injection'], cache=True):
                    if not fe == 'skipped' and fe == 'SYNC':
                        self.chip.inject_analog_single(send_ecr=True, repetitions=n_injections)
                    elif not fe == 'skipped':
                        self.chip.inject_analog_single(repetitions=n_injections)
                    pbar.update(1)

        pbar.close()
        self.log.success('Scan finished')

    def analyze(self, create_pdf=True):
        with analysis.Analysis(raw_data_file=self.output_filename + '.h5') as a:
            a.analyze_data()
            if not os.path.exists('crosstalkNPY'):
                os.mkdir('crosstalkNPY')
                print('crosstalkNPY directory created.')
            fileName = 'crosstalkNPY/'+chip + '_PA-BIAS-' + str(self.PA) + '_' + self.injectionPattern
            np.save(fileName,a.threshold_map)
            
            if create_pdf:
                with plotting.Plotting(analyzed_data_file=a.analyzed_data_file) as p:
                    p.create_standard_plots()
            
            return a.threshold_map


def performThrScan(PA,injectionPattern):
    if injectionPattern is injectionPatterns[0]:
        scan_parameters['VCAL_HIGH_start'] = mainThrVcalRange[0]
        scan_parameters['VCAL_HIGH_stop'] = mainThrVcalRange[1]
        scan_parameters['VCAL_HIGH_step'] = mainThrVcalRange[2]
    else:
        scan_parameters['VCAL_HIGH_start'] = secondThrVcalRange[0]
        scan_parameters['VCAL_HIGH_stop'] = secondThrVcalRange[1]
        scan_parameters['VCAL_HIGH_step'] = secondThrVcalRange[2]
    with ThresholdScan(scan_config=scan_parameters) as scan:
        scan.VTH = vthLin
        scan.krum = krumenacher
        scan.PA = PA
        scan.injectionPattern = injectionPattern
        scan.start()
        interpretedFile = scan.output_filename + '_interpreted'
        renamedFile = 'output_data/' + chip + '_PA' + str(PA) + '_xtalkScan_' + injectionPattern
        thrMap = scan.analyze()
        os.rename(interpretedFile+'.h5',renamedFile+'.h5')
        os.rename(interpretedFile+'.pdf',renamedFile+'.pdf')
        if injectionPattern is injectionPatterns[2]:
            print('                    ----    ----    ----    ----')
            print('Threshold scans of PA-BIAS',PA,'are finished.')
            print('                    ----    ----    ----    ----')
        return thrMap


def getThrMaps(PA):
    fileNameBase = 'crosstalkNPY/'+chip + '_PA-BIAS-' + PA + '_'
    return [np.load(fileNameBase+'firstPeak.npy'),np.load(fileNameBase+'secondPeakEven.npy'),np.load(fileNameBase+'secondPeakOdd.npy')]

def xtalkCalc(patternThrMaps):
    fP = 10.4 * patternThrMaps[0] + 180
    sPE = 10.4 * patternThrMaps[1] + 180
    sPO = 10.4 * patternThrMaps[2] + 180
    xtalk = np.zeros((192,400))
    nPair = int((264-128)/2)
    nBins = maxXtalk*binsPerPercent
    pairCorrelation = np.zeros((nBins,nBins))
    for row in range(0,192):
        for pair in range(nPair):
            col1 = 2*pair + 128
            col2 = col1 + 1
            r1 = fP[col1,row]/sPO[col1,row]
            r2 = fP[col2,row]/sPE[col2,row]
            xtalk[row,col1] = 100*(r2-r1*r2)/(1-r1*r2)
            xtalk[row,col2] = 100*(r1-r1*r2)/(1-r1*r2)
            if (xtalk[row,col1] > maxXtalk):
                #xtalk[row,col1] = maxXtalk
                xtalk[row,col1] = 0
            if (xtalk[row,col2] > maxXtalk):
                #xtalk[row,col2] = maxXtalk
                xtalk[row,col2] = 0
            if (math.isnan(xtalk[row,col1]) or math.isnan(xtalk[row,col2])):
                xtalk[row,col1] = 0
                xtalk[row,col2] = 0
            col1Bin = max(0,min(nBins-1,int(xtalk[row,col1]*binsPerPercent)))
            col2Bin = max(0,min(nBins-1,int(xtalk[row,col2]*binsPerPercent)))
            pairCorrelation[col2Bin,col1Bin] += 1

    return xtalk,pairCorrelation

def plot(xtalk,pairCorrelation,PA):
    if not os.path.exists('plots'):
        os.mkdir('plots')
        print('plots directory created.')
    if not os.path.exists('plots/crosstalk'):
        os.mkdir('plots/crosstalk')
        print('crosstalk plot directory created.')
    plt.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(13, 11))
    ax = fig.add_subplot(111)
    ax.set_title('crosstalk distribution for '+ chip +' with PA bias '+ PA +'\n')
    plt.xlabel('crosstalk (%)')
    plt.ylabel('number of pixels in the linear FE')
    plt.hist(xtalk[:,128:264].reshape(((264-128)*192)),100,range=(0.0001, maxXtalk))
    plt.savefig('plots/crosstalk/crosstalkDist_'+ chip +'_PA'+ PA)
    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111)
    ax.set_title('crosstalk map for '+ chip +' with PA bias '+ PA +'\n')
    plt.xlabel('RD53A column')
    plt.ylabel('RD53A row')
    plt.imshow(xtalk[:,128:264],extent=[128,264,192,0])
    ax.set_aspect('equal')
    cax = fig.add_axes([0.12, 0.11, 0.79, 0.77])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    maskedMap = xtalk[:,128:264].reshape(((264-128)*192))[np.nonzero(xtalk[:,128:264].reshape(((264-128)*192)))]
    plt.clim(np.quantile(maskedMap,0.01), np.quantile(maskedMap,0.99))
    plt.savefig('plots/crosstalk/crosstalkMap_'+ chip +'_PA'+ PA)
    fig = plt.figure(figsize=(14, 11))
    ax = fig.add_subplot(111)
    ax.set_title('crosstalk correlation for '+ chip +' with PA bias '+ PA +'\n')
    plt.imshow(pairCorrelation,origin='lower',norm=LogNorm(),extent=[0,maxXtalk,0,maxXtalk])
    ax.set_aspect('equal')
    plt.xlabel('xtalk from top to bottom (%)')
    plt.ylabel('xtalk from bottom to top (%)')
    cax = fig.add_axes([0.12, 0.11, 0.82, 0.77])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig('plots/crosstalk/crosstalkCor_'+ chip +'_PA'+ PA)
    print(' **** **** '+chip+' with PA bias '+PA+' is plotted. **** ****')

if __name__ == '__main__':
    for amp in range(len(PA_BIAS)):
        if thrMapsAvailable:
            patternThrMaps = getThrMaps(str(PA_BIAS[amp]))
        else:
            patternThrMaps = []
            for pat in range(len(injectionPatterns)):
                patternThrMaps.append(performThrScan(PA_BIAS[amp],injectionPatterns[pat]))
        xtalkMap,pairCorrelation = xtalkCalc(patternThrMaps)
        plot(xtalkMap,pairCorrelation,str(PA_BIAS[amp]))
    if showPlots:
        plt.show()
        
