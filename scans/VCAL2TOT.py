# Author: Arash Jofrehei - University of Zurich
# Pixel-wise crosstalk scan to be used for BDAQ53

from tqdm import tqdm
from bdaq53.scan_base import ScanBase
from bdaq53.analysis import analysis
from bdaq53.analysis import plotting
from bdaq53.analysis import online as oa
import bdaq53.analysis.analysis_utils as au
import numpy as np
import h5py
import tables as tb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# The parameters below will be written over the scan parameters and the chip parameters.

masks = ['output_data/meta_tuning_HLL_irr_Bern_linear_372_thr1123_noise_occupancy_conservative.masks.h5','output_data/meta_tuning_HLL_irr_Bern_linear_380_thr1400_noise_occupancy_conservative.masks.h5']
krumList = [20,24,28]
vthList = [372,380]
thrList = ['372_1123','380_1400']
ROC = 'irr-HLL-Bern' # Arbitrary chip name. To be used in the name of the saved plots and numpy maps.
thrFile = False
#thrFile = 'output_data/long-threshold-scan_irr-HLL-Bern_VTH372_thr1123_HV300_krum20_interpreted.h5'
stepVcal = 10
nInjections = 100
precisionPerTot = 5
availableMap = True
showPlots = False    # plots are saved anyways


scan_parameters = {
    'start_column': 128,
    'stop_column': 264,
    'start_row': 0,
    'stop_row': 192,
    'maskfile': None,
    'use_good_pixels_diff': False,
    'VCAL_MED': 500,
    'VCAL_HIGH_start': 500,
    'VCAL_HIGH_stop': 2850,
    'VCAL_HIGH_step': 10
}


class ThresholdScan(ScanBase):
    scan_id = 'threshold_scan'
    krum = 28
    VTH = 400
    def configure(self, start_column=0, stop_column=400, start_row=0, stop_row=192, TDAC=None, **_):
        self.chip.masks['enable'][start_column:stop_column, start_row:stop_row] = True
        self.chip.masks['injection'][start_column:stop_column, start_row:stop_row] = True
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
            return a.analyzed_data_file



def performThrScan(mask,VTH,krum,stepVcal):
    scan_parameters['VCAL_HIGH_step'] = stepVcal
    scan_parameters['maskfile'] = mask
    with ThresholdScan(scan_config=scan_parameters) as scan:
        scan.VTH = VTH
        scan.krum = krum
        scan.start()
        return scan.analyze()

def getOcc(inputFile):
    with h5py.File(inputFile, "r") as f:
        HistOcc = f['HistTot'][()]
        return HistOcc

def getTotMean(occ,nInjections):
    noiseORxtalkHITS = np.sum(occ) - nInjections
    if (noiseORxtalkHITS > 20):
        minTOT = np.min(occ[occ.astype(bool)])
        noise_arg = np.where(occ==minTOT)[0][0]
        occ[noise_arg] = max(0,occ[noise_arg]-noiseORxtalkHITS)
    meanTot = np.sum(np.multiply(np.arange(16),occ))*1.0/max(nInjections,np.sum(occ))
    return meanTot

def getNpyFiles(VCALFile,totFile):
    return np.load(VCALFile),np.load(totFile)

def plotting(VCALs,meanTotMap,thr,krum,precisionPerTot):
    nSteps = len(VCALs)
    plotLog = np.zeros((precisionPerTot*15,nSteps))
    for col in range(128,264):
        for row in range(0,192):
            #if (row == 80 and col%80):
                #plt.plot(VCALs,meanTotMap[col,row,:])
            for v in range(nSteps):
                plotLog[int(precisionPerTot*meanTotMap[col,row,v]),v] += 1
    fig = plt.figure(figsize=(13,11))
    ax = fig.add_subplot(111)
    ax.set_title('TOT vs VCAL for chip '+ROC+' at Vthreshold_threshold '+thr+' and krumenacher '+krum)
    plt.xlabel('delta VCAL')
    plt.ylabel('TOT')
    plt.imshow(plotLog,origin='lower',norm=LogNorm(),extent=[0,VCALs[-1],0,15])
    ax.set_aspect(int(VCALs[-1]/15))
    cax = fig.add_axes([0.12, 0.11, 0.85, 0.77])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig('plots/calibrationVCAL/TOTvsVCAL_'+ ROC +'_VTH-thr_'+ thr + 'krum_' + krum)
    if showPlots:
        plt.show()

def v2tAnalysis(thrFile,stepVcal,nInjections):
    occMap = getOcc(thrFile)
    nSteps = occMap.shape[2]
    VCALs = stepVcal*np.arange(nSteps)
    meanTotMap = np.zeros((400,192,nSteps))
    for col in range(128,264):
        for row in range(192):
            for step in range(nSteps):
                meanTotMap[col,row,step] = getTotMean(occMap[col,row,step,:],nInjections)
    return VCALs, meanTotMap

    

if __name__ == '__main__':
    THR = thrList[0]
    KRUM = str(krumList[0])
    VCALFile = 'VCALs_'+ROC+'_'+THR+'_'+KRUM+'.npy'
    totFile = 'mean_tot_'+ROC+'_'+THR+'_'+KRUM+'.npy'
    if not thrFile:
        takeThreshold = True


    if availableMap:
        for thr in range(len(vthList)):
            for krum in range(len(krumList)):
                THR = thrList[thr]
                KRUM = str(krumList[krum])
                VCALFile = 'calibrationVCAL/VCALs_'+ROC+'_'+THR+'_'+KRUM+'.npy'
                totFile = 'calibrationVCAL/mean_tot_'+ROC+'_'+THR+'_'+KRUM+'.npy'
                VCALs,meanTotMap = getNpyFiles(VCALFile,totFile)
                plotting(VCALs,meanTotMap,thrList[thr],str(krumList[krum]),precisionPerTot)
    else:
        for thr in range(len(vthList)):
            for krum in range(len(krumList)):
                if takeThreshold:
                    thrFile = performThrScan(masks[thr],vthList[thr],krumList[krum],stepVcal)
                THR = thrList[thr]
                KRUM = str(krumList[krum])
                VCALFile = 'calibrationVCAL/VCALs_'+ROC+'_'+THR+'_'+KRUM+'.npy'
                totFile = 'calibrationVCAL/mean_tot_'+ROC+'_'+THR+'_'+KRUM+'.npy'
                VCALs,meanTotMap = v2tAnalysis(thrFile,stepVcal,nInjections)
                # saving files
                np.save(totFile,meanTotMap)
                np.save(VCALFile,VCALs)
                print('                    ----    ----    ----    ----')
                print('Calibration files for chip <<',ROC,'>> at Vthreshold_threshold',thrList[thr],'and krumenacher',krumList[krum],'are saved.')
                print('                    ----    ----    ----    ----')
                plotting(VCALs,meanTotMap,thrList[thr],str(krumList[krum]),precisionPerTot)
