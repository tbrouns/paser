# PASER

MATLAB toolbox for processing and analyzing extracellular recordings, including local field potential (LFP) and spiking data.  

Currently, PASER can only be used with data saved by Open Ephys GUI [Ref. 4], specifically `.continuous` files. 
See https://github.com/open-ephys/plugin-GUI or http://www.open-ephys.org/ for more information.

PASER contains the following types of processing routines:

* Spike detection and sorting
* Cluster quality control
* Spike analysis 
* LFP detection and analysis
* Data visualization 
* Artifact removal

# Getting started

## Install third-party dependencies

In order to use the toolbox, a couple of third-party software packages need to be installed on your system.

### FieldTrip 

We use the FieldTrip toolbox [Ref. 3] for LFP processing and analysis, which can be downloaded or cloned from:

https://github.com/fieldtrip/fieldtrip

If you are not planning on using FieldTrip for anything else, then do not add the FieldTrip toolbox to your MATLAB path. 
This will be done at a later step. Otherwise, follow the directions here:

http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path

### KiloSort

The default spike sorting method in PASER is KiloSort [Ref. 1], which can be downloaded or cloned from: 

https://github.com/cortex-lab/KiloSort

It is highly recommended to use KiloSort with a CUDA enabled GPU. Attempting to run KiloSort on the CPU is errorprone and not guaranteed to result in satisfactory cluster quality. 
Therefore, please follow the installation instructions given in the KiloSort README file and "Docs" folder. 

Once again, if you are not using KiloSort for anything else, then do not add the KiloSort directory to the MATLAB path. 
We will instead load the toolbox the moment it is needed in the data processing pipeline. 

## PASER toolbox installation

Clone or download PASER and add it to your MATLAB path: 

```
addpath(genpath('C:\Path\To\paser-master'));
```

## Quick start

### Default processing script

Create a MATLAB script containing the following lines of code:

```
parameters = [];
parameters.subject      = 'SubjectName'; % Name of subject used in output MAT filename (no spaces)
parameters.loadPath     = 'C:\PathToLoadFolder\LoadFolder\'; % Where the data folders are
parameters.savePath     = 'D:\PathToSaveFolder\SaveFolder\'; % Where you want to save the output MAT files
parameters.configPath   = 'E:\PathToConfigFile\ConfigFile; % Where the parameters are loaded from
parameters.txtfile      = [];            % Folders to process given in text file
parameters.patterns     = [];            % Used to differentiate between experimental sessions (string cell array)
parameters.type         = 'all';         % Which session type to process ('all' or one of the chosen patterns)
parameters.process      = 'new';         % Which specific sessions to process ('new', 'given', 'from' or 'all')
parameters.folders      = [];            % Sessions that you wish to process, if 'given' is chosen (string cell array)
parameters.extension    = 'continuous';  % File extension of raw data files
parameters.filepattern  = 'CH';          % Pattern to look for in data files
parameters.blockpattern = [];            % Used to differentiate between blocks within session
parameters.nelectrodes  = 4;             % Number of electrodes per polytrode (e.g. tetrode: 4)

psr_batch_processing(parameters); % Process raw data files
```

This script will load and then batch process `continuous` files in the directory given by `parameters.loadPath`. 
Processed data is saved to the folder specified by `parameters.savePath`, where a matching directory tree will be created. 

Not all of the parameters in this script are set correctly, so what follows is an explanation of how to select the right settings.

### Initial parameters

#### Config file

You should create a `ConfigFile.m` that contains parameter settings for the various processing functions and then point `parameters.configPath` to this file. 
This script should at least contain the following lines of code:

```
psr_parameters_general;
parameters.path.kst = 'C:\PathToKiloSort\KiloSort';   % Path to the KiloSort repository
parameters.path.ft  = 'D:\PathToFieldTrip\FieldTrip'; % Path to the FieldTrip repository
```

More `parameters` fields can be changed by adding more lines to the script. For example, if you do not want to process the LFP, you add:
`parameters.process.lfp = false;`

See `psr_parameters_general` for all other parameter fields. Comments next to each parameter explain its purpose. 

#### Path settings

To avoid clogging up the MATLAB path, we add the third-party toolboxes from the path whenever we need them and remove them afterwards. 
In order for the program to know where to look for the toolboxes, set the path parameters in the following way:

* `parameters.path.kst` should point to the KiloSort main directory, which contains e.g. `preprocessData.m`
* `parameters.path.ft` should point to the FieldTrip main directory, which contains e.g. `ft_defaults.m`

#### Load directory

`parameters.loadPath` should point to a directory tree like the one given below. 

```
loadPath
│
├───Session_1
│   │
│   ├───Block_0M
│   │   │   100_ADC1.continuous
│   │   │   100_ADC2.continuous
│   │   │   ...
│   │   │
│   │   │   100_CH1.continuous
│   │   │   100_CH2.continuous
│   │   │   ...
│   │
│   │...
│
├───Session_1_Condition
│   │   
│   ├───Block_0M
│   │   │   100_ADC1.continuous
│   │   │   100_ADC2.continuous
│   │   │   ...
│   │   │
│   │   │   100_CH1.continuous
│   │   │   100_CH2.continuous
│   │   │   ...
│   │
│   ├───Block_1M
│   │   │   100_ADC1.continuous
│   │   │   100_ADC2.continuous
│   │   │   ...
│   │   │
│   │   │   100_CH1.continuous
│   │   │   100_CH2.continuous
│   │   │   ...
│   │
│   │...
│
│...
```

The `loadPath` folder should include one or more folders for specific experimental sessions. Here, we have two such folders called `Session_1` and `Session_1_Condition`.
Each session folder then contains one or more folders for specific blocks that hold the raw data files. 
In the case of `Session_1_Condition`, we have the `Block_0M` and `Block_1M`, which contain the `continuous` files. 
To be clear, the names `Session` and `Block` are arbitrary here, you can use any other name you desire.

However, we must note that we have two types of `continuous` files here. We only wish to load the `100_CH*.continuous` ones. 
We can differentiate between the two types using `parameters.filepattern`, which should be set to `parameters.filepattern = 'CH'` in this case, because that is the common pattern between them. 

The pattern specified by `parameters.filepattern` must immediately be followed by an integer in the filename. In turn, the integer should immediately be followed by the file extension or by an underscore. 
This integer is used to determine which channels belong to which probe. For example, if we are using a tetrode, then the channels `{*CH1,*CH2,*CH3,*CH4}.continuous` are the channels for the first tetrode.

Depending on the experimental conditions, you may wish to vary some kind of experimental variable across different blocks. 
If you indicate the value of the variable in the block folder name, then it will be extracted and saved by setting `parameters.blockpattern`. 
Any value between the underscore and the specified pattern is recorded. If we set `parameters.blockpattern = 'M'`, we would get 0 and 1 for the blocks in `Session_1_Condition`. 

Furthermore, you can also differentiate between different session types and only process one particular type of session. 
First set the different session types in `parameters.patterns` and in `parameters.type` which type you wish to process. 
If we only want to process the 'Condition' sessions, we would have to specify `parameters.patterns = {'condition'}` and `parameters.type = 'condition'` (not case sensitive). 

The session folders you wish to process can be directly specified as well by setting `parameters.process` to `'given'` 
and giving the session names as a cell array in `parameters.folders` (e.g. `parameters.folders = {'Session_1_Condition'}`). 
You can also set `parameters.process` to `'from'` to process the data starting from a specific session, where the starting session is specified in `parameters.folders` as a cell. 
Otherwise, `parameters.process` can be set to `'new'` when you only want to process sessions for which no output MAT files exist in the `savePath`, or `'all'` if you want to process every single session and overwrite any existing data.

Lastly, ... [INCLUDE INFO ABOUT TXTFILE FIELD]

#### Note on blocks

We assume that blocks occur one after another during the experiment, which means that the electrodes will remain in the same position. 
Therefore, we always perform spike sorting across different blocks for each probe within a session.

#### Polytrode channels

The number of channels of the polytrode should be indicated by `parameters.nelectrodes`.
Each block folder should then contain a number of `continuous` files equal to the number of polytrodes multiplied by `parameters.nelectrodes`.
 
#### Save directory

For `parameters.savePath` select the folder where you want to save the output MAT files. Folders are automatically created to match the `loadPath` directory tree. 
So, in the example given above, the folders `Session_1` and `Session_1_Condition` will be created in `savePath`, which will contain the output MAT files for the corresponding sessions. 

### Output files

As mentioned earlier, a MAT file will be saved for each probe to the `savePath` for the current session. These files have the following naming convention:

`PSR_%SubjectName%_%Session%_P%ProbeNumber%.mat`

Where `%SubjectName%` is specified by `parameters.subject`, `%Session%` by the name of the session folder and `%ProbeNumber%` by the probe number.

The data in every MAT file is organized in the following way:

```
Structures:     Fields:	       Type:    Size:            Description:

freq                           struct                    Local field potential (LFP) data                         
                artifacts      double   [Na x  2]        Start and end point of LFP artifacts [sec]
                cfg            struct                    FieldTrip configuration parameters
				fsample        double   scalar           Sampling frequency
				hdr            struct                    FieldTrip parameters
                label          cell     [ 1 x Nc]        Cell array of label for each channel, where each cell contains:
                               string                    Channel label
				sampleinfo     double   [Nt x  2]        Start and end points of trials 
				time           cell     [ 1 x Nt]        Cell array of timestamps for each trial, where each cell contains:			   
							   double   [ 1 x Ns]        Timestamps of specific trial [sec]
                trial          cell     [ 1 x Nt]        Cell array of voltage values for each trial, where each cell contains:
                               double   [Nc x Ns]        Voltage values of specific trial [muV]
				
                [ see FieldTrip docs for details: http://www.fieldtriptoolbox.org/ ]

                Na: number of artifacts	
                Nc: number of channels
				Ns: number of data points 
                Nt: number of trials
				
metadata                       struct                    General experimental data
                duration       double   scalar           Duration of session [sec]
				probe          integer  scalar           Probe number
				session        string   [ 1 x Ns]        Name of session
				stimtimes      cell     [Nb x  2]        Cell array with stimulus on- and offset timings, where each cell contains:
                               double   [Nt x  2]        Stimulus on- and offset times [sec]
                               string                    Stimulus type
                stimulus       cell     [Nt x  1]        Vector of trial conditions
                subject        string                    Name of subject
                trialonset     double   [Nb x  1]        Onset time of every block [sec]
				
				Nb: number of blocks
				Ns: number of sessions
                Nt: number of trials
				
parameters                     struct                    Parameters for all data processing functions
                
                [ see files in "parameters" folder for details ]
			
spikes                         struct                    Neural spiking data
                assigns        int16    [ 1 x Ns]        Cluster index of each detected spike after merging
                assigns_prior  int16    [ 1 x Ns]        Cluster index of each detected spike before merging
                clusters       struct                    See further below
				delete         struct                    Logical arrays indicating spikes that are tagged for removal
				features       single   [Nd x Ns]        Array of principle component scores
                info           struct                    See further below
                spiketimes     single   [ 1 x Ns]        Spike time of each detected spike [sec]
                trials         int16    [ 1 x Ns]        Trial index of each detected spike
                waveforms      int16    [Ns x Np x Nc]   Waveform of each detected spike for each channel 
				Fs             double   scalar           Sampling frequency of raw extracellular recording
				
                Ns: number of spikes
                Np: number of data points
                Nc: number of channels
				
spikes.info                    struct                    Session information
                std            double   [1 x Nc]         Standard deviation of signal for each channel 
				mad            double   [1 x Nc]         Median absolute deviation of signal for each channel
				rms            double   [1 x Nc]         Root-mean-square of signal for each channel
				env            double   [1 x Nc]         Signal envelope for each channel
				bgn            double   [1 x Nc]         Measure of background noise for each channel
                dur            double   [1 x Nb]         Duration of each experimental block [sec]				
                thresh         double   [1 x Nc]         Spike detection thresholds for each channel
				detected       logical                   Indicates whether the spikes were detected using a threshold
				
                Optional:
                kst                                      KiloSort output variables 
                                                         [ see KiloSort docs for details: https://github.com/cortex-lab/KiloSort ]
				
				Nb: number of blocks
                Nc: number of channels
				
spikes.clusters                struct                   Metrics for spike clusters                  	
                metrics        struct                   
                zeta           double   [Nc x Nc]       Zeta distance between clusters in feature space (see [Ref. 7])
				
                Nc: number of clusters

spikes.clusters.metrics                                 Metrics for each spike cluster
                id             integer  scalar          Cluster identifier
				nspikes        integer  scalar          Number of spikes in cluster
				fspikes        double   scalar          Fraction of total spike number
				frate          double   scalar          Mean firing rate (Hz)
                rpv            double   scalar          Fraction of refractory period violations 
                sub            double   scalar          Fraction of sub-threshold spikes 
                co             double   scalar          Fraction of spikes that are expected to coincide with spikes from other clusters
                xc             double   vector          Cross-correlations of channel pair with greatest lag for cross-correlation peak
                xcLag          double   scalar          Greatest lag of pairwise cross-correlation peaks between channels
                cAuc           double   scalar          Area under curve of empirical and theoretically Poisson distributions
                cxDist         double   vector          Empirical Poisson distribution (x-values)
                cyDist         double   vector          Empirical Poisson distribution (y-values)
                amp            double   scalar          Mean absolute amplitude of waveform
                ampRel         double   scalar          Mean relative amplitude of waveform, normalized by threshold
				p2p            double   scalar          Mean peak-to-peak amplitude of waveform
                chans          integer  vector          Channels IDs that have above-threshold mean amplitude
                artifact       double   scalar          Ratio between actual and expected number of spikes in LFP artifact region
                snr            double   scalar          Signal-to-noise ratio
				quality        integer  scalar          Quality measure of cluster (see "psr_sst_cluster_thresholds")
                Lratio         double   scalar          L-ratio [Ref. 5]
                IsoDis         double   scalar          Isolation distance [Ref. 5]
                FP_t           double   scalar          False positive rate, based on fitting mixture of drifting t-distributions [Ref. 6]
                FN_t           double   scalar          False negative rate, based on fitting mixture of drifting t-distributions [Ref. 6]
                FP_g           double   scalar          False positive rate, based on fitting mixture of drifting Gaussians 
                FN_g           double   scalar          False negative rate, based on fitting mixture of drifting Gaussians
```       

### Temporary files

While running the processing pipeline, a number of temporary MAT files will be created in `savePath`. 
We do this to avoid having to keep the data stored in memory, allowing us to only load what we need at any one time.
These temporary MAT files are deleted after they are no longer needed. 
It is highly recommended to make sure that MATLAB deletes these files permanently, so no manual clean-up is needed. 
You can select this option by going to the `Home` tab in MATLAB, selecting `Preferences` (the cogwheel), then the `General` menu and clicking on the `Delete permanently` option under `Deleting files`. 

## Visualization and Quality Control

## Analysis

## Troubleshooting

Problem: `'Undefined function 'ft_senstype' for input arguments of type 'cell'.'`
Solution: Restart MATLAB or execute `restoredefaultpath` followed by `startup`.

## Third-party dependencies

The toolbox uses a number of third-party software packages, which are listed below. Some of these are included with the toolbox. 

```
Name:                    Website:                                                   License:
bhattacharyya            https://nl.mathworks.com/matlabcentral/fileexchange/18662  See included LICENSE file
CBPSpikesortDemo  [CBP]  https://github.com/chinasaur/CBPSpikesortDemo              GitHub Terms of Service
export_fig               https://github.com/altmany/export_fig                      BSD 3-clause "New" or "Revised" License
FieldTrip         [FT]   http://www.fieldtriptoolbox.org/                           GNU General Public License v2.0
FMMSpikeSorter    [FMM]  https://github.com/decarlson/FMMSpikeSorter                GNU General Public License v2.0
ISO-SPLIT         [ISO]  https://github.com/magland/isosplit_old                    See included COPYRIGHT file
KiloSort          [KST]  https://github.com/cortex-lab/KiloSort                     GNU General Public License v2.0
MClust (SPC only) [SPC]  http://redishlab.neuroscience.umn.edu/MClust/MClust.html   See included LICENSE file
MoDT              [MDT]  https://github.com/kqshan/MoDT                             MIT License
opass             [OPS]  https://github.com/decarlson/opass                         MIT License
Open Ephys               https://github.com/open-ephys/analysis-tools               See included LICENSE file
UMS2K             [UMS]  https://github.com/danamics/UMS2K                          See included LICENSE file
varycolor                https://nl.mathworks.com/matlabcentral/fileexchange/21050  See included LICENSE file
```

## References

[1] Pachitariu, Marius, et al. "Fast and accurate spike sorting of high-channel count probes with KiloSort." Advances in Neural Information Processing Systems. 2016.

[2] Hill, Daniel N., Samar B. Mehta, and David Kleinfeld. "Quality metrics to accompany spike sorting of extracellular signals." Journal of Neuroscience 31.24 (2011): 8699-8705.

[3] Oostenveld, Robert, et al. "FieldTrip: open source software for advanced analysis of MEG, EEG, and invasive electrophysiological data." Computational intelligence and neuroscience 2011 (2011): 1.

[4] Siegle, Joshua Handman, et al. "Open Ephys: An open-source, plugin-based platform for multichannel electrophysiology." Journal of Neural Engineering (2017).

[5] Schmitzer-Torbert, N., et al. "Quantitative measures of cluster quality for use in extracellular recordings." Neuroscience 131.1 (2005): 1-11.

[6] Shan, K. Q., Lubenov, E. V., & Siapas, A. G. (2017). Model-based spike sorting with a mixture of drifting t-distributions. bioRxiv, 109850.

[7] Yger, Pierre, et al. "Fast and accurate spike sorting in vitro and in vivo for up to thousands of electrodes." bioRxiv (2016): 067843.