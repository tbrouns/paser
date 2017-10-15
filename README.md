# PASER

MATLAB toolbox for processing and analyzing extracellular recordings. Intended for local field potential (LFP) and spike data processing, analysis and visualization.  

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
parameters.configPath   = 'E:\PathToConfigFile\ConfigFile.m; % Where the parameters are loaded from
parameters.patterns     = [];            % Used to differentiate between experimental sessions (string cell array)
parameters.type         = 'all';         % Which session type to process ('all' or one of the chosen patterns)
parameters.process      = 'new';         % Which specific sessions to process ('new', 'given' or 'all')
parameters.folders      = [];            % Sessions that you wish to process, if 'given' is chosen (string cell array)
parameters.extension    = 'continuous';  % File extension of raw data files
parameters.filepattern  = 'CH';          % Pattern to look for in data files
parameters.trialpattern = [];            % Used to differentiate between trials within session
parameters.nelectrodes  = 4;             % Number of electrodes per polytrode (e.g. tetrode: 4)

parameters.path.kst = 'C:\PathToKiloSort\KiloSort';   % Path to the KiloSort repository
parameters.path.ft  = 'D:\PathToFieldTrip\fieldtrip'; % Path to the FieldTrip repository

psr_batch_processing(parameters); % Process raw data files
```

This script will load and then batch process `continuous` files in the directory given by `parameters.loadPath`. 
Processed data is saved to the folder specified by `parameters.savePath`, where a matching directory tree will be created. 

Not all of the parameters in this script are set correctly, so what follows is an explanation of how to select the right settings.

### Initial parameters

#### Load directory

`parameters.loadPath` should point to a directory tree like the one given below. 

```
loadPath
│
├───Session_1
│   │
│   ├───Trial_0M
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
│   ├───Trial_0M
│   │   │   100_ADC1.continuous
│   │   │   100_ADC2.continuous
│   │   │   ...
│   │   │
│   │   │   100_CH1.continuous
│   │   │   100_CH2.continuous
│   │   │   ...
│   │
│   ├───Trial_1M
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
Each session folder then contains one or more folders for specific trials that hold the raw data files. 
In the case of `Session_1_Condition`, we have the `Trial_0M` and `Trial_1M`, which contains the `continuous` files. 
To be clear, the names `Session` and `Trial` are arbitrary here, you can use any other name you desire.

However, we must note that we have two types of `continuous` files here. We only wish to load the `100_CH*.continuous` ones. 
We can differentiate between the two types using `parameters.filepattern`, which should be set to `parameters.filepattern = 'CH'` in this case, because that is the common pattern between them. 

Depending on the experimental conditions, you may wish to vary some kind of experimental variable across different trials. 
If you indicate the value of the variable in the trial folder name, then it will be extracted and saved by setting `parameters.trialpattern`. 
Any value between the underscore and the specified pattern is recorded. If we set `parameters.trialpattern = 'M'`, we would get 0 and 1 for the trials in `Session_1_Condition`. 

Furthermore, you can also differentiate between different session types and only process one particular type of session. 
First set the different session types in `parameters.patterns` and in `parameters.type` which type you wish to process. 
If we only want to process the 'Condition' sessions, we would have to specify `parameters.patterns = {'condition'}` and `parameters.type = 'condition'` (not case sensitive). 

Lastly, the session folders you wish to process can be directly specified as well by setting `parameters.process` to `'given'` 
and giving the session names as a cell array in `parameters.folders` (e.g. `parameters.folders = {'Session_1_Condition'}`).
Otherwise, you can set `parameters.process` to `'new'` when you only want to process sessions for which no output MAT files exist in the `savePath` or `'all'` if you want to do every single session.

#### Note on trials

We assume that trials occur one after another during the experiment, which means that the electrodes will remain in the same position. 
Therefore, we always perform spike sorting across different trials for each probe within a session.

#### Polytrode channels

The number of channels of the polytrode should be indicated by `parameters.nelectrodes`.
Each trial folder should then contain a number of `continuous` files equal to the number of polytrodes multiplied by `parameters.nelectrodes`.
 
#### Save directory

For `parameters.savePath` select the folder where you want to save the output MAT files. Folders are automatically created to match the `loadPath` directory tree. 
So, in the example given above, the folders `Session_1` and `Session_1_Condition` will be created in `savePath`, which will contain the output MAT files for the corresponding sessions. 

#### Config file

You should create a `ConfigFile.m` that contains parameter settings for the various processing functions and then point `parameters.configPath` to this file. 
The file `psr_parameter_default.m` gives typical values for the parameters, which you can copy to create `ConfigFile.m` and make changes as you see fit. 
Comments next to each parameter explain its purpose. 

#### Path settings

To avoid clogging up the MATLAB path, we add the third-party toolboxes from the path whenever we need them and remove them afterwards. 
In order for the program to know where to look for the toolboxes, set the path parameters in the following way:

* `parameters.path.kst` should point to the KiloSort main directory, which contains e.g. `preprocessData.m`
* `parameters.path.ft` should point to the FieldTrip main directory, which contains e.g. `ft_defaults.m`

### Output files

As mentioned earlier, a MAT file will be saved for each probe to the `savePath` for the current session. These files have the following naming convention:

`Spikes_%SubjectName%_%Session%_P%ProbeNumber%_%Method%.mat`

Where `%SubjectName%` is specified by `parameters.subject`, `%Session%` by the name of the session folder, 
`%ProbeNumber%` by the probe number and `%Method%` by `parameters.sorting.method` (see `psr_parameter_default.m`).

The data in every MAT file is organized in the following way:

```
Structures:     Fields:	       Type:   Size:            Description:

freq                           struct                   Local field potential (LFP) data                         
                artifacts      double  [Na x 2]         Start and end point of LFP artifacts [sec]
                cfg            struct                   Configuration parameters
                dimord         string                   Defines how the numeric data should be interpreted
                freq           double  [1 x Nf]         Frequency bins [Hz]
                label          cell    [1 x Nc]         Cell array of label for each channel (see next line)
                               string                   Channel label
                powspctrm      double  [Nc x Nf x Nt]   Power per channel, frequency bin and time bin
                time           double  [1 x Nt]         Time bins [sec]
				
                [ see FieldTrip docs for details: http://www.fieldtriptoolbox.org/ ]

                Na: number of artifacts	
                Nc: number of channels
                Nf: number of frequency bins
                Nt: number of time bins
				
metadata                       struct                   General experimental data
                subject        string                   Name of subject
                session        string                   Name of session
                dir            string                   Location of raw data files
                stimulus       double  [1 x Nt]         Vector of trial conditions
                probe          integer [1 x Np]         Probe number
				
                Np: number of probes
                Nt: number of trials
				
parameters                     struct                   Parameters for all data processing functions
                
                [ see "psr_parameter_default" for details ]
			
spikes                         struct                   Neural spiking data
                assigns        integer  [1 x Ns]        Cluster index of each detected spike after merging
                assigns_prior  integer  [1 x Ns]        Cluster index of each detected spike before merging
                clusters       struct                   See "psr_sst_clusterfeatures" for details
                info           struct                   See further below
                params         struct                   Needed for compatibility
                removed        logical  [1 x Ns]        Spikes that are designated for removal
                spiketimes     double   [1 x Ns]        Spike time of each detected spike [sec]
                trials         integer  [1 x Ns]        Trial index of each detected spike
                waveforms      double   [Ns x Np x Nc]  Waveform of each detected spike for each channel 
				
                Ns: number of spikes
                Np: number of data points
                Nc: number of channels
				
spikes.info                    struct                   Session information                                         
                detect.dur     double   scalar          Length of session [sec]
                detect.thresh  double   [1 x Nc]        Spike detection thresholds for each channel
                stimtimes      cell     [1 x Nt]        Contains arrays of stimulus onset times for each trial (see next line)
                               double   [1 x Ns]        Stimulus onset times [sec]
                trialonset     double   [1 x Nt]        Trial onset times [sec]
				
                Optional:
                kst                                     KiloSort output variables 
                                                        [ see KiloSort docs for details: https://github.com/cortex-lab/KiloSort ]
				
                Nc: number of channels
                Ns: number of stimuli
                Nt: number of trials
```

## Visualization and Quality Control

## Analysis

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

[1] Pachitariu, Marius, et al. "Kilosort: realtime spike-sorting for extracellular electrophysiology with hundreds of channels." BioRxiv (2016): 061481.

[2] Hill, Daniel N., Samar B. Mehta, and David Kleinfeld. "Quality metrics to accompany spike sorting of extracellular signals." Journal of Neuroscience 31.24 (2011): 8699-8705.

[3] Oostenveld, Robert, et al. "FieldTrip: open source software for advanced analysis of MEG, EEG, and invasive electrophysiological data." Computational intelligence and neuroscience 2011 (2011): 1.

[4] Siegle, Joshua Handman, et al. "Open Ephys: An open-source, plugin-based platform for multichannel electrophysiology." Journal of Neural Engineering (2017).