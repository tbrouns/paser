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

In order to use the toolbox, a few third-party software packages need to be installed on your system.

### FieldTrip 

We use the FieldTrip toolbox [Ref. 3] for LFP processing and analysis, which can be downloaded or cloned from:

https://github.com/fieldtrip/fieldtrip

In case of problems with FieldTrip, please download the last version that was confirmed to be compatible with PASER: 

https://github.com/fieldtrip/fieldtrip/tree/af6871348043f8c912b0c9c24552f9bb8db4b412

If you are not planning on using FieldTrip for anything else, then do not add the FieldTrip toolbox to your MATLAB path. 
This will be done at a later step. Otherwise, follow the directions here:

http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path

### KiloSort

The default spike sorting method in PASER is KiloSort [Ref. 1], which can be downloaded or cloned from: 

https://github.com/cortex-lab/KiloSort

In case of problems with KiloSort, please download the last version that was confirmed to be compatible with PASER: 

https://github.com/cortex-lab/KiloSort/tree/0ea839e33527891a379e29ff9a4512d89f27bf60

It is highly recommended to use KiloSort with a CUDA enabled GPU. Attempting to run KiloSort on the CPU is errorprone and not guaranteed to result in satisfactory cluster quality. 
Therefore, please follow the installation instructions given in the KiloSort README file and "Docs" folder. 

Once again, if you are not using KiloSort for anything else, then do not add the KiloSort directory to the MATLAB path. 
We will instead load the toolbox the moment it is needed in the data processing pipeline. 

### OpenEphys

The OpenEphys toolbox is used to load the raw data:

https://github.com/open-ephys/analysis-tools

In case of problems with OpenEphys, please download the last version that was confirmed to be compatible with PASER: 

https://github.com/open-ephys/analysis-tools/tree/f66b83f09e1896b1b5874daabadde3cff9424e9c

As with the FieldTrip and KiloSort toolboxes, it is not required to add it to the MATLAB path straight away. 

## PASER toolbox installation

Clone or download PASER and add it to your MATLAB path: 

```
addpath(genpath('C:\Path\To\paser-master'));
```

## Quick start

To verify that the toolbox is working, you can download a test data set at:

https://drive.google.com/open?id=1pelaK9NgXjJh_bOGas_ujpRpmjUcpYB9

[You should download the "loadPath" folder]

The data set is a relatively short extracellular recording by just two tetrodes. 

Then create the following script:

```
psr_parameters_general;
parameters.loadPath   = 'C:\path\to\loadPath\;  % Location of the downloaded "loadPath" folder
parameters.savePath   = 'D:\path\to\savePath\;  % Where you want to save the output data 
parameters.path.ft    = 'E:\path\to\FieldTrip'; % Path to the FieldTrip main directory (folder that contains e.g. ft_defaults.m)
parameters.path.kst   = 'F:\path\to\KiloSort';  % Path to the KiloSort  main directory (folder that contains e.g. preprocessData.m)
parameters.path.ephys = 'G:\path\to\OpenEphys'; % Path to the OpenEphys main directory (folder that contains e.g. load_open_ephys_data.m)

psr_example_batch_processing(parameters)
```

Where you have to specify the five different path variables yourself. When you run the script and everything is working correctly, PASER should begin to load data from the `loadPath` directory and generate data in the `savePath` directory.
As soon as PASER is finished processing the data, a few figures will be created showing some quality metrics of a detected neuron, which are saved in the `savePath` directory. 
Furthermore, we perform a basic analysis on the data and plot some figures in the `analysis_figures` folder in the `savePath` directory. Note that these figures are only for illustrative purposes and won't show any interesting results. 

## Further details

### Default processing script

To process data using PASER, you should create a MATLAB script containing the following lines of code:

```
parameters = [];

parameters.configPath   = 'C:\path\to\ConfigFile';  % Where the parameters are loaded from
parameters.loadPath     = 'D:\path\to\LoadFolder\'; % Where the data folders are
parameters.savePath     = 'E:\path\to\SaveFolder\'; % Where you want to save the output MAT files
parameters.subject      = 'SubjectID';   % Name of subject used in output MAT filename (no spaces)
parameters.nelectrodes  = 4;             % Number of electrodes per polytrode (e.g. tetrode: 4)
parameters.extension    = 'continuous';  % File extension of raw data files
parameters.rawpattern   = 'CH';          % Pattern to look for in data files
parameters.blockpattern = [];            % Used to differentiate between blocks within session
parameters.stimpatterns = [];            % Which session type to process
parameters.process      = 'all';         % Which specific sessions to process ('all', 'given' or 'from')
parameters.folders      = [];            % Sessions that you wish to process, if 'given' or 'from' is chosen above (string cell array)
parameters.overwrite    = false;         % Whether to overwrite data from existing processed sessions
parameters.txtfile      = [];            % Folders to process given in text file

psr_batch_processing(parameters); % Process raw data files
```

This script will load and then batch process `continuous` files in the directory given by `parameters.loadPath`. 
Processed data is saved to the folder specified by `parameters.savePath`, where a matching directory tree will be created. 

Not all of the parameters in this script are set correctly, so what follows is an explanation of how to select the right settings.

### Initial parameter definitions

#### `parameters.configPath`

You should create a `ConfigFile.m` that contains parameter settings for the various processing functions and then point `parameters.configPath` to this file. 
This script should at least contain the following lines of code:

```
psr_parameters_general;
parameters.path.ft    = 'C:\path\to\FieldTrip'; % Path to the FieldTrip repository
parameters.path.kst   = 'D:\path\to\KiloSort';  % Path to the KiloSort  repository
parameters.path.ephys = 'E:\path\to\OpenEphys'; % Path to the OpenEphys repository
```

To avoid clogging up the MATLAB path, we add the third-party toolboxes given above whenever we need them and remove them afterwards. 
In order for the program to know where to look for the toolboxes, set the path parameters in the following way:

* `parameters.path.ft`    should point to the FieldTrip main directory, which contains e.g. `ft_defaults.m`
* `parameters.path.kst`   should point to the KiloSort  main directory, which contains e.g. `preprocessData.m`
* `parameters.psth.ephys` should point to the OpenEphys main directory, which contains e.g. `load_open_ephys_data.m`

More `parameters` fields can be changed by adding more lines to the script. For example, if you do not want to process the LFP, you add:
`parameters.process.lfp = false;`

See `psr_parameters_general` for all other parameter fields. Comments next to each parameter explain its purpose. 

#### `parameters.loadPath`

`parameters.loadPath` should point to a directory tree like the one given below. This is the directory tree of the test data. 

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
Each session folder then contains one or more folders for specific experimental blocks that hold the raw data files. 
In the case of `Session_1_Condition`, we have the `Block_0M` and `Block_1M`, which contain the `continuous` files. 
To be clear, the names `Session` and `Block` are arbitrary, you can use any other name you want.

##### Note on experimental "sessions" and "blocks" 

Experimental sessions are processed separately from one-another, unless specified otherwise (see `parameters.txtfile`).
Experimental blocks, on the other hand, are processed together for each probe within an experimental session.

This is based on the assumption that the electrodes will drift between experimental sessions, but roughly remain in the same position between experimental blocks in the same session.

#### `parameters.savePath`

For `parameters.savePath` select the folder where you want to save the output MAT files. Folders are automatically created to match the `loadPath` directory tree. 

##### Example for the directory tree under `parameters.loadPath`:

The folders `Session_1` and `Session_1_Condition` will be created in `savePath`, which will contain the output MAT files for the corresponding sessions. 

#### `parameters.subject`

All data that is loaded from `parameters.loadPath` is assumed to come from the same animal. In `parameters.subject` you should specify a unique character string for identification of the subject animal. 
When processing more than one animal, you should always change the `parameters.loadPath` and `parameters.subject` fields for each animal. 

#### `parameters.nelectrodes`

The number of channels of the polytrode should be indicated by `parameters.nelectrodes`.
Each block folder should then contain a number of `continuous` files equal to the number of polytrodes multiplied by `parameters.nelectrodes`.

#### `parameters.extension`

The extension of the raw data files. Right now this should always be set to `continuous`. 

#### `parameters.rawpattern`

We use this field to decide which files should be loaded from the `parameters.loadPath` directory and to determine which raw data files go with which probe (see example below). 
The pattern specified by `parameters.rawpattern` must immediately be followed by an integer in the filename. In turn, the integer should immediately be followed by the file extension or by an underscore. 

##### Example for the directory tree under `parameters.loadPath`:

We have two types of `continuous` files (`*ADC*.continuous` and `*CH*.continuous`). We only wish to load the `100_CH*.continuous` ones. 
We can differentiate between the two types using `parameters.rawpattern`, which should be set to `parameters.rawpattern = 'CH'` in this case, because that is the common pattern between them. 
If we are using a tetrode, then the channels `{*CH1,*CH2,*CH3,*CH4}.continuous` are the channels for the first tetrode.

#### `parameters.blockpattern`

Depending on the experimental conditions, you may wish to vary some kind of experimental variable across different experimental blocks. 
If you indicate the value of the variable in the block folder name, then it will be extracted and saved by setting `parameters.blockpattern`. 
Any value between the underscore and the specified pattern is recorded. 

##### Example for the directory tree under `parameters.loadPath`:

If we set `parameters.blockpattern = 'M'`, we would get `0` and `1` for the blocks in the `Session_1_Condition` folder. 

#### `parameters.stimpatterns`

Furthermore, you can also differentiate between different session types and only process one particular type of session. 
The session type that you wish to process should be set by `parameters.stimpatterns`. 

##### Example for the directory tree under `parameters.loadPath`:

If we only want to process the "*Condition" sessions (e.g. for the directory tree under `parameters.loadPath`, we would have to specify `parameters.patterns = {'condition'}` (not case sensitive). 

#### `parameters.process`

The session folders you wish to process can be directly specified as well by setting `parameters.process` to `'given'` 
and giving the session names as a cell array in `parameters.folders`. You can also set `parameters.process` to `'from'` to process the data starting from a specific session, where the starting session is specified in `parameters.folders` as a single cell. 

##### Example for the directory tree under `parameters.loadPath`:

We could set `parameters.folders = {'Session_1_Condition'}` in order to only process the `Session_1_Condition` folder and thus ignore the `Session_1` folder. 

#### `parameters.overwrite`

Set `parameters.overwrite = false` if you only want to process sessions for which no output MAT files exist in `parameters.savePath`. Alternatively, set `parameters.overwrite = true` if you want to process every single session and overwrite any existing data.

#### `parameters.txtfile`

Multiple sessions can also be processed together by creating a TXT file and specifying on each line the sessions that should be combined. A whitespace should be left between each session name. 
This is useful when sessions are recorded immediately after one another, so we can assume that the electrodes have not drifted significantly. 

##### Example for the directory tree under `parameters.loadPath`:

We can create a TXT file containing just the following line to process the two sessions together:
```
Session_1 Session_1_Condition
```

More lines can be added for more session combinations, e.g: 
```
Session_1 Session_1_Condition
Session_2 Session_2_Condition_1 Session_2_Condition_2
```

### Output files

As mentioned under `parameters.savePath`, a MAT file will be saved for each probe to the `savePath` for the current session. These files have the following naming convention:

`PSR_%SubjectID%_%Session%_P%ProbeNumber%.mat`

Where `%SubjectID%` is specified by `parameters.subject`, `%Session%` by the name of the session folder and `%ProbeNumber%` by the probe number.

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
                blocks         int16    [ 1 x Ns]        Experimental block index of each detected spike
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
                mse            double   scalar          Mean squared-error between empirical and theoretical spike count distribution (temporal stability measure)
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

## Visualization & Analysis

### Default analysis script 

After processing the data, we would like to visualize and analyse it. Similar to the data processing section, we should create a default analysis script as follows:

```
cfg = [];
cfg.loadpath = 'C:\path\to\LoadFolder'; % Location of output files from LFP and spike processing
cfg.savepath = 'D:\path\to\SaveFolder'; % Where to save the output analysis files
cfg.subject  = 'SubjectID';             % The subject we want to analysis 

cfg.analysis.run   = true; % Whether to do the analysis
cfg.analysis.fpath = 'E\path\to\AnalysisFunction'; % Location of your custom analysis function

cfg.plot.quality = true; % Plot unit quality figures
cfg.plot.merges  = true; % Plot unit merges

psr_batch_analysis(cfg);
```

Again, we explain the various fields below.

#### `cfg.loadpath`

Directory to load the processed data from. This field should probably be the same as the `parameters.savePath` field in the data processing section. 

#### `cfg.savepath`

Directory to save your analysis output files. 

#### `cfg.subject`

Which subject's data we are going to analyse. If `cfg.subject = 'SubjectID'`, then we only load processed data from folders containing the `SubjectID` string at the start of their names. 

#### `cfg.analysis.run`

Boolean indicating whether we want to run the analysis or not, i.e. whether we call the analysis function given by `cfg.analysis.fpath`. 

#### `cfg.analysis.fpath`

The `cfg.analysis.fpath` field points to an M-file that is supplied by the user to carry out data analysis. 
The M-file must contain a function, which is called in the `psr_batch_analysis` routine at the end of your analysis script, given above. 

For an analysis function example, you can look at the `psr_example_analysis.m` file in the `examples` folder. 
You should be able to run this analysis function on your own data, by setting:
```
cfg.analysis.fpath = 'C:\path\to\examples\psr_example_analysis'; % Location of the "psr_example_analysis.m" file
```
For most experiments, however, you want to create a more advanced analysis, but the `psr_example_analysis.m` file can be used as a starting point. The `psr_example_analysis.m` file has therefore been heavily annotated to make it easier for people to make their own custom analysis function. 

#### `cfg.plot.quality`

Boolean indicating whether we want to plot unit quality figures. 

#### `cfg.plot.merges`

Boolean indicating whether we want to plot the cluster merges. 

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
