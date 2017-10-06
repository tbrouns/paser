function ops = ept_kst_convert2binary(ops)

nTrials = size(ops.files,2);

clear fs
fidout = fopen(ops.fbinary, 'w');

for iBlock = 1:nTrials
    for iChan = 1:ops.Nchan
        fs{iChan,iBlock} = ops.files(iChan);
    end
end

nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks)) > 1
    error('different number of blocks for different channels!')
end
%
nBlocks  = unique(nblocks);
nSamples = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, nTrials);

data_int16 = []; % TEMP

tic
for iTrial = 1:nTrials
    for iBlock = 1:nBlocks
        for iChan = 1:ops.Nchan
            fid{iChan,iTrial} = fopen(fullfile(ops.root{iTrial}, char(fs{iChan,iTrial}(iBlock))));
            % discard header information
            fseek(fid{iChan,iTrial}, 1024, 0);
        end
        nsamps = 0;
        flag = 1;
        while 1
            samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
            for iChan = 1:ops.Nchan
                collectSamps = zeros(nSamples * 1000, 1, 'int16');
                rawData      = fread(fid{iChan,iTrial},  1000 * (nSamples + 6), '1030*int16', 10, 'b');
%                 rawData      = fread(fid{iChan,iTrial},  1000 * (nSamples + 6), '1030*int16');
                nbatches     = ceil(numel(rawData) / (nSamples + 6));
                for iBatch = 1:nbatches
                    rawSamps = rawData((iBatch - 1) * (nSamples + 6) + 6 + (1:nSamples));
                    collectSamps((iBatch - 1) * nSamples + (1:nSamples)) = rawSamps;
                end
                samples(:,iChan) = collectSamps;
            end
            
            if nbatches < 1000
                flag = 0;
            end
            if flag == 0
                samples = samples(1:iBatch*nSamples, :);
            end
            
            samples = samples';
            fwrite(fidout, samples, 'int16');
            
            
            nsamps = nsamps + size(samples,2);
            
            if flag == 0
                break;
            end
            %%%%%%%%%%%%%% TEMP
            data_int16 = [data_int16,samples]; % TEMP
            %%%%%%%%%%%%%% TEMP
            
        end
        
        %%%%%%%%%%%%%% TEMP
        for iChan = 1:ops.Nchan
            file = fullfile(ops.root{iTrial}, char(fs{iChan,iTrial}(iBlock))); % TEMP
            [data, ~, ~] = load_open_ephys_data(file);  % TEMP
            if (iChan == 1); data_double = data'; % TEMP
            else,            data_double = [data_double;data']; % TEMP
            end % TEMP
        end
        %%%%%%%%%%%%%% TEMP
        
        ops.nSamplesBlocks(iBlock) = nsamps;
        
        for iChan = 1:ops.Nchan
            fclose(fid{iChan,iTrial});
        end
        
    end
end

fclose(fidout);

toc