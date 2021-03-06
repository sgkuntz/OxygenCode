clearvars *

% the directory with the time series of interest is identified
timelapsefile = fopen('/Users/eisenlab/Desktop/Imaging/TimeLapseInstructions_notRun_trial15.txt','r');
timelapseinfo = textscan(timelapsefile, '%s %s %s %s %s %s %s %s %s %s %s');
fclose(timelapsefile);

for w = 1:length(timelapseinfo{1})
	startingDate = timelapseinfo{1}{w}
	startingPosition = timelapseinfo{2}{w}
	starterImage = timelapseinfo{3}{w};
	fileNameScheme = timelapseinfo{5}{w};
	zPosition = timelapseinfo{11}{w}		
	% need to call up the aligned images only, not the raw images
	% cd /Users/eisenlab/Desktop/Imaging
	cd /Volumes/Coluber/Imaging/
	cd(startingDate)
	cd(startingPosition)
	% cd aligned
	p = which(starterImage);
	filelist = dir([fileparts(p) filesep fileNameScheme]);
	fileNames = {filelist.name}'; nImages = length(fileNames);

	%%% STAGE 1B - DETERMINE OPTIMAL FOCUS
	% only need to go through the images, not each event
	FscoresList = []; F_normedVarList = []; F_autocorrList = []; F_stdevcorrList = [];
	height = 705; width = 250;

	try
		
		% cycle through the images
		for k = 1:nImages
			I = imread(fileNames{k}); 

			% Normalized Variance
			runningSum = 0.0; normalizedSum = 0.0;
			for h = 1:height
				for w = 1:width
					runningSum = runningSum + double(I(h,w));
				end
			end
			mean = runningSum/(height*width);

			for h = 1:height
				for w = 1:width
					normalizedSum = normalizedSum + (double(I(h,w)) - mean)^2;
				end
			end
			F_normedVar = 1/(height*width*mean)*normalizedSum;

			% Auto-Correlation
			product1 = 0.0; product2 = 0.0;
			for h = 1:height
				for w = 1:width-1
					product1 = product1 + double(I(h,w)) * double(I(h,w+1));
				end
				for w = 1:width-2
					product2 = product2 + double(I(h,w)) * double(I(h,w+2));
				end
			end
			F_autocorr = (product1-product2)/1000000;

			% Standard Deviation Correlation

			% use the mean from Normalized Variance
			F_stdevcorr = 0.0;
			for h = 1:height
				for w = 1:width-1
					F_stdevcorr = F_stdevcorr + (double(I(h,w)) * double(I(h,w+1)));
				end
			end

			F_stdevcorr = (F_stdevcorr - (height*width*mean^2))/100000000;
			% record the 3 F scores for each image
			Fscores = horzcat(k, F_normedVar, F_autocorr, F_stdevcorr);
			F_normedVarList = vertcat(F_normedVarList, F_normedVar);
			F_autocorrList = vertcat(F_autocorrList, F_autocorr);
			F_stdevcorrList = vertcat(F_stdevcorrList, F_stdevcorr);
			FscoresList = vertcat(FscoresList, Fscores);

		end
		% can do the smoothing here
		F_normedVar_SGF = sgolayfilt(F_normedVarList, 4, 41);
		F_autocorr_SGF = sgolayfilt(F_autocorrList, 4, 41);
		F_stdevcorr_SGF = sgolayfilt(F_stdevcorrList, 4, 41);
		FscoresList_SGF = [];
		for (l = 1:nImages)
			Fscores_SGF = horzcat(F_normedVar_SGF(l), F_autocorr_SGF(l), F_stdevcorr_SGF(l));
			FscoresList_SGF = vertcat(FscoresList_SGF, Fscores_SGF);
		end
		cd /Users/eisenlab/Desktop/Imaging/
		savefile_st = horzcat('FocusScoresUnaligned_', char(startingDate), '_', char(startingPosition),'_', char(zPosition),'.txt'); dlmwrite(savefile_st, FscoresList_SGF, 'delimiter', '\t', 'precision','%.6f');

	catch
		FscoresList = []; F_normedVarList = []; F_autocorrList = []; F_stdevcorrList = [];
		height = 250; width = 705;
		
		% cycle through the images
		for k = 1:nImages
			I = imread(fileNames{k}); 

			% Normalized Variance
			runningSum = 0.0; normalizedSum = 0.0;
			for h = 1:height
				for w = 1:width
					runningSum = runningSum + double(I(h,w));
				end
			end
			mean = runningSum/(height*width);

			for h = 1:height
				for w = 1:width
					normalizedSum = normalizedSum + (double(I(h,w)) - mean)^2;
				end
			end
			F_normedVar = 1/(height*width*mean)*normalizedSum;

			% Auto-Correlation
			product1 = 0.0; product2 = 0.0;
			for h = 1:height
				for w = 1:width-1
					product1 = product1 + double(I(h,w)) * double(I(h,w+1));
				end
				for w = 1:width-2
					product2 = product2 + double(I(h,w)) * double(I(h,w+2));
				end
			end
			F_autocorr = (product1-product2)/1000000;

			% Standard Deviation Correlation

			% use the mean from Normalized Variance
			F_stdevcorr = 0.0;
			for h = 1:height
				for w = 1:width-1
					F_stdevcorr = F_stdevcorr + (double(I(h,w)) * double(I(h,w+1)));
				end
			end

			F_stdevcorr = (F_stdevcorr - (height*width*mean^2))/100000000;
			% record the 3 F scores for each image
			Fscores = horzcat(k, F_normedVar, F_autocorr, F_stdevcorr);
			F_normedVarList = vertcat(F_normedVarList, F_normedVar);
			F_autocorrList = vertcat(F_autocorrList, F_autocorr);
			F_stdevcorrList = vertcat(F_stdevcorrList, F_stdevcorr);
			FscoresList = vertcat(FscoresList, Fscores);

		end
		% can do the smoothing here
		F_normedVar_SGF = sgolayfilt(F_normedVarList, 4, 41);
		F_autocorr_SGF = sgolayfilt(F_autocorrList, 4, 41);
		F_stdevcorr_SGF = sgolayfilt(F_stdevcorrList, 4, 41);
		FscoresList_SGF = [];
		for l = 1:nImages
			Fscores_SGF = horzcat(F_normedVar_SGF(l), F_autocorr_SGF(l), F_stdevcorr_SGF(l));
			FscoresList_SGF = vertcat(FscoresList_SGF, Fscores_SGF);
		end
		cd /Users/eisenlab/Desktop/Imaging/
		savefile_st = horzcat('FocusScoresUnaligned_', char(startingDate), '_', char(startingPosition),'_', char(zPosition),'.txt'); dlmwrite(savefile_st, FscoresList_SGF, 'delimiter', '\t', 'precision','%.6f');
		% failedRun = startingDate
		% failedRun = startingPosition
	end
end
