function process_folder_of_OCT_scans(directory, oct_extension)

    %% HOUSEKEEPING

        close all
    
        if nargin == 0
            directory = fullfile('.', 'data');
            oct_extension = 'img'; % TODO! update if you have mixed exts
        end
        
        curr_path = mfilename('fullpath');
        [curr_dir,filename,ext] = fileparts(curr_path);
        if isempty(curr_dir); curr_dir = pwd; end
        cd(curr_dir)
                    
        % add folders to Matlab path so that their functions are found
        addpath(fullfile('.', 'BM4D'))
        addpath(fullfile('.', 'utils'))
        addpath(fullfile('.', 'L0smoothing'))
        addpath(fullfile('.', 'matlabsoftware-2013', 'ImageQualityFeatures'))
        addpath(fullfile('.', 'matlabsoftware-2013', 'DespeckleFilters'))
        
        config = read_config();
    
    %% Get file listing
    
        s = dir(fullfile(directory, ['*.', oct_extension]));
        file_list = {s.name}';
        no_of_files = length(file_list);
        disp([' - found ', num2str(no_of_files), ' image files'])
        
        % get also the list of A scan coordinates from the text file
        formatSpec = '%s %d %d %d %d %d %d';
        file_specs = readtable(fullfile(directory, config.file_listing_txt),...
                               'Delimiter','\t', 'Format',formatSpec);
        
        
    %% Process all of the files from the folder
    
        for file = 1 : no_of_files
            if(~strcmp('AO_1mo_AO_P80999_Macular Cube 512x128_1-13-2016_14-26-40_OD_sn122513_cube_z.img' ,file_list{file}))
            %continue; % XXX 
            end
            
            % check first if the file has been already denoised, so that we
            % don't have to do it again (remove the old denoised file from
            % the folder if you want to re-denoise with other method or
            % with other parameters
            ext_full = ['.', oct_extension]; % well just with the dot
            stripped_filename = strrep(file_list{file}, ext_full, '');
            denoised_filename = fullfile(directory, [stripped_filename, config.denoised_ending]);
            
            disp(['Reading the file "', file_list{file}, '"'])
            im = importZeissIMG(directory, file_list{file}); % Read the .img file in
            
            if exist(denoised_filename, 'file') == 2
                                
                disp(['The file "', file_list{file}, '" has already been denoised, skipping the denoising part'])
                % Read the .tif file in
                Z = import_tiff_stack(denoised_filename);
                
            else                
                
                % Pre-processing Image Restoration part
                % Denoises the whole cube (takes some time)
                % e.g. "1024x512x128px" Basic estimate completed (1265.4s) 
                disp([' Denoising the file "', file_list{file}, '"'])
                output = 'double';
                [Z, sigmaEst, PSNR, SSIM] = BM4D_cubeWrapper(im, output); % Z is denoised im
               % Z = imgaussfilt( im,2);
                save_in_parfor_loop(Z, denoised_filename)                
                
                % Save to disk as .TIFF
                % TODO! if you want .nrrd, .MAT,  .hdf5, or OME-TIFF, etc.
                %bitDepth = 16; % you can save some disk space with 8-bit                                
                %export_stack_toDisk(denoised_filename, Z, bitDepth)

            end
            
       % The faster part that extracts A-scans based on your desired
            % coordinates
            try
                coords_ind = check_if_coords(file_list{file}, file_specs.filename);
            catch err                
                if strcmp(err.identifier, 'MATLAB:table:UnrecognizedVarName')
                    % for now some reason, the headers are gone and we
                    % cannot reference with the field names?
                    error(['header fields are just as var1 / var2 / var3 / etc.,', ...
                           'they would need to have actually the variable names', ...
                           'You probably have error when writing the text file with too many columns, extra tabs or something'])                    
                else
                    err
                end
            end
            
             % check the filename for which eye
            if ~isempty(strfind(lower(file_list{file}), 'od'))
                disp('OD found from filename, RIGHT EYE')
                eye = 'right';
            elseif ~isempty(strfind(lower(file_list{file}), 'os'))
                disp('OS found from filename, LEFT EYE')
                eye = 'left';
            else
                warning('OD nor OS was found from the filename so the program now does not know which eye the scan was from, using RIGHT EYE coordinates')
                eye = 'right';
            end
            
            coords_ind = find_coords_from_cube(Z, coords_ind, file_specs, eye);

            
            if ~isempty(coords_ind)
                [Z_crop, A_scan, x, z_min, z_max] = crop_the_cube(Z, coords_ind, file_specs, eye);
            else
                warning(['No custom cropping coordinates found for file = "', ...
                        file_list{file}, '", how is that possible?'])
                [Z_crop, A_scan, x, z_min, z_max] = crop_the_cube(Z, coords_ind, file_specs, eye);
            end
            
            
for z = z_min:z_max
           
    
    if strcmp(stripped_filename, 'FM_1mo_FM_P71768_Macular Cube 512x128_7-30-2014_11-23-26_OS_sn94680_cube_z') &&  z <= 59
        disp('skippping because patient blinked')
        continue
    end
            
            % Save the cropped cube, you could do further
            % denoising/smoothing at this point a lot faster than on the
            % first step with the full cube.
            bitDepth = 16; % you can save some disk space with 8-bit                                
            export_stack_toDisk(strrep(denoised_filename, '.tif', '_crop.tif'), Z_crop, bitDepth)
            
            % Denoise the A-scan (from 1D line)     
            %size(A_scan);  
            
%            A_scan_denoised = smooth(A_scan(:,z), 0.01, 'loess');            
            
            % Test to denoise just the single frame and see the diff.
            frame = double(Z(:,:,z));
            
            % Despeckle, https://doi.org/10.1016/j.cmpb.2014.01.018
            frame_denoised = DsFlsmv(frame, [3 3], 5);
            
            % L0 Smoothing            
            logOffset = 0.1; lambdaS = 1e-4;
            [frame_denoised, ~, ~, ~] = enhance_logSmooothingFilter(frame_denoised, logOffset, lambdaS);                        
            A_scan_denoised_frame = frame_denoised(:,x);
            
            % Further smooth with LOESS
            % Note! Now the smoothing parameters are hand set which is
            % quite suboptimal in the end
            A_scan_denoised_2D_1D = smooth(A_scan_denoised_frame, 0.01, 'loess');
            
            % Normalize, export, compare peak ratios
            [A_scan,  ...
                A_Scan_denoised_gauss, ... 
                A_scan_denoised_frame, A_scan_denoised_2D_1D, ...                
                 GCL_PEAKS, GCL_VALUES, RPE_PEAKS, RPE_VALUES GCL_RPE_RATIOS] = ...
                compare_A_scans(A_scan,  ...
                            A_scan_denoised_frame, A_scan_denoised_2D_1D, ...
                            directory, stripped_filename);                                            

                        % OPTIONAL VISUALIZATION
                        
            gcl_rpe_ratios_filename = sprintf("%s_%d.txt", strrep(denoised_filename, '.tif', '_gcl_rpe_ratios'), z);
            gcl_rpe_ratios_file = fopen(gcl_rpe_ratios_filename,'w')
            fprintf(gcl_rpe_ratios_file,'%12.8f\n',GCL_RPE_RATIOS);
            
            
            
            gcl_values_filename = sprintf("%s_%d.txt", strrep(denoised_filename, '.tif', '_gcl_values'), z);
            gcl_values_file = fopen(gcl_values_filename,'w')
            fprintf(gcl_values_file,'%12.8f\n',GCL_VALUES);
            
            
            
            rpe_values_filename = sprintf("%s_%d.txt", strrep(denoised_filename, '.tif', '_rpe_values'), z);
            rpe_values_file = fopen(rpe_values_filename,'w')
            fprintf(rpe_values_file,'%12.8f\n',RPE_VALUES);
            
            
            
            z_of_interest = z;
            visualize_ON = 1;
            if visualize_ON == 1 
                
                scrsz = get(0,'ScreenSize');
                fig = figure('Color', 'w', 'Name', 'OCT A-Scan');
                    set(fig, 'Position', [0.1*scrsz(3) 0.1*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)])
                
                    
                rows = 2; cols = 3;
                
                % normalize for visualization
                in = double(im(:,:,z_of_interest));
                in = in - min(in(:));
                in = in / max(in(:));
                
                denoised = Z(:,:,z_of_interest);
                denoised = denoised - min(denoised(:));
                denoised = denoised / max(denoised(:));
                
                frame_denoised = frame_denoised - min(frame_denoised(:));
                frame_denoised = frame_denoised / max(frame_denoised(:));
                
                a = subplot(rows,cols,1); imshow(in, [])
                
                
                for x_index = 1:length(RPE_PEAKS) 
                    if GCL_PEAKS(x_index) == -1 || RPE_PEAKS(x_index) == -1
                        continue
                    end
                    
                    lh = line([x(x_index) x(x_index)], [GCL_PEAKS(x_index) (GCL_PEAKS(x_index) + 15)], 'Color', 'r');
                    %lh.Color=[1,0,0,0.5];

                    lh = line([x(x_index) x(x_index)], [(RPE_PEAKS(x_index) - 5) (RPE_PEAKS(x_index) + 5)], 'Color', 'b');
                    %lh.Color=[1,0,0,0.5];

                end
                
 %               line([x x], [GCL_peak_index (GCL_peak_index + 20)], 'Color', 'b')
             %   line([x x], [(RPE_peak_index - 5) (RPE_peak_index + 5)], 'Color', 'b')

                title('Input'); colorbar
                subplot(rows,cols,2); imshow(denoised, [])
                title('BM4D Denoised'); colorbar
                subplot(rows,cols,3); 
                gaussed_denoised = imgaussfilt(denoised, 12);
                imshow(gaussed_denoised, [])
                
                
                for x_index = 1:length(RPE_PEAKS) 
                    if GCL_PEAKS(x_index) == -1 || RPE_PEAKS(x_index) == -1
                        continue
                    end                   
                    lh = line([x(x_index) x(x_index)], [GCL_PEAKS(x_index) (GCL_PEAKS(x_index) + 15)], 'Color', 'r');
                    lh.Color=[1,0,0,0.25];
                    lh = line([x(x_index) x(x_index)], [(RPE_PEAKS(x_index) - 5) (RPE_PEAKS(x_index) + 5)], 'Color', 'b');
                    lh.Color=[0,0,1 ,0.25];

                end
                
                                title('Gaussian used for peak finding'); colorbar
                
                subplot(rows,cols,4); 
                
                histogram(GCL_RPE_RATIOS,40);
                xlim([0 1.5])

                title('Distribution of GCL/RPE Ratios'); colorbar; hold on
               % line([x x], [1 size(frame_denoised,1)], 'Color', 'r')
                
                subplot(rows,cols,5); 
                bplot(GCL_RPE_RATIOS);
                title('Interquartile Range of GCL/RPE Ratios'); colorbar
                
                % A-Scan Comparison
                subplot(rows,cols,6); hold on
                y = linspace(1, length(A_scan), length(A_scan));
                
                
                
                p = plot(A_Scan_denoised_gauss);
                
                for x_index = 1:length(RPE_PEAKS)
                    p_peaks = plot(RPE_PEAKS(x_index), A_Scan_denoised_gauss(RPE_PEAKS(x_index),x_index) , '^', GCL_PEAKS(x_index), A_Scan_denoised_gauss(GCL_PEAKS(x_index),x_index), '^');
                    set(p_peaks, 'MarkerFaceColor', 'k')
                end

                
                xlim([0 length(y)])
                                
  
                hold off
                 
                saveOn = 1;
                if saveOn == 1
                    filename = strrep(file_list{file}, '.img', '')
                    filename = sprintf('%s_%d',filename, z);
                    directory
                    filename
                    
                    filename = fullfile(directory, filename);
                    filename
                    try 
                                        saveas(fig,  filename, 'png')
                    catch
                        warning('failed to save');
                        warning(filename)
                    end
                    
                    close all
                end
            end
            %break % XXX stopping from going to all zs
        end
        end   


        
        
   function coords_ind = find_coords_from_cube(Z, coords_ind, file_specs, eye)
        
        disp('   Placeholder here if you want to do automated ROI localization from image')
                
    function [A_scan,A_Scan_denoised_gauss, A_scan_denoised_frame, A_scan_denoised_2D_1D, ...
            GCL_PEAKS, GCL_VALUES,  RPE_PEAKS, RPE_VALUES, GCL_RPE_RATIOS] = ...
            compare_A_scans(A_scan,  ...
                    A_scan_denoised_frame, A_scan_denoised_2D_1D, ...
                    directory, stripped_filename)

        % normalize the A_scan frame
        A_scan_denoised_frame = A_scan_denoised_frame - min(A_scan_denoised_frame(:));
        A_scan_denoised_frame = A_scan_denoised_frame / max(A_scan_denoised_frame(:));    

        % normalize the A_scan        
        A_scan_denoised_2D_1D = A_scan_denoised_2D_1D - min(A_scan_denoised_2D_1D(:));
        A_scan_denoised_2D_1D = A_scan_denoised_2D_1D / max(A_scan_denoised_2D_1D(:));
        
        % save the A-scan
        %dlmwrite(fullfile(directory, [stripped_filename, '_Ascan_raw.txt']), A_scan)
        %dlmwrite(fullfile(directory, [stripped_filename, '_Ascan_2D_denoised.txt']), A_scan_denoised_frame)
            % TODO! Maybe save some metadata if wanted at some point?    
        
        % [peak_1, peak_2, locs_peaks, ratio] = find_intensity_peaks(A_scan_denoised_frame, A_scan_denoised_2D_1D);
        [A_Scan_denoised_gauss, GCL_PEAKS, GCL_VALUES, RPE_PEAKS, RPE_VALUES, GCL_RPE_RATIOS] = find_intensity_peaks(A_scan, A_scan_denoised_frame,stripped_filename);
        
            % TODO! You could try different combinations, and you could do
            % systematic sensitivity analysis of how these parameters
            % actually affect your final estimates of the intensity ratio
            
            
    function [A_Scan_denoised_gauss,  GCL_PEAKS, GCL_VALUES, RPE_PEAKS, RPE_VALUES, GCL_RPE_RATIOS] = find_intensity_peaks(A_scan, A_Scan_denoised,stripped_filename)
        
        % Quite dumb algorithm in the end working for your canonical OCT
        % profile for sure. Think of something more robust if this start
        % failing with pathological eyes or something
        
        % Use the denoised version for peak location, and get the values
        % from the raw non-denoised version
        
        %A_scan is 1024 tall ( first index) by 151 pixels wide ( second index)
        % for each 151 pixels wide, we want to do our calculation across
        % the 1024 pixels tall segment
        
        %we will fill up an array for peak 1 location and peak 2 location
        %separately create a function for averaging and creating the ratios       
        
        x_space = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x
        GCL_PEAKS = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x 
        GCL_VALUES = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x 

        RPE_PEAKS = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x
        RPE_VALUES = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x 
        GCL_RPE_RATIOS = linspace(1, size(A_scan,2), size(A_scan,2)); % Along x
        gauss_blur_amount = 10;
        if strcmp('BN_1Month_BN_P71161_Macular Cube 512x128_7-11-2014_16-53-42_OD_sn71731_cube_z', stripped_filename) || ...
            strcmp('FM_1mo_FM_P71768_Macular Cube 512x128_7-30-2014_11-23-26_OS_sn94680_cube_z', stripped_filename) 
                gauss_blur_amount = 14;
        end
        A_Scan_denoised_gauss = imgaussfilt(A_Scan_denoised, gauss_blur_amount);

        
        for x=1:length(x_space)
            found_left_peak = 0;
            A_scan_denoised_slice = A_Scan_denoised(:,x);
            A_scan_denoised_gauss_slice = A_Scan_denoised_gauss(:,x);
            start_point =20;
            

             for y_index = start_point:size(A_scan_denoised_slice,1) % Along y
                 my_std = std(A_scan_denoised_gauss_slice);
                 delta = A_scan_denoised_gauss_slice(y_index) - mean(A_scan_denoised_gauss_slice(start_point:y_index));
 
                if delta > my_std*1
                    found_left_peak=1;
                   break;
                end
             end
             if found_left_peak == 0                               
            %     plot(A_scan_denoised_gauss_slice);
                     GCL_PEAKS(x) = -1;
                     RPE_PEAKS(x) = -1;
                     GCL_RPE_RATIOS(x) = -1;
                     
            continue
           %      error('Could not effectively find a GCL peak');
             end
            GCL_peak_index = y_index+ 20;
            GCL_PEAKS(x) = GCL_peak_index;
            GCL_peak = mean(A_scan_denoised_gauss_slice(GCL_peak_index:GCL_peak_index+15)); % Or does it make more sense to use the original values?

            
                            
                            
                
            distance_threshold = 69;
            if strcmp('AO_Pres_AO_P80999_Macular Cube 512x128_12-18-2015_15-46-18_OD_sn121527_cube_z', stripped_filename) || ...
                    strcmp('MB_Pres_MB_P72514_Macular Cube 512x128_8-6-2014_16-8-38_OD_sn95286_cube_z', stripped_filename)  ||  ...
                    strcmp('RS_1Month_RS_P71242_Macular Cube 512x128_6-27-2014_9-25-8_OS_sn71105_cube_z', stripped_filename)  || ...
                    strcmp('VR_Pres_VR_P80744_Macular Cube 512x128_12-16-2015_14-45-48_OS_sn121425_cube_z', stripped_filename)
                
                distance_threshold = 120;
            end
             
            
            slice_after_left_peak = A_scan_denoised_gauss_slice(GCL_peak_index+distance_threshold:length(A_scan_denoised_gauss_slice))  ;          
            [pks,locs] = findpeaks(slice_after_left_peak);
            [pks,I] = sort(pks, 'descend');
            locs = locs(I);
            RPE_peak_index =  locs(1) + GCL_peak_index + distance_threshold;
            
%             offset = GCL_peak_index + 15;
            %[M,RPE_peak_index] = max(A_scan_denoised_gauss_slice(offset+1:length(A_scan_denoised_slice)));
             
           % RPE_peak_index = peak_1
           % RPE_peak_index = RPE_peak_index + offset;
            RPE_PEAKS(x) = RPE_peak_index;
            RPE_peak = mean(A_scan_denoised_slice(RPE_peak_index-5:RPE_peak_index+5));
            peak_1 = RPE_peak;
            peak_2=GCL_peak;
            ratio = GCL_peak / RPE_peak;
            locs_peaks = [GCL_peak_index RPE_peak_index ];
            GCL_VALUES(x) = GCL_peak;
            RPE_VALUES(x) = RPE_peak;
            GCL_RPE_RATIOS(x) = ratio;

        end
            
        % TODO! If you feel like, you could add some uncertainty estimation with
        % Monte Carlo sampling or something similar as if your final
        % statistical analysis is done on the ratio, then if you take this
        % ratio as the "gold standard" it may lead to problems if this is
        % rather sensitivite to the actual denoising of the A-scan
        
        
    function [Z_crop, A_scan, x, z_min, z_max] = crop_the_cube(Z, coords_ind, file_specs, eye)
        
       
        if isempty(coords_ind)
            disp('Using the default crop coordinates')
            config = read_config();
            z_min = config.crop_z_window(1);
            z_max = config.crop_z_window(2);
            left_min = config.crop_left_eye(1);
            left_max = config.crop_left_eye(2);
            right_min = config.crop_right_eye(1);
            right_max = config.crop_right_eye(2);
        else
            z_min = file_specs.z_min(coords_ind);
            z_max = file_specs.z_max(coords_ind);
            left_min = file_specs.left_min(coords_ind);
            left_max = file_specs.left_max(coords_ind);
            right_min = file_specs.right_min(coords_ind);
            right_max = file_specs.right_max(coords_ind);
        end
        
        if strcmp(eye, 'left')
            x_min = left_min;
            x_max = left_max;
            
        elseif strcmp(eye, 'right')
            x_min = right_min;
            x_max = right_max;
            
        else
            error(['Typo with your eye, it is now: "', eye, '"'])
        end
        
        % The actual crop
        Z_crop = Z(:,x_min:x_max,z_min:z_max);
        
        % One A-scan that is the one in the middle of the x and z-range
        x = round((x_min + x_max)/2);
        x = linspace(double(x_min),double(x_max),x_max - x_min+1); 
        x
        

        %loop here 
        %TODO
        
        
        z = round((z_min + z_max)/2);
        z = z_max
        A_scan = Z(:,x,z);
        
        % normalize the A_scan
        A_scan = A_scan - min(A_scan(:));
        A_scan = A_scan / max(A_scan(:));
        
        % TODO! There is possibly quite a lot of empty space left still on
        % y-axis direction (note! in Matlab the first dimension is y)
        
     
        
    function coords_ind = check_if_coords(filename, list_of_filenames_with_coords)
        
        IndexC = strfind(list_of_filenames_with_coords, filename);
        coords_ind = find(not(cellfun('isempty', IndexC)));
        
    function save_in_parfor_loop(Z, filename)
        
        disp(['  Saving as double-precision .mat: "', strrep(filename, '.tif', '.mat'), '"'])
        save(strrep(filename, '.tif', '.mat'), 'Z')
        
    
    
