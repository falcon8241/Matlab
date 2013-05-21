%% Parser of the CAMERA output results to assist analysis of ion
%% annotation.
%% Note: the monoisotpic mass means the mass of [M+H] type of ion with
%% the mono isotopes.
%% This version works with XLSX files.
%% Author: Bin Zhou
%% Date: 09-21-10
clear all; close all; clc;
%% Input
annotation_filename = 'C:\Users\bz37\Desktop\test.xls';
polarity = 'neg';
%% Read the annotation file
[num, txt, raw] =xlsread(annotation_filename);
PROTO_MASS = 1.00727646677;
%% Read isotope, adduct, pcgroup and m/z information
raw_header = raw(1, :);
raw = raw(2:end, :);
isotopes=raw(:, end-2);
adducts=raw(:, end-1);
mz = cell2mat(raw(:, 2));
pcgroup=cell2mat(raw(:, end));
for i = 1: length(adducts),
    if(~all(isnan(adducts{i})))
        adduct1 = adducts{i};
        idx = find(strcmp(adduct1, adducts)&(pcgroup(i)==pcgroup));
        if(length(idx)>1)
            disp_txt = sprintf('Find same annotation in the same group, both of them are removed\n%s, group %d\n', adducts{i}, pcgroup(i));
            disp(disp_txt)
            adducts(idx)={NaN};
        elseif(isempty(idx))
            error('Fatal Error')
        end
    end
end
raw_cleaned = raw;
raw_cleaned(:, end-1) = adducts;
raw = raw_cleaned;
raw_cleaned = [raw_header; raw_cleaned];
for j = 1: length(isotopes)
    if(isnan(isotopes{j}))
        isotopes{j}='';
    end
    if(isnan(adducts{j}))
        adducts{j}='';
    end
end
pcgroup_no = sort(unique(pcgroup));
%% Find and split those ions with more than one annotations
raw_single_annotation = [];
ambiguity_flg = [];
for j = 1: length(adducts),
    if(~isempty(regexp(adducts{j}, '\[.+\].+\[.+\]', 'match')))
        disp_txt = sprintf('Following conflict needs to be resolved in pc_group %d\n', pcgroup(j));
        disp(disp_txt)
        single_annotations = regexp(adducts{j}, '(?<=.\d+) ', 'split');
        % If the multiple adduct annotation coincide with isotope
        % annotation, the isotope annotation only preserve for the first
        % adduct annotation. Otherwise will cause error. One possible
        % solution in the future is to add new isotope group and split
        % related isotope annotations.
        for k = 1: length(single_annotations);
            single_line = raw(j, :);
            single_line(1, end - 1) = single_annotations(k);
            if(k > 1)
                single_line(1, end-2) = {NaN};
            end
            raw_single_annotation = [raw_single_annotation; single_line];
            ambiguity_flg = [ambiguity_flg; 1];
        end
    else
        raw_single_annotation = [raw_single_annotation; raw(j, :)];
        ambiguity_flg = [ambiguity_flg; 0];
    end
end
%% Re-extract the annotations
isotopes=raw_single_annotation(:, end-2);
adducts=raw_single_annotation(:, end-1);
for j = 1: length(isotopes)
    if(isnan(isotopes{j}))
        isotopes{j}='';
    end
    if(isnan(adducts{j}))
        adducts{j}='';
    end
end
mz = cell2mat(raw_single_annotation(:, 2));
pcgroup=cell2mat(raw_single_annotation(:, end));
pcgroup_no = sort(unique(pcgroup));
%% Prepare the flags and monoisotopic m/zs. If an ion is not processed, the
%% flag remains NaN.
ion_info = zeros(1, 3);
adduct_flg=nan.*ones(length(pcgroup),1);
iso_flg=nan.*ones(length(pcgroup),1);
mono_flg=nan.*ones(length(pcgroup),1);
ion_grp=nan.*ones(length(pcgroup),1);
iso_grp=nan.*ones(length(pcgroup),1);
mono_mass=nan.*ones(length(pcgroup),1);
for i = 1: length(pcgroup_no),
    curr_group = pcgroup_no(i);
    ion_idx = find(pcgroup==curr_group);
    group_iso = isotopes(ion_idx, :);
    group_add = adducts(ion_idx, :);
    group_mz = mz(pcgroup==curr_group);
%     group_iso=regexprep((group_iso),'"', '');
%     group_add=regexprep((group_add),'"', '');
%% Parse the adducts for a pcgroup. If it is [M+H], it is a monoisotopic
%% ion. If it is annotated otherwise, it is an adduct ion. If there is no
%% annotation, it is NOT an adduct ion at least, but not necessarily a
%% monoisotopic ion. It should be noted if an ion is annotated as an
%% adduct, only the one of mono isotopes is annotated. Although the adduct
%% ion can have its set of isotopic ions. That's the reason why adducts are
%% first processed here.
    for j = 1: length(group_add)
        if(~isempty(group_add{j}))
            rule = regexp(group_add{j}, '\[.+\]', 'match');
            %mass_tmp is the equivalent m/z of [M+H] or [M-H] type of ion
            if(strcmpi(polarity, 'pos'))
                mass_tmp = str2double(cell2mat(regexp(group_add{j}, '\d+.\d+', 'match'))) + PROTO_MASS;
            elseif(strcmpi(polarity, 'neg'))
                mass_tmp = str2double(cell2mat(regexp(group_add{j}, '\d+.\d+', 'match'))) - PROTO_MASS;
            else
                error('Wrong Polarity Setting')
            end
            if(strcmp(rule, '[M+H]')==true)
                mono_flg(ion_idx(j)) = true;
                mono_mass(ion_idx(j)) = group_mz(j);
                adduct_flg(ion_idx(j)) = false;
            elseif(strcmp(rule, '[M-H]')==true)
                mono_flg(ion_idx(j)) = true;
                mono_mass(ion_idx(j)) = group_mz(j);
                adduct_flg(ion_idx(j)) = false;
            else
                adduct_flg(ion_idx(j)) = true;
                mono_flg(ion_idx(j)) = false;
            end
            %% If the ion has the same monoisotopic value and pcgroup as a
            %% previos ion, they are assigned to the same ion group with
            %% same monoisotopic mass. If not, a new ion_info entry is
            %% created
            if(any((mass_tmp==ion_info(:, 2))&(curr_group==ion_info(:, 3))))
                ion_grp_idx = find((mass_tmp==ion_info(:, 2))&(curr_group==ion_info(:, 3)));
                ion_grp(ion_idx(j)) = ion_info(ion_grp_idx, 1);
            else
                ion_grp(ion_idx(j)) = max(ion_info(:, 1)) + 1;
                ion_info_temp = [max(ion_info(:, 1))+1, mass_tmp, curr_group];
                ion_info = [ion_info; ion_info_temp];
            end
        else
            adduct_flg(ion_idx(j)) = false;
        end
    end
    %% Parse the isotopes of a pcgroup. If an ion is [M]+ and not an adduct of
    %% the other ion, it is a monoisotopic ion, otherwise it is an adduct ion
    %% with monoisotopic mass (but not considered as monoisotopic ion here). If
    %% an ion is multiple-charged [M]n+(it is labeled under "isotopes" rather than
    %% "adduct"), its monoisotopic mass with single charge is calculated (as
    %% n*m/z - (n-1)*mass-of-Proton). If an ion is labeled but neither [M]+ nor
    %% [M]n+, it is an isotopic ion. If an ion is not labeled, it is not an
    %% isotopic ion.
    for j = 1: length(group_iso)
        if(~isempty(group_iso{j}))
            iso_grp(ion_idx(j)) = str2double(cell2mat(regexp(group_iso{j}, '(?<=\[)\d+(?=\])', 'match')));
            rule = regexp(group_iso{j}, '\[M(\+\d)?\](\d)?\+', 'match');
            if(strcmp(rule, '[M]+'))
                if(adduct_flg(ion_idx(j)) == false)
                    mono_flg(ion_idx(j)) = true;
                    mono_mass(ion_idx(j)) = group_mz(j);
                    iso_flg(ion_idx(j)) = false;
                else
                    mono_flg (ion_idx(j)) = false;
                    iso_flg(ion_idx(j)) = false;                    
                end
            elseif(~isempty(str2double(regexp(rule{1}, '(?<=\[M\])\d+(?=\+)', 'match'))))
                charge_state = str2double(regexp(rule{1}, '(?<=\[M\])\d+(?=\+)', 'match'));
                mono_flg(ion_idx(j)) = false;
                if(strcmpi(polarity, 'pos'))
                    mono_mass(ion_idx(j)) = charge_state*group_mz(j) - (charge_state-1)*PROTO_MASS;
                elseif(strcmpi(polarity, 'neg'))
                    mono_mass(ion_idx(j)) = charge_state*group_mz(j) + (charge_state-1)*PROTO_MASS;
                else
                    error('Wrong Polarity Setting')
                end
                iso_flg(ion_idx(j)) = false;
            else
                iso_flg(ion_idx(j)) = true;
                mono_flg(ion_idx(j)) = false;
            end
        else
            iso_flg(ion_idx(j)) = false;
        end
    end
end
%% For ions which are neither adducts nor isotopes and has not been assigned
%% a monoisotopic m/z yet, they are assumed to be monoisotopic and the
%% monoisotopic m/z are assigned.
for i = 1: length(pcgroup),
    if((~adduct_flg(i))&&(~iso_flg(i))&&(isnan(mono_mass(i))))
        mono_flg(i) = true;
        mono_mass(i) = mz(i);
    end
end
ion_info = ion_info(2: end, :);
%% If an ion group has only one member, there is a mistake.
for k = 1: size(ion_info, 1),
    if(sum(ion_grp == k)==1)
%         ion_grp(ion_grp==k) = nan;
        disp('Check following entry')
        disp(raw_single_annotation(ion_grp==k,end-2:end));
        error('Ion group %d has only one ion member', k);
    end
end
%% If multiple ions are in the same ion group, they are adducts for each
%% other and share the same monoisotopic mass. If one of them is previously
%% assigned a mono isotopic mass, use the assigned one. If none of them has
%% monoisotopic mass (e.g. they are [M+Na] and [M+K]), use the one
%% calculated by CAMERA and stored in ion_info.
for k = 1: size(ion_info, 1),
    grp_ion_idx = find(ion_grp == ion_info(k, 1));
    grp_mass = mono_mass(grp_ion_idx);    
    if(sum(~isnan(grp_mass))>0)
        if(sum(~isnan(grp_mass))==1)
            mono_mass(grp_ion_idx) = grp_mass(~isnan(grp_mass));
        elseif((sum(~isnan(grp_mass)) > 1)&&(sum(mono_flg(grp_ion_idx))==1))
            % If there are both [M+2H]2+ and [M]2+ type of annotation, that
            % means there are both adducts and isotopes of this ion. In
            % this case, there could be two valid monoisotopic mass from
            % previous calculation. One from [M+2H] and the other from
            % [M+H]([M+Na] will not give a mono mass). In this case, we use
            % the monoisotopic m/z from [M+H]+
            mono_mass(grp_ion_idx) = grp_mass((~isnan(grp_mass))&mono_flg(grp_ion_idx));
        else
            disp(raw_single_annotation(grp_ion_idx,end-2:end))
            error('Error: more than one mono mass value for adducts')
        end
    else
        mono_mass(grp_ion_idx) = ion_info(k, 2);
    end
end
%% If multiple ions are in the same isotope group, they are isotopes for
%% each other and have the same monoisotopic mass. For one isotopic ions to
%% appear, its monoisotopic ion must present.
iso_grp_no = unique(iso_grp(~isnan(iso_grp)));
for k = 1: length(iso_grp_no),
    iso_grp_idx = find(iso_grp == iso_grp_no(k));
    iso_grp_mass = mono_mass(iso_grp_idx);
    if(sum(~isnan(iso_grp_mass)) > 0)
        if(sum(~isnan(iso_grp_mass))==1)
            mono_mass(iso_grp_idx) = iso_grp_mass(~isnan(iso_grp_mass));
        else
            error('Error: more than one mono mass value for isotopes')
        end
    else
        error('Error: there is no monoisotopic ion found')
    end
end
%% Group the ions into metabolites
label_flg = zeros(1, length(pcgroup));
metabolite_grp = zeros(1, length(pcgroup));
metabolite_idx = 1;
%% The adducts are first grouped into metabolites
for k = 1: size(ion_info, 1),
    ion_grp_idx = find(ion_grp == ion_info(k, 1));
    if(~isempty(ion_grp_idx)),
        label_flg(ion_grp_idx) = true;
        metabolite_grp(ion_grp_idx) = metabolite_idx;
        metabolite_idx = metabolite_idx+1;
    else
        error('Error: ion group does not exist for ion group No. %d', ion_info(k, 1));
    end
end
%% The isotopes are then grouped into metabolites. If a metabolite number
%% already exists from the adducts of the isotopics ions, use the existing
%% one. Otherwise, use a new one.
for k = 1: length(iso_grp_no),
    iso_grp_idx = find(iso_grp == iso_grp_no(k));
    curr_metabolite_grp_no = metabolite_grp(iso_grp_idx);
    if(any(curr_metabolite_grp_no > 0))
        if(sum(curr_metabolite_grp_no > 0)==1)
            metabolite_grp(iso_grp_idx) = curr_metabolite_grp_no(curr_metabolite_grp_no > 0);
            label_flg(iso_grp_idx) = true;
        else
            error('Error: more than one metabolite group number for isotopes')
        end
    else
        metabolite_grp(iso_grp_idx) = metabolite_idx;
        label_flg(iso_grp_idx) = true;
        metabolite_idx = metabolite_idx + 1;
    end
end
%% For the remaining ions, assign a metabolite number
metabolite_grp(~label_flg) = metabolite_idx: metabolite_idx + sum(~label_flg)-1;
%% Output the parsed results into a .xlsx file in the same folder
Annotation_content = [mono_mass, metabolite_grp.', mono_flg, adduct_flg, iso_flg];
Annotation_header = {'mono_mz' ,'metabolite_group', 'monoisotpic_flg', 'adduct_flg', 'isotope_flg'};
Annotation_output = [[raw_header; raw_single_annotation], ...
    [Annotation_header;num2cell(Annotation_content)], ['ambiguity_flg'; num2cell(ambiguity_flg)]];
[folder, filename] = fileparts(annotation_filename);
output_filename = sprintf('%s_parsed.xlsx', fullfile(folder, filename));
xlswrite(output_filename, Annotation_output)

uni_metabolite = unique(metabolite_grp);
mono_selector = zeros(1, length(mono_flg));
for i=1: length(uni_metabolite),
    idx = find(metabolite_grp==uni_metabolite(i));
    idx_mono = find(mono_flg(idx));
    if(length(idx_mono)==1),
        mono_selector(idx(idx_mono))=1; %If there is monoisotopic (de)protonated ion, use that one gives a higher accuracy
    elseif(length(idx_mono)>=2)
        error('More than one monoisotopic peak')
    else
        mono_selector(idx(1))=1; %If there is no monoisotopic (de)protonated ion, use the first one, the mass is calculated from CAMERA annotation +/- PROTON_MASS. Accuracy is lower.
    end
end
mono_selector = logical(mono_selector);
Monoisotope_output = [[raw_header; raw_single_annotation(mono_selector,:)],...
    [Annotation_header;num2cell(Annotation_content(mono_selector, :))], ...
    ['ambiguity_flg'; num2cell(ambiguity_flg(mono_selector))]];
output_filename = sprintf('%s_Monoisotopic.xlsx', fullfile(folder, filename));
xlswrite(output_filename, Monoisotope_output)
cleaned_filename = sprintf('%s_Cleaned.xlsx', fullfile(folder, filename));
xlswrite(cleaned_filename, raw_cleaned);