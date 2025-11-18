function [D_reordered] = bml_reorder_channels(cfg, D)
%BML_REORDER_CHANNELS Reorders channels in a fieldtrip object 
%   ft_selectdata() does not apply reordering--it just takes a subset of
%   channels. This function is meant to subset and reorder channels
%
% cfg.idxs - array of new indexes of ft_datatype_raw
%
% Use as
%   [D] = bml_rereference(cfg, D)
%
% cfg.label - Nx1 cellstr with labels. Defaults to raw.label
%
% D - ft_datatype_raw to be re-referenced
% 
% Returns
% D_reordered - ft_datatype_raw with subsetted and reordered channels


D_reordered = D; 
D_reordered.trial = cellfun(@(eeg) eeg(cfg.idxs, :), D_reordered.trial, 'UniformOutput', false);
D_reordered.label = D_reordered.label(cfg.idxs); 

end

