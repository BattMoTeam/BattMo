mne = MLXnotebookExporter();
mne.outputdir = '/home/xavier/Matlab/Projects/battmo/Documentation/pynbnotebooks';
mne.exportMLX('part_1_battery_modeling_guide.mlx', 'run', false, 'format', 'ipynb');
return

mne.outputdir = '/home/xavier/Matlab/Projects/battmo/Documentation/pynbnotebooks';

for inote = 1 : numel(mne.notebooknames)
    notebookname = mne.notebooknames{inote};
    mne.exportMLX(notebookname, 'format', 'ipynb');
end
