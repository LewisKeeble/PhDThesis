function out = model

import com.comsol.model.*
import com.comsol.model.util.*

model = mphopen('D:\\DNAISFET_Benchmark')

model.modelPath('C:\Users\lk3618\OneDrive - Imperial College London\Sims\COMSOL\DNAISFET')

ISFET_Nx = str2double(model.param.get('ISFET_Nx'));
ISFET_Ny = str2double(model.param.get('ISFET_Ny'));
%{
ISFET_Width = str2double(model.param.get('ISFET_Width'));
ISFET_Length = str2double(model.param.get('ISFET_Length'));
ElectrodeElectrode_Separation = str2double(model.param.get('ElectrodeElectrode_Separation'));
ISFET_StartSeparation = str2double(model.param.get('ISFET_StartSeparation'));
ISFETISFET_Separation = str2double(model.param.get('ISFETISFET_Separation'));
Electrode_Width = str2double(model.param.get('Electrode_Width'));
ISFETElectrode_Separation = str2double(model.param.get('ISFETElectrode_Separation'));
%}
ISFET_Entity = [6:45, 54:63, 72:121];
%Create surface integration of H+ concentration for each ISFET
model.result.table.create('tbl1', 'Table');
for i = 1:ISFET_Nx*ISFET_Ny
    %model.result.numerical.create(sprintf('av%d', i), 'AvSurface');
    %model.result.numerical(sprintf('av%d', i)).set('data', 'dset2');
    %model.result.numerical(sprintf('av%d', i)).set('intvolume', true);
    %ISFET_EntityNumber = mphselectbox(model, 'geom1', [(ISFETElectrode_Separation/2+(i-1)*(ISFET_Width+2*ISFETElectrode_Separation)) (ISFET_StartSeparation/2 +(j-1)*(ISFETISFET_Separation/2 + ISFET_StartSeparation + ISFET_Length)) -1e-6;  (3*ISFETElectrode_Separation/2+ISFET_Width+(i-1)*(2*ISFETElectrode_Separation+Electrode_Width)) (ISFET_StartSeparation+ ISFETISFET_Separation/4 + ISFET_Length+(j-1)*(ISFETISFET_Separation/2 + ISFET_StartSeparation + ISFET_Length)) 1e-6]', 'boundary');
    %model.result.numerical(sprintf('av%d', i)).selection.set([ISFET_Entity(i)]);
    %model.result.numerical(sprintf('av%d', i)).setIndex('expr', 'c_SiO', 0);
    model.result.numerical(sprintf('av%d', i)).set('table', 'tbl1');
    if i==1
        model.result.numerical(sprintf('av%d', i)).setResult;
    else
        model.result.numerical(sprintf('av%d', i)).appendResult;
    end
end

%Import table data
ISFET_cSiO = mphtable(model, 'tbl1');
mphlaunch

out = model