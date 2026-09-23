function [OCP, dUdT] = computeOCPSilicon_Delithiation_Chandresekaran2011(c, T, cmax)
%
% References
% ----------
% .. [1] Chandrasekaran, R., and Fuller, T. F. (2011). Analysis of the Lithium-Ion Insertion
%    Silicon Composite Electrode/Separator/Lithium Foil Cell. Journal of The Electrochemical
%    Society, 158(8), A859-A871. DOI: 10.1149/1.3589301.


    theta = c/cmax;

    %fname = fullfile('ParameterData','BatteryCellParameters',...
    %                'LithiumIonBatteryCell','lithium_ion_battery_lco_silicon.json');
    %jsonstruct = parseBattmoJson(fname);
    %R_delith = jsonstruct.NegativeElectrode.ActiveMaterial.SolidDiffusion.rp;
    %molarVolumeSi = 1.2e-05;
    %molarVolumeLi = 8.8e-06;
%
    %Q = (3.75.*molarVolumeLi)./(molarVolumeSi);
%
    %radius = computeRadius(c,cmax,R_delith);
%
    %theta = c_ratio .* ((radius ./ R_delith).^3) ./(1+Q);


data = {{0.0010162601626016452,1.1031664964249233},{0.02032520325203246,1.1031664964249233},{0.035569105691056924,1.1031664964249233},{0.05386178861788615,1.0824584568624032},{0.06504065040650403,1.039584942325419},{0.07113821138211379,0.9856673061112632},{0.07723577235772355,0.9231694860360248},{0.0843495934959349,0.8545454545454545},{0.08943089430894313,0.8030768080918806},{0.09959349593495936,0.7516206183512294},{0.11077235772357724,0.7161015471237449},{0.11585365853658534,0.6621814195670046},{0.13617886178861788,0.6033956999426991},{0.1615853658536585,0.5495253992376491},{0.19207317073170727,0.5177308851740202},{0.2164634146341463,0.4932758663643837},{0.2357723577235772,0.4700341313934079},{0.2540650406504065,0.4529186078377637},{0.2865853658536585,0.4321607414235532},{0.3180894308943089,0.4040459403572584},{0.3302845528455284,0.3954956526071899},{0.3424796747967479,0.3954956526071899},{0.35873983739837395,0.3954956526071899},{0.3800813008130081,0.3954956526071899},{0.4308943089430893,0.3954956526071899},{0.4989837398373984,0.3954956526071899},{0.5497967479674796,0.3954956526071899},{0.5945121951219512,0.3954956526071899},{0.6361788617886177,0.3954956526071899},{0.6849593495934958,0.3926879095144372},{0.7266260162601625,0.3793069084929869},{0.7662601626016259,0.3597922220284511},{0.8008130081300813,0.33291063554149314},{0.8302845528455283,0.2864047435162809},{0.8587398373983739,0.24112210070006726},{0.875,0.20316400508233867},{0.8892276422764227,0.16520092677944145},{0.9065040650406502,0.13827698746854675},{0.9186991869918696,0.11134059144472963},{0.9390243902439024,0.08932708836792136},{0.95630081300813,0.06485463016019333},{0.9684959349593494,0.0440469368942924},{0.9817073170731704,0.026918956625725517},{0.9959349593495934,0.00856772714816012}};


    theta_vals = [];
    OCP_vals = [];

    N = length(data);

    for i = 1:N
        theta_vals(end+1) = (data{i}{1});
        OCP_vals(end+1) = (data{i}{2});
    end


    OCP = interpTable(theta_vals, OCP_vals, theta);
    dUdT = 0;

end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
