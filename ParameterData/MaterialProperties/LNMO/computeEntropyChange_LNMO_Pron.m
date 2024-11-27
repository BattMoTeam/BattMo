function dUdT = computeEntropyChange_LNMO_Pron(theta)
% Data from
%  @article{Pron_2019, title={Electrochemical Characterization and Solid Electrolyte Interface Modeling of LiNi0.5Mn1.5O4-Graphite Cells}, volume={166}, ISSN={1945-7111}, url={http://dx.doi.org/10.1149/2.0941910jes}, DOI={10.1149/2.0941910jes}, number={10}, journal={Journal of The Electrochemical Society}, publisher={The Electrochemical Society}, author={Pron, Vittorio Giai and Versaci, Daniele and Amici, Julia and Francia, Carlotta and Santarelli, Massimo and Bodoardo, Silvia}, year={2019}, pages={A2255–A2263} }
    
    dUdT_table = [[0e-2 ,  -0.06e-3]; ...
                 [10e-2 , -0.061e-3]; ...
                 [20e-2 , -0.095e-3]; ...
                 [30e-2 , -0.123e-3]; ...
                 [40e-2 , -0.324e-3]; ...
                 [50e-2 , -0.178e-3]; ...
                 [60e-2 , -0.168e-3]; ...
                 [70e-2 , -0.191e-3]; ...
                 [80e-2 , -0.249e-3]; ...
                 [85e-2 , -2.788e-3]; ...
                 [90e-2 , -1.297e-3]; ...
                 [100e-2, -0.954e-3]];
    
    dUdT = interpTable(dUdT_table(:, 1), dUdT_table(:, 2), theta);    
    
end


























