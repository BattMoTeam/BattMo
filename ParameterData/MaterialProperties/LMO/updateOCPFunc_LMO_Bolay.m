function [OCP, dUdT] = updateOCPFunc_LMO_Bolay(c, T, cmax)
% Data from 
% @article{Mao_2020, title={Entropy Change Characteristics of the LiNi0.5Mn1.5O4 Cathode Material for Lithium-Ion Batteries}, volume={5}, ISSN={2470-1343}, url={http://dx.doi.org/10.1021/acsomega.9b03794}, DOI={10.1021/acsomega.9b03794}, number={8}, journal={ACS Omega}, publisher={American Chemical Society (ACS)}, author={Mao, Jing and Zhang, Peng and Liu, Xin and Liu, Yanxia and Shao, Guosheng and Dai, Kehua}, year={2020}, month=feb, pages={4109–4114} }
    
    Tref = 293.15;  % [K]

    soc = c./cmax;
    
    OCP = 289.99 − 336.28*soc − 164.73*tanh((soc + 0.5302)*6.824) − 0.0768*tanh((soc − 0.442)*7.617) − 0.171*tanh((soc − 0.9051)*13.16) − 0.3126*tanh((soc − 0.9908)*96.14) + 3896.375*tanh((soc − 0.36182)*0.08632);

    dUdT = 0;
    
end


























