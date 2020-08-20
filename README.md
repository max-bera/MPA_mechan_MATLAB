# MPA_mechan_MATLAB
MATLAB functions to calculate the Rheological properties of materials in a Micropipette Aspiration experiment

fitLinZhou.m calculates the elastic modulus of a sample from a pressure/aspirated length array.
             the sample is assumed to be incompressible and having finite size. The equation used was originally presented in: https://doi.org/10.1039/c5sm00125k
             
fitDMA.m calculates storage/loss moduli for a DMA-like aspiration test. The equations are an extension of the model used in fitLinZhou, and the same hypotheses hold true here.
