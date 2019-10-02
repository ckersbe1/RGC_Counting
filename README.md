# RGC_Counting
Semi-automated MATLAB cell counting algorithm for densely packed, overlapping retinal ganglion cells from fixed retinal flat mounts

Calvin J. Kersbergen, 20171117

ckersbe1@jhmi.edu

Calabresi Lab

Johns Hopkins University School of Medicine, Department of Neurology

Matlab R2017b

This is a semi-automated cell counting algorithm optimized for 
nuclear (Brn3a) staining of Retinal Ganglion Cells (RGCs) in flat mount 
preparations.
Selections can be from central, middle, or peripheral retina.
It entails a low-pass filter to reduce noise, grayscale automated or 
manual thresholding, and identification of cell boundaries. Based on the 
average cell size, it divides up large clusters of overlapping cells into
estimated numbers of cells within the cluster (Size Segmentation). 
Output is the number of cells and cell density in the loaded .tif file. 
Option for user to self-adjust contrast is available. 
Can be modified for other applications, but currently optimized for the
described analysis of RGCs only. 

This script is published here:
Jing, J., Smith, M. D., Kersbergen, C. J., Kam, T., Viswanathan, M., 
Martin, K., Dawson, T. M., Dawson, V. L., Zack, D. J., Whartenby, K. A. 
& Calabresi, P. A. Glial Pathology and Retinal Neurotoxicity in the 
Anterior Visual Pathway in Experimental Autoimmune Encephalomyelitis. 
Acta Neuropathologica Communications 7:125 (2019). 
https://doi.org/10.1186/s40478-019-0767-6

Additional details available upon request to the Calabresi Lab.
