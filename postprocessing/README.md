A Python version of the post-processing is under development. In the meantime, you can find the MATLAB files used to apply the post-processing.
These are (3 files): applyCellWatershedFreq, obtain_freqAnalysis, obtain_freqAnalysis.

These functions require the DIP toolbox (https://diplib.org/) because it uses the function 'watershed' from that toolbox.

<b>STEPS</b>: 

1. Load the output image of the CNN (8-bit image, range 0-255).
2. Compute the characteristic frequency (related to the size of the most common cell).
3. Apply watershed to the CNN output. This includes: smoothing the image based on the characteristic frequency and then applying watershed.

<b>Matlab code</b>:
```
   imDense = imread([dir imEdge_name]);
   fr_edge_estm = obtain_freqAnalysis(imDense);
   flagBorder = true;
   imEdge = applyCellWatershedFreq(imDense, fr_edge_estm, flagBorder);
```

The output (imEdge) is the binary, skeletonized segmented image, from where the corneal parameters can be obtained easily.

<b>SOME CLARIFICATIONS</b>:
* The 'flagBorder' option is a new add-on (not included in the original results from the TVST manuscript). 
* The last part of the post-processing, where we combine the output of the CNN-ROI with imEdge, is not done here, but it is straightforward to check whether a superpixel (cell) in imEdge has 75% (or more) of its pixels within the ROI. If not, discard that cell when computing the biomarkers.
* These functions were originally designed for the manuscript "Vigueras-Guillén JP, Andrinopoulou ER, Engel A, Lemij HG, van Rooij J, Vermeer KA, van Vliet LJ. Corneal endothelial cell segmentation by classifier-driven merging of oversegmented images. IEEE Transactions on Medical Imaging. 2018; 37(10):2278-2289.". This simply means that some options within the functions were designed to deal with different types of images (for instance, images from the output of a stochastic watershed, or the coneal endothelial image itself), whereas now it simply deals with the output of the CNNs. It is much more robust to apply post-processing to the output of the CNN, so many of the comments (within the functions) are not relevant for our current case.
