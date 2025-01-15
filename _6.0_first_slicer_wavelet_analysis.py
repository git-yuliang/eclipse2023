'''
---
### python code step 7**
#### Wavelet analysis of the slicers 
- input: mrdata.npy   # the slicer data
- output: *.png 
'''
'''

# Segment 1: Import libraries and set initialization parameters
# Segment 2: Wavelet analysis, plot the result graph
# Segment 3: Calculate characteristic parameters and save the characteristic parameters
'''

#from __future__ import division    
import numpy
from matplotlib import pyplot
import pycwt as wavelet
from pycwt.helpers import find
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits
import skimage.io as io 
import struct
import imageio as iio
from skimage import data
from skimage.registration import phase_cross_correlation
from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift
import cv2


plt.style.use('seaborn-v0_8-white')
# Set global font settings
plt.rcParams.update({
    'font.family': 'serif',   # Use serif font, which is in line with Nature's style
    'font.serif': ['Times New Roman'],  # Specify font as Times New Roman
    'axes.titlesize': 14,  # Set the font size for axis titles
    'axes.labelsize': 12,  # Set the font size for axis labels
    'xtick.labelsize': 10,  # Set the font size for x-axis tick labels
    'ytick.labelsize': 10,  # Set the font size for y-axis tick labels
    'legend.fontsize': 10,  # Set the font size for the legend
    'figure.dpi': 300,  # Set image resolution to 300 DPI, which meets high-quality paper standards
    'savefig.dpi': 300,  # Use 300 DPI for saved images
})


# Step 1st: Load data and Initialize parameters for the wavelet analysis.
mdata = np.load('./output/mrdata.npy')

outputdir = './output/'
# Pre-defined radii multipliers
MULTIPLIERS = np.arange(1.1, 8.6, 0.1)

# start frame No.= 480
# end frame No. = 12480
sno, eno = 480, 12960
ts = (eno - sno) / 240

# def find_transition_index(a, b):
#     # Ensure a and b have the same length
#     if len(a) != len(b):
#         raise ValueError("a and b must have the same length.")
    
#     # Iterate through the arrays to find the first index where a[i] > b[i]
#     for i in range(len(a)):
#         if a[i] > b[i]:
#             # Check if from this index onwards, all subsequent a[i] > b[i] hold true
#             for j in range(i, len(a)):
#                 if a[j] <= b[j]:
#                     return -1  # If any a[j] <= b[j] is found, it doesn't satisfy the condition
#             return i  # Return the first index that satisfies the condition
    
#     # If no index is found that satisfies the condition, return -1
#     return -1

def find_transition_index(a, b):
    """
    Finds the first index `i` such that for all j >= i, a[j] > b[j].
    Returns -1 if no such index exists.
    """
    # Ensure a and b have the same length
    if len(a) != len(b):
        raise ValueError("a and b must have the same length.")
    
    # Start from the first index and track if the condition is met
    for i in range(len(a)):
        if a[i] > b[i]:
            # Check if all subsequent elements satisfy a[j] > b[j]
            if all(a[j] > b[j] for j in range(i, len(a))):
                return i  # Return the first index where condition holds for all subsequent elements
    
    # If no such index is found, return -1
    return -1

   
# Centroid calculation function
def calculate_centroid(period, power):
    # Calculate centroid: weighted average along the period direction
    centroid = np.sum(power * period) / np.sum(power)
    return centroid


     
def slicer_wavelet(mdata, frn):
    ############################################
    # Section 1: Initialize parameters for the wavelet analysis.
    ############################################
    """
    Initialize parameters for the wavelet analysis.
    Args:
        frn (int): Frame number.
    """
    import numpy
    from matplotlib import pyplot
    import pycwt as wavelet
    from pycwt.helpers import find
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    import numpy as np
    from astropy.io import fits
    import skimage.io as io 
    import struct
    import imageio as iio
    from skimage import data
    from skimage.registration import phase_cross_correlation
    from skimage.registration._phase_cross_correlation import _upsampled_dft
    from scipy.ndimage import fourier_shift
    import cv2
    import pywt
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy.signal import find_peaks

    # Corresponding Slicer No.
    MULTIPLIERS = np.arange(1.1, 8.6, 0.1)
    
    # Corresponding Slicer No.
    nn = int(MULTIPLIERS[frn]*10)              # save name of the wavelet result.
    snn = str(nn)           # STRING format of the save name of the wavelet result.
    lnm = (nn-1)/10         # slicer radial height (lower part)
    lnm = str(lnm)
    hn = nn/10              # slicer radial height (higher part)
    hnm = str(hn)
    
    dat = mdata[frn,:]

    ############################################
    # Section 2: Wavelet analysis, plot results
    ############################################
    
    ctf = eno - sno# 12000 # ctf: cut frame
    dt = (ts)/(ctf-1) # In seconds
    dat = dat[sno:eno]
    #dat = dat / np.abs(np.max(dat))

    title = 'Intensity of Ridx_Rlim1'
    label = 'Ridx_Rlim1'
    units = 'counts'
    t0 = 0

    # We also create a time array in years.
    N = dat.size
    t = np.arange(0, N) * dt + t0

    mother = wavelet.DOG(6)
    mother = wavelet.Morlet(6)
    s0 = 2 * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
    dj = 1 / 12  # Twelve sub-octaves per octaves
    #J = 7 / dj  # Seven powers of two with dj sub-octaves
    J = (np.log2(N * dt / s0)) / dj
    #alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise
    #alpha, _, _ = 0,0,0
    alpha = 0

    # We write the following code to detrend and normalize the input data by its
    # standard deviation. Sometimes detrending is not necessary and simply
    # removing the mean value is good enough. However, if your dataset has a well
    # defined trend, such as the Mauna Loa CO\ :sub:`2` dataset available in the
    # above mentioned website, it is strongly advised to perform detrending.
    # Here, we fit a one-degree polynomial function and then subtract it from the
    # original data.
    p = np.polyfit(t - t0, dat, 1)
    dat_notrend = dat - np.polyval(p, t - t0)
    std = dat_notrend.std()         # Standard deviation
    var = std ** 2                  # Variance
    dat_norm = dat_notrend / std    # Normalized dataset

    #dat_norm = dat#_notrend / std    # Normalized dataset

    # intensity normaolized during the sample time series.
    plt.figure(figsize=[16,5])
    plt.subplot(131)
    plt.plot(dat)
    plt.title('Intensity of the slicer during the totally')
    plt.subplot(132)
    plt.plot(dat_notrend)
    plt.title('Detrend intensity')
    plt.subplot(133)
    plt.plot(dat_norm)
    plt.title('Normalized - detrend intensity')
    # plt.show()


    # https://qastack.cn/stats/134104/why-do-we-divide-by-the-standard-deviation-and-not-some-other-standardizing-fact
    #dat_norm = dat_notrend 
    # wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J,
    #                                                       mother)
    # iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std   

    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J,
                                                        mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std   

    # We calculate the normalized wavelet and Fourier power spectra, as well as
    # the Fourier equivalent periods for each wavelet scale.
    power = (numpy.abs(wave)) ** 2
    fft_power = numpy.abs(fft) ** 2
    period = 1 / freqs

    # We could stop at this point and plot our results. However we are also
    # interested in the power spectra significance test. The power is significant
    # where the ratio ``power / sig95 > 1``.
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                            significance_level=0.95,
                                            wavelet=mother)
    sig95 = numpy.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    # Then, we calculate the global wavelet spectrum and determine its
    # significance level.
    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha,
                                            significance_level=0.95, dof=dof,
                                            wavelet=mother)

    # We also calculate the scale average between 2 years and 8 years, and its
    # significance level.
    sel = find((period >= s0) & (period < 6))
    Cdelta = mother.cdelta
    scale_avg = (scales * numpy.ones((N, 1))).transpose()
    scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
    scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
    scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha,
                                                significance_level=0.95,
                                                dof=[scales[sel[0]],
                                                    scales[sel[-1]]],
                                                wavelet=mother)

    # Finally, we plot our results in four different subplots containing the
    # (i) original series anomaly and the inverse wavelet transform; (ii) the
    # wavelet power spectrum (iii) the global wavelet and Fourier spectra ; and
    # (iv) the range averaged wavelet spectrum. In all sub-plots the significance
    # levels are either included as dotted lines or as filled contour lines.

    # Prepare the figure
    pyplot.close('all')
    pyplot.ioff()


    plt.style.use('seaborn')

    figprops = dict(figsize=(11, 8), dpi=72)
    fig = pyplot.figure(**figprops)

    # First sub-plot, the original time series anomaly and inverse wavelet
    # transform.
    ax = pyplot.axes([0.1, 0.75, 0.65, 0.2])
    #ax.plot(t, iwave, '-', linewidth=1, color=[0.5, 0.5, 0.5])
    ax.plot(t, dat_norm, 'k', linewidth=1.5)
    ax.set_title('a) {}'.format(title))
    ax.set_ylabel(r'{} [{}]'.format(label, units))
    #plt.ylim(140,148)

    # Second sub-plot, the normalized wavelet power spectrum and significance
    # level contour lines and cone of influece hatched area. Note that period
    # scale is logarithmic.
    bx = pyplot.axes([0.1, 0.37, 0.65, 0.28], sharex=ax)
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    bx.contourf(t, numpy.log2(period), numpy.log2(power), numpy.log2(levels),
                extend='both', cmap=pyplot.cm.viridis)
    extent = [t.min(), t.max(), 0, max(period)]
    bx.contour(t, numpy.log2(period), sig95, [-99, 1], colors='k', linewidths=2,
            extent=extent)
    bx.fill(numpy.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                            t[:1] - dt, t[:1] - dt]),
            numpy.concatenate([numpy.log2(coi), [1e-9], numpy.log2(period[-1:]),
                            numpy.log2(period[-1:]), [1e-9]]),
            'k', alpha=0.3, hatch='x')
    bx.set_title('b) {} Wavelet Power Spectrum ({})'.format(label, mother.name))
    bx.set_xlabel('Time (seconds)')
    bx.set_ylabel('Period (s)')
    #fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    #
    Yticks = 2 ** numpy.arange(numpy.ceil(numpy.log2(period.min())),
                            numpy.ceil(numpy.log2(period.max())))
    bx.set_yticks(numpy.log2(Yticks))
    bx.set_yticklabels(Yticks)

    # Third sub-plot, the global wavelet and Fourier power spectra and theoretical
    # noise spectra. Note that period scale is logarithmic.
    cx = pyplot.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
    cx.plot(glbl_signif, numpy.log2(period), 'k--')
    cx.plot(var * fft_theor, numpy.log2(period), '--', color='#cccccc')
    #cx.plot(var * fft_power, numpy.log2(1./fftfreqs), '-', color='#cccccc',
            #linewidth=1.)
    cx.plot(var * glbl_power, numpy.log2(period), 'k-', linewidth=1.5)
    cx.set_title('c) Global Wavelet Spectrum')
    cx.set_xlabel(r'Power [({})^2]'.format(units))
    cx.set_xlim([0, glbl_power.max() *0.8])
    cx.set_ylim(numpy.log2([period.min(), period.max()]))
    cx.set_yticks(numpy.log2(Yticks))
    cx.set_yticklabels(Yticks)
    pyplot.setp(cx.get_yticklabels(), visible=False)

    # Fourth sub-plot, the scale averaged wavelet spectrum.
    dx = pyplot.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
    dx.axhline(scale_avg_signif, color='k', linestyle='--', linewidth=1.)
    dx.plot(t, scale_avg, 'k-', linewidth=1.5)
    dx.set_title('d) {}-{} s / UTC 11:25 Apr 20, 2023 scale-averaged power'.format(15 , 27))
    dx.set_xlabel('Time (year)')
    dx.set_ylabel(r'Average variance [{}]'.format(units))
    #plt.ylim(0,0.1)
    ax.set_xlim([t.min(), t.max()])

    # Save the figure
    pathfigs = './output/'
    fig.savefig(pathfigs + 'R' + hnm +'_full.png', bbox_inches='tight')
    plt.close()
    pyplot.close()



    # Initialize figure
    units = 'counts'
    figprops = dict(figsize=(4, 3))
    fig = plt.figure(**figprops)

    # Create the second subplot bx: Plot wavelet power spectrum
    bx = plt.axes([0.1, 0.18, 0.6, 0.7], sharex=ax)
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    bx.contourf(t, np.log2(period), np.log2(power), np.log2(levels), extend='both', cmap=plt.cm.viridis)
    extent = [t.min(), t.max(), 0, max(period)]
    bx.contour(t, np.log2(period), sig95, [-99, 1], colors='k', linewidths=2, extent=extent)
    bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),
            np.concatenate([np.log2(coi), [1e-9], np.log2(period[-1:]),
                        np.log2(period[-1:]), [1e-9]]), 'k', alpha=0.3, hatch='x')
    #bx.set_title('Within ' + lnm +  '$R_{\odot} \sim ' + hnm +  ' R_{\odot}$', fontsize=16)
    bx.set_title('Within ' + '1.006 $R_{\odot} \sim ' + hnm +  ' R_{\odot}$', fontsize=16)
    bx.set_xlabel('Time (seconds)', fontsize=16)
    bx.set_ylabel('Period (seconds)', fontsize=16)

    Yticks = 2 ** np.arange(np.ceil(np.log2(period.min())), np.ceil(np.log2(period.max())))
    bx.set_yticks(np.log2(Yticks))
    bx.set_yticklabels(Yticks)

    # Filter data within the 0 to 26 second period range
    period_range = (period >= 0) & (period <= 26)
    period_selected = period[period_range]
    glbl_power_selected = glbl_power[period_range]
    sig95_selected = sig95[period_range]

    # Filter out glbl_power data above the 95% confidence level
    time_index = 0  # Select the index for a specific time point, here assumed to be the first time point

    # Extract sig95 data at that time point
    sig95_at_time = sig95_selected[:, time_index]  # Take the confidence data for that time point

    # Filter out glbl_power greater than sig95
    valid_data = glbl_power_selected > sig95_at_time
    period_valid = period_selected[valid_data]
    glbl_power_valid = glbl_power_selected[valid_data]

    # Calculate the centroid for valid data
    centroid = calculate_centroid(period_valid, glbl_power_valid)
    print(f'Calculated centroid: {centroid:.2f} seconds')
    sumpower = np.sum(glbl_power_valid)
    print(f'Calculated sum power: {sumpower:.2f}')

    # Create cx panel, sharing the y-axis with bx
    cx = plt.axes([0.71, 0.18, 0.18, 0.7], sharey=bx)

    # Plot global power spectrum
    cx.plot(glbl_signif, np.log2(period), 'k--')
    #cx.plot(var * fft_theor, np.log2(period), '--', color='#cccccc')
    cx.plot(var * glbl_power, np.log2(period), 'k-', linewidth=1.5)

    # Draw centroid line on cx panel
    #cx.axhline(y=np.log2(centroid), color='r', linestyle='--', label='Centroid')

    # Set the x-axis of cx panel to logarithmic scale
    cx.set_xscale('log')

    # Set title and labels for cx panel
    cx.set_title('Global', fontsize=16)
    cx.set_xlabel('Power\n' + r'[({})^2]'.format(units), x=0.7, y=0.02, fontsize=8)
    cx.set_ylim(np.log2([period.min(), period.max()]))
    cx.set_yticks(np.log2(Yticks))
    cx.set_yticklabels(Yticks)
    plt.setp(cx.get_yticklabels(), visible=False)

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    pathfigs = './output/'
    fig.savefig(pathfigs + 'R' + hnm +'.png', bbox_inches='tight')
    plt.close()

    # Save centroid position
    np.save('./output/centroid' + hnm + '.npy', centroid)  # save the centroid of the global power
    np.save('./output/sumpower' + hnm + '.npy', np.sum(glbl_power_valid))  # save the total power value the global power

    ############################################
    # Section 3: Calculate characteristic parameters and save the characteristic parameters
    ############################################
    # Ensure the power array is correctly computed
    print(f"Power array shape: {power.shape}")

    if power.ndim == 1:
        glbl_power = power.mean()  # For 1D array
    elif power.ndim == 2:
        glbl_power = power.mean(axis=1)  # For 2D array
    else:
        raise ValueError(f"Unexpected power array dimensions: {power.ndim}")

    # Continue with the rest of your code


    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha,
                                            significance_level=0.95, dof=dof,
                                            wavelet=mother)

    # Consider only periods less than 26 seconds
    valid_indices = numpy.where(period < 26)[0]
    filtered_glbl_power = (var * glbl_power)[valid_indices]
    filtered_glbl_signif = glbl_signif[valid_indices]
    filtered_periods = period[valid_indices]

    # Find the period where var * glbl_power starts consistently exceeding glbl_signif


    
    
    matching_periods_idx = find_transition_index(filtered_glbl_power, filtered_glbl_signif)
    matching_periods = filtered_periods[matching_periods_idx]
    matching_periods.astype('float32')
    # greater_indices = numpy.where(filtered_glbl_power > filtered_glbl_signif)[0]

    # if len(greater_indices) > 0:
    #     start_period = filtered_periods[greater_indices[0]]  # Get the corresponding period value
    #     print(f"The period where var * glbl_power starts exceeding glbl_signif (period < 26 s): {start_period:.4f} s")
    # else:
    #     print("No matching period found (period < 26 s)")

    # # Save all matching period values
    print('start_period' + hnm + '----------->' +str(matching_periods))#matching_periods = filtered_periods[greater_indices]
    np.save('./output/' + 'start_period' + hnm + '.npy',matching_periods)
    # np.savetxt('./output/' + 'matching_periods' + hnm + '.txt', matching_periods, fmt="%.4f")
    # print("Matching period values (period < 26 s) have been saved to matching_periods.txt")

    plt.figure(num = 'start_periods', figsize=[16,5])
    plt.subplot(121)
    sn = 128
    plt.plot(np.log2(period[0:sn]) ,var * glbl_power[0:sn])
    plt.plot(np.log2(period[0:sn]), glbl_signif[0:sn], 'k--')
    plt.xlabel('Period ($log_{2}$ seconds)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.yscale('log')
    plt.subplot(122)
    plt.plot(period[0:sn], var * glbl_power[0:sn])
    plt.plot(period[0:sn], glbl_signif[0:sn], 'k--')
    plt.yscale('log')
    plt.xlabel('Period (seconds)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.tight_layout()

    pathfigs = './output/'
    plt.savefig(pathfigs + 'R' + hnm +'_start_periods.png', bbox_inches='tight')
    #plt.show()
    plt.close()


    plt.figure(num = 'periods_vs_power1', figsize=[16,5])
    plt.subplot(121)
    plt.plot(var * glbl_power)
    plt.plot(glbl_signif, 'k--')
    plt.xlabel('Period (index)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.yscale('log')
    plt.subplot(122)
    plt.plot(period, var * glbl_power)
    plt.plot(period, glbl_signif, 'k--')
    plt.yscale('log')
    plt.xlabel('Period (seconds)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.tight_layout()
    pathfigs = './output/'
    plt.savefig(pathfigs + 'R' + hnm +'_periods_vs_power1.png', bbox_inches='tight')
    #plt.show()
    plt.close()



    print(glbl_power.shape, glbl_signif.shape, period.shape)

    plt.figure(num = 'periods_vs_power2', figsize=[16,5])
    plt.subplot(121)
    plt.plot(np.log2(period[0:sn]) ,var * glbl_power[0:sn])
    plt.plot(np.log2(period[0:sn]), glbl_signif[0:sn], 'k--')
    plt.xlabel('Period ($log_{2}$ seconds)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.yscale('log')
    plt.subplot(122)
    plt.plot(period, var * glbl_power)
    plt.plot(period, glbl_signif, 'k--')
    plt.xlabel('Period (seconds)')
    plt.ylabel('Global power (Normalized $counts^{2}$)')
    plt.tight_layout()
    # plt.show()
    # save figure
    pathfigs = './output/'
    plt.savefig(pathfigs + 'R' + hnm +'_periods_vs_power2.png', bbox_inches='tight')
    #plt.show()
    plt.close()


    # Consider only periods less than 26 seconds
    valid_indices = numpy.where(period < 13)[0]
    filtered_glbl_power = (var * glbl_power)[valid_indices]
    filtered_glbl_signif = glbl_signif[valid_indices]
    filtered_periods = period[valid_indices]

    # Find where filtered_glbl_power > filtered_glbl_signif
    greater_indices = numpy.where(filtered_glbl_power > filtered_glbl_signif)[0]

    if len(greater_indices) > 0:
        # Extract the corresponding periods and global power values
        matching_periods = filtered_periods[greater_indices]
        matching_power = filtered_glbl_power[greater_indices]

        # Find local peaks in the matching power values
        peaks_indices, _ = find_peaks(matching_power)
        peaks_periods = matching_periods[peaks_indices]
        peaks_power = matching_power[peaks_indices]

        # Combine peaks and power into a single array
        peaks_data = numpy.column_stack((peaks_periods, peaks_power))

        # Print the results
        print("Local peaks (period < 26s, filtered_glbl_power > filtered_glbl_signif):")
        for period, power in peaks_data:
            print(f"Period: {period:.4f}s, Power: {power:.4f}")
            
        print(peaks_periods, peaks_power)
        # Save the results to a file
        #np.save('./output/local_peaks'+ hnm +'.npy', peaks_data, fmt="%.4f", header="Period(s) Power")
        np.save('./output/peaks_periods'+ hnm +'.npy', peaks_periods)
        np.save('./output/peaks_power'+ hnm +'.npy', peaks_power)
        
    else:
        print("No matching periods found (period < 26s)")

m,n = mdata.shape
for i in range(0,1):
    # corresponding the average value of the polar form image: [ridx:r11,:]
    slicer_wavelet(mdata, i)
    print('-----------------> slicer:',i+1,'analysis over. <-----------------')