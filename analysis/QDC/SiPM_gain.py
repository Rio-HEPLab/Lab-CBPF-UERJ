#!/usr/bin/python

def find_distances(peaks):

    diffs_peaks = []
    for idx in range(0,(len(peaks)-1)):
        val = peaks[idx+1] - peaks[idx]
        diffs_peaks.append(val)

    return diffs_peaks

if __name__ == '__main__':

    # In pC
    peaks_LED_1 = (2.15066e+01, 3.34556e+01, 4.46849e+01, 5.57120e+01, 6.66103e+01) 
    
    diffs_peaks_LED_1 = find_distances(peaks_LED_1)

    mean = 0.
    N = 0  
    for val in diffs_peaks_LED_1:
        mean += val
        N += 1
        print "Distance: %.3f" % (val)

    # In pC
    peaks_LED_2 = (2.17635e+01, 3.28682e+01, 4.34672e+01)  
   
    diffs_peaks_LED_2 = find_distances(peaks_LED_2)

    for val in diffs_peaks_LED_2:
        mean += val
        N += 1
        print "Distance: %.3f" % (val)

    mean /= N
    print "Overall mean: %.3f from %d peak distances." % (mean,N)
  
    gain_amp = pow(10,(32./20.))
    e_C = 1.60217662e-19
    e_pC = e_C*1e+12
    
    gain_SiPM = (mean/gain_amp)/e_pC
    print "SiPM gain: %.1f" % (gain_SiPM)

