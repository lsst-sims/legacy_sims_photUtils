"""
  Questions or comments, email : ljones.uw@gmail.com

$Id$

 The point of this class is mostly for convenience when dealing with
 sets of Seds and Bandpasses. Often this convenience is needed when
 dealing with these sets from the python interpreter (when figuring out
 if a group of SEDS looks appropriate, etc.)
 
 So, a lot of these functions actually deal with plotting. 
 Still, particularly in SedSet.py you may find the methods to calculate
 magnitudes or colors of a large group of seds 
 (with a set of Bandpasses, defined in this class) useful.
 
 Many of the functions defined here are useful for testing the set
 of LSST filters (i.e. do they meet the filter leak requirements?)
 or plotting the filters (i.e. plotFilters). 
""" 


import os
import numpy as n
import pylab as pyl
import Bandpass
import Sed


# wavelength range parameters for calculations.
WAVELEN_MIN = 300   # minimum wavelength for transmission/source (nm)
WAVELEN_MAX = 1200  # maximum wavelength for transmission/source (nm)
WAVELEN_STEP = 0.1  # step size in wavelength grid (nm)

# figure format to save output figures, if desired. (can choose 'png' or 'eps' or 'pdf' or a few others). 
figformat = 'png'

class BandpassSet:
    """ Set up a dictionary of a set of bandpasses (multi-filters).
    Run various engineering tests or visualizations."""
    
    def __init__(self):
        """Initialize the class but don't do anything yet."""
        return

    def setThroughputs_SingleFiles(self, filterlist=('u', 'g', 'r', 'i', 'z', 'y'), 
                                   rootdir="./", rootname="total_", rootsuffix=".dat", verbose=True):
        """Read bandpass set with filters in filterlist, from directory rootdir with base name rootname."""
        # Set up dictionary to hold bandpass information.
        bandpass = {}
        # Loop through filters: 
        for f in filterlist:
            # Build full filename.
            filename = os.path.join(rootdir, rootname+f+rootsuffix)
            # read filter throughput and set up Sb/Phi and zeropoint
            if verbose:
                print "Reading throughput file %s" %(filename)
            # Initialize bandpass object.
            bandpass[f] = Bandpass.Bandpass()
            # Read the throughput curve, sampling onto grid of wavelen min/max/step.
            bandpass[f].readThroughput(filename, wavelen_min=WAVELEN_MIN,
                                            wavelen_max=WAVELEN_MAX,
                                            wavelen_step=WAVELEN_STEP)
            # Calculate phi as well. 
            bandpass[f].sbTophi()
        # Set data in self.
        self.bandpass = bandpass
        self.filterlist = filterlist        
        return 

    def setThroughputs_ComponentList(self, filterlist=('u', 'g', 'r', 'i', 'z', 'y'),
                                     separate_filter_complist = ('filter_u.dat', 'filter_g.dat', 'filter_r.dat',
                                                                 'filter_i.dat', 'filter_z.dat', 'filter_y.dat'),
                                     all_filter_complist = ('detector.dat', 'lens1.dat', 'lens2.dat',
                                                            'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                                                            'atmos.dat'),
                                     rootdir = "./", verbose=True):
        """Read and build bandpass set from all_filter_complist + (for each filter in filterlist) an element
        from separate_filter_complist, using data from directory rootdir. """
        # Check that inputs for filterlist and separate_filter_complist are the same length
        #  as there should be one separate_filter_complist item for each filter.
        if len(filterlist) != len(separate_filter_complist):
            raise Exception("filterlist and separate_filter_complist should be the same length.")
        # Set up dictionary to hold final bandpass information.
        bandpass = {}
        # Loop through filters.
        i = 0
        for f in filterlist:
            # Set up full filenames in a list containing all elements of final throughput curve.  
            complist = []
            # Join all 'all-filter' items.
            for cp in all_filter_complist:
                complist.append(os.path.join(rootdir, cp))
            # Add in filter-specific item.
            complist.append(os.path.join(rootdir, separate_filter_complist[i]))
            i += 1
            if verbose:
                print "Reading throughput curves ", complist, " for filter ", f
            # Initialize bandpass object.
            bandpass[f] = Bandpass.Bandpass()
            bandpass[f].readThroughputList(complist, wavelen_min=WAVELEN_MIN,
                                           wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP)
            bandpass[f].sbTophi()
        self.bandpass = bandpass
        self.filterlist = filterlist
        return
    
    def writePhis(self):
        """Write all phi values and wavelength to stdout"""
        # This is useful for getting a data file with only phi's, as requested by some science collaborations.
        # Print header.
        headerline = "#Wavelen(nm) "
        for filter in self.filterlist:
            headerline = headerline + "  phi_"  + filter
        print headerline
        # print data
        for i in range(0, len(self.bandpass[self.filterlist[0]].wavelen), 1):
            outline = "%.2f " %(self.bandpass[self.filterlist[0]].wavelen[i])
            for f in self.filterlist:
                outline = outline + " %.6g " %(self.bandpass[f].phi[i])
            print outline
        return

    def calcFilterEffWave(self, verbose=True):
        """Calculate the effective wavelengths for all filters."""
        # Set up dictionaries for effective wavelengths, as calculated for Transmission (sb) and Phi (phi).
        effsb = {}
        effphi = {}
        # Calculate values for each filter. 
        for f in self.filterlist:
            effphi[f], effsb[f] = self.bandpass[f].calcEffWavelen()
        self.effsb = effsb
        self.effphi = effphi
        if verbose:
            print "Filter  Eff_Sb   Eff_phi"
            for f in self.filterlist:
                print " %s      %.3f  %.3f" %(f, self.effsb[f], effphi[f])
        return

    def calcZeroPoints(self, gain=1.0, verbose=True):
        """Calculate the theoretical zeropoints for the bandpass, in AB magnitudes."""
        exptime = 15   # Default exposure time.
        effarea = n.pi*(6.5*100/2.0)**2   # Default effective area of primary mirror. 
        zpt = {}
        print "Filter Zeropoint"
        for f in self.filterlist:
            zpt[f] = self.bandpass[f].calcZP_t(expTime=exptime, effarea=effarea, gain=gain)
            print " %s     %.3f" %(f, zpt[f])
        return

    def calcFilterEdges(self, drop_peak=0.10, drop_percent=0.10, verbose=True):
        """Calculate the edges of each filter for Sb, at values of 'drop_*'.
        
        Values for drop_peak are fraction of max throughput, drop_percent is where the
        filter throughput drops to an absolute percent value. """
        bandpass = self.bandpass
        filterlist = self.filterlist
        try:
            effsb = self.effsb
            effphi = self.effphi
        except AttributeError:
            self.calcFilterEffWave()
            effsb = self.effsb
            effphi = self.effphi
        # Set up dictionary for effective wavelengths and X% peak_drop wavelengths.
        drop_peak_blue = {}
        drop_peak_red = {}
        drop_perc_blue = {}
        drop_perc_red = {}
        maxthruput = {}
        # Calculate values for each filter.
        for f in filterlist:
            # Calculate minimum and maximum wavelengths for bandpass.
            minwavelen = bandpass[f].wavelen.min()
            maxwavelen = bandpass[f].wavelen.max()
            # Set defaults for dropoff points.
            drop_peak_blue[f] = maxwavelen
            drop_peak_red[f] = minwavelen
            drop_perc_blue[f] = maxwavelen
            drop_perc_red[f] = minwavelen
            # Find out what current wavelength grid is being used.
            wavelenstep=(bandpass[f].wavelen[1] - bandpass[f].wavelen[0])
            # Find peak throughput.
            maxthruput[f] = bandpass[f].sb.max()
            # Find the nearest spot on the wavelength grid used for filter, for edge lookup.
            sbindex = n.where(abs(bandpass[f].wavelen - effsb[f]) < wavelenstep/2.0)
            sbindex = sbindex[0]
            # Now find where Sb drops below 'drop_peak_thruput' of max thruput for the first time.
            # Calculate wavelength where dropoff X percent of max level.
            # Start at effective wavelength, and walk outwards. 
            for i in range(sbindex, len(bandpass[f].wavelen)):
                if bandpass[f].sb[i] <= (drop_peak*maxthruput[f]):
                    drop_peak_red[f] = bandpass[f].wavelen[i]
                    break
            for i in range(sbindex, 0, -1):
                if bandpass[f].sb[i] <= (drop_peak*maxthruput[f]):
                    drop_peak_blue[f] = bandpass[f].wavelen[i]
                    break
            # Calculate wavelength where dropoff X percent,  absolute value
            for i in range(sbindex, len(bandpass[f].wavelen)):
                if bandpass[f].sb[i] <= (drop_percent):
                    drop_perc_red[f] = bandpass[f].wavelen[i]
                    break
            for i in range(sbindex, 0, -1):
                if bandpass[f].sb[i] <= (drop_percent):
                    drop_perc_blue[f] = bandpass[f].wavelen[i]
                    break
        # Print output to screen.
        if verbose:
            print "Filter  MaxThruput EffWavelen  %s_max(blue)  %s_max(red)  %s_abs(blue)  %s_abs(red)" \
                %(drop_peak, drop_peak, drop_percent, drop_percent)
            for f in self.filterlist:
                print "%4s   %10.4f %10.4f  %12.2f  %12.2f  %12.2f  %12.2f" \
                    % (f, maxthruput[f], 
                       effsb[f],  
                       drop_peak_blue[f], 
                       drop_peak_red[f],
                       drop_perc_blue[f], 
                       drop_perc_red[f])
        # Set values (dictionaries keyed by filterlist).
        self.drop_peak_red = drop_peak_red
        self.drop_peak_blue = drop_peak_blue
        self.drop_perc_red = drop_perc_red
        self.drop_perc_blue = drop_perc_blue
        return 

    def calcFilterLeaks(self, makeplot=True, savefig=False, figroot = "bandpass"):
        """ Calculate throughput leaks beyond peak_droplo and peak_drophi wavelengths. 
        
        According to SRD these leaks must be below 1% of peak value in any 10nm interval,
        and less than 5% of total transmission over all wavelengths beyond where thruput<0.1peak.
        Assumes wavelength is in nanometers! (because of nm requirement).
        Generates plots for each filter, as well as calculation of fleaks. """
        # Go through each filter, calculate filter leaks.
        filterlist = self.filterlist
        bandpass = self.bandpass
        # Check that data are already defined, or try defining it now.
        try:
            effsb = self.effsb
            drop_peak_red = self.drop_peak_red
            drop_peak_blue = self.drop_peak_blue
        except AttributeError:
            self.calcFilterEffWave(verbose=False)
            effsb = self.effsb
            self.calcFilterEdges(drop_peak=0.10, verbose=False)
            drop_peak_red = self.drop_peak_red
            drop_peak_blue = self.drop_peak_blue
        # Set up plot colors.
        colors = ('m', 'b', 'g', 'y', 'r', 'k', 'c')
        colorindex = 0
        for f in filterlist: 
            print "====="
            print "Analyzing %s filter" %(f)
            # find wavelength range in use.
            minwavelen = bandpass[f].wavelen.min()
            maxwavelen = bandpass[f].wavelen.max()
            # find out what current wavelength grid is being used
            wavelenstep= bandpass[f].wavelen[1] - bandpass[f].wavelen[0]
            # find the wavelength in the wavelength grid which is closest to effsb
            condition = (abs(bandpass[f].wavelen - effsb[f]) < wavelenstep/2.0)
            waveleneffsb = bandpass[f].wavelen[condition]
            # calculate peak transmission
            peaktrans = bandpass[f].sb.max()
            # calculate total transmission withinin proper bandpass
            condition = ((bandpass[f].wavelen>drop_peak_blue[f]) & 
                         (bandpass[f].wavelen<drop_peak_red[f]))
            temporary = bandpass[f].sb[condition]
            totaltrans = temporary.sum()
            # calculate total transmission outside drop_peak wavelengths of peak
            condition = ((bandpass[f].wavelen>=drop_peak_red[f]) | 
                         (bandpass[f].wavelen<=drop_peak_blue[f]))
            temporary = bandpass[f].sb[condition]
            sumthruput_outside_bandpass = temporary.sum()
            print "Total transmission through filter: %s" %(totaltrans)
            print "Transmission outside of filter edges (drop_peak): %f" %(sumthruput_outside_bandpass)
            print "Ratio of total out-of-band to in-band transmission: %f%s" \
                %(100.*sumthruput_outside_bandpass / totaltrans, "%")
            infotext = "Out-of-band/in-band transmission %.2f%s" \
                %(sumthruput_outside_bandpass/totaltrans*100., '%')
            if sumthruput_outside_bandpass > 0.05 * totaltrans: 
                print " Does not meet SRD-This is more than 5%s of throughput outside the bandpass %s" %('%', f)
            else:
                print " Meets SRD - This is less than 5%s of total throughput outside bandpass" %('%')
            # calculate transmission in each 10nm interval.
            sb_10nm = n.zeros(len(bandpass[f].sb), dtype='float')
            gapsize_10nm = 10.0 # wavelen gap in nm
            meet_SRD=True
            maxsb_10nm = 0.
            maxwavelen_10nm = 0.
            for i in range(0, len(sb_10nm), 1):
                # calculate 10 nm 'smoothed' transmission
                wavelen = bandpass[f].wavelen[i]
                condition = ((bandpass[f].wavelen >= wavelen-gapsize_10nm/2) & 
                             (bandpass[f].wavelen<wavelen+gapsize_10nm/2) & 
                             ((bandpass[f].wavelen<drop_peak_blue[f]-gapsize_10nm/2) 
                              | (bandpass[f].wavelen>drop_peak_red[f]+gapsize_10nm/2)))
                sb_10nm[i]  = bandpass[f].sb[condition].mean()
                # now see if in relevant 'non-bandpass' area
                if ((bandpass[f].wavelen[i]<drop_peak_blue[f]) 
                    | (bandpass[f].wavelen[i]>drop_peak_red[f])):
                    value = sb_10nm[i] / peaktrans *100.
                    if value > maxsb_10nm:
                        maxsb_10nm = value
                        maxwavelen_10nm = bandpass[f].wavelen[i]
                    if (sb_10nm[i] > 0.01 * peaktrans):
                        meet_SRD = False
            if meet_SRD==False:
                print "Does not meet SRD - %s has at least one region not meeting the 10nm SRD filter leak requirement (max is %f%s of peak transmission at %.1f A)" %(f, maxsb_10nm, "%", maxwavelen_10nm)
            if makeplot:
                # make plot for this filter
                pyl.figure()
                # set colors for filter in plot 
                color = colors[colorindex]
                colorindex = colorindex + 1
                if colorindex == len(colors):
                    colorindex = 0
                # Make lines on the plot. 
                pyl.plot(bandpass[f].wavelen, bandpass[f].sb, color=color, linestyle="-")
                pyl.plot(bandpass[f].wavelen, sb_10nm, 'r-',linewidth=2)
                pyl.axvline(drop_peak_blue[f], color='b', linestyle=':')
                pyl.axvline(drop_peak_red[f], color='b', linestyle=':')
                pyl.axhline(0.01*peaktrans, color='b', linestyle=':')
                legendstring = f + " filter thruput, 10nm average thruput in red\n"
                legendstring = legendstring + "  Peak throughput is %.2f%s\n" %(peaktrans, '%')
                legendstring = legendstring + "  Total throughput (in band) is %.2f\n" %(totaltrans)
                legendstring = legendstring + "  " + infotext
                pyl.figtext(0.25, 0.76, legendstring)
                pyl.xlabel("Wavelength (nm)")
                pyl.ylabel("Throughput (0-1)")
                pyl.yscale('log')
                pyl.title(f)
                pyl.ylim(1e-6, 1)
                pyl.xlim(xmin=300, xmax=1100)
                if savefig:
                    figname = figroot + "_" + f + "_fleak."+ figformat            
                    pyl.savefig(figname, format=figformat)
        # end of loop through filters
        return

    def plotFilters(self, rootdir=".", throughput=True, phi=False, 
                    plotdropoffs=False, ploteffsb=True, compare=None, savefig=False, 
                    figroot='bandpass', xlim=(300, 1100), ylimthruput=(0, 1), ylimphi=(0, 0.002), 
                    filter_tags='normal', leg_tag='LSST', compare_tag='', atmos=True, 
                    linestyle='-', linewidth=2, newfig=True):
        """ Plot the filter throughputs and phi's, with limits xlim/ylimthruput/ylimphi. 
        
        Optionally add comparison (another BandpassSet) throughput and phi curves.
        and show lines for % dropoffs ; filter_tags can be side or normal. """
        # check that all self variables are set up if needed
        bandpass = self.bandpass
        filterlist = self.filterlist
        try: 
            self.effsb
            self.effphi
            if plotdropoffs:
                self.drop_peak_red
                self.drop_peak_blue
        except AttributeError:
            self.calcFilterEffWave()
            if plotdropoffs:
                self.calcFilterEdges(verbose=False)
        effsb = self.effsb
        effphi = self.effphi
        if plotdropoffs:
            drop_peak_red = self.drop_peak_red
            drop_peak_blue = self.drop_peak_blue
        # read files for atmosphere and optional comparison throughputs
        if atmos:
            atmosfile = os.path.join(rootdir, 'atmos.dat')
            atmosphere = Bandpass.Bandpass()
            atmosphere.readThroughput(atmosfile)
        Xatm=1.3
        # set up colors for plot output
        colors = ('k', 'b', 'g', 'y', 'r', 'm', 'burlywood', 'k') 
        #colors = ('r', 'b', 'r', 'b', 'r', 'b', 'r', 'b')
        if (throughput):
            if newfig:
                pyl.figure()
            # plot throughputs
            colorindex = 0
            for f in filterlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                    colorindex=0
                pyl.plot(bandpass[f].wavelen, bandpass[f].sb, 
                         color=color, linestyle=linestyle, linewidth=linewidth)
            # add effective wavelengths (optional)
            if ploteffsb:
                vertline = n.arange(0, 1.2, 0.1)
                temp = vertline*0.0 + 1.0
                colorindex = 0
                for f in filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1            
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(effsb[f]*temp, vertline, color=color, linestyle='-')        
            # add dropoff limits if desired (usually only good with reduced x/y limits) (optional)
            if (plotdropoffs): 
                colorindex = 0 
                for filter in filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex+1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(drop_peak_red[f]*temp, vertline, color=color, linestyle='--') 
                    pyl.plot(drop_peak_blue[f]*temp, vertline, color=color, linestyle='--')
            # plot atmosphere (optional)
            if atmos:
                pyl.plot(atmosphere.wavelen, atmosphere.sb, 'k:')
            # plot comparison throughputs (optional)
            if compare!=None:
                colorindex = 0
                for filter in compare.filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(compare.bandpass[f].wavelen, compare.bandpass[f].sb, 
                             color=color, linestyle='--')
            # add line legend (type of filter curves)
            legendtext = "%s = solid" %(leg_tag)
            if leg_tag==None:
                legendtext = ""
            if compare!=None:
                legendtext= legendtext + "\n%s = dashed" %(compare_tag)
            if atmos: 
                legendtext = legendtext + "\nAirmass %.1f" %(Xatm)
            pyl.figtext(0.15, 0.8, legendtext)
            # add names to filter throughputs
            if filter_tags == 'side':
                xtags = n.zeros(len(filterlist), dtype=float)
                xtags = xtags + 0.15
                spacing = (0.8 - 0.1) / len(filterlist)
                ytags = n.arange(0.8, 0.1, -1*spacing, dtype=float)
                print ytags
                ytags = ytags 
            else: # 'normal' tagging
                xtags = (0.16, 0.27, 0.42, 0.585, 0.68, 0.8, 0.8, 0.8)
                ytags = (0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.69, 0.65)
            index= 0
            colorindex = 0
            for f in filterlist: 
                pyl.figtext(xtags[index], ytags[index], f, color=colors[colorindex], 
                            va='top', size='x-large')
                index = index+1
                colorindex = colorindex + 1
                if colorindex == len(colors):
                    colorindex = 0
            # set x/y limits
            pyl.xlim(xmin=xlim[0], xmax=xlim[1])
            pyl.ylim(ymin=ylimthruput[0], ymax=ylimthruput[1])
            pyl.xlabel("Wavelength (nm)")
            pyl.ylabel("Throughput (0-1)")
            pyl.grid()
            if savefig:
                figname = figroot + "_thruputs." + figformat
                pyl.savefig(figname, format=figformat)
        if (phi):
            if newfig:
                pyl.figure()
            # plot LSST 'phi' curves
            colorindex = 0
            for f in filterlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                        colorindex = 0
                pyl.plot(bandpass[f].wavelen, bandpass[f].phi, color=color, 
                         linestyle=linestyle, linewidth=linewidth)
            # add effective wavelengths for main filter set (optional)
            if ploteffsb:
                vertline = n.arange(0, .1, 0.01)
                temp = vertline*0.0 + 1.0
                colorindex = 0
                for filter in filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1            
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(effphi[f]*temp, vertline, color=color, linestyle='-')        
            # plot comparison throughputs (optional)
            if compare!=None:
                colorindex = 0
                for filter in compare.filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(compare.bandpass[f].wavelen, compare.bandpass[f].phi,
                             color=color, linestyle='--')
            # add line legend
            legendtext = "%s = solid" %(leg_tag)
            if leg_tag ==None:
                legendtext = " "
            if compare!=None:
                legendtext = legendtext + "\n%s = dashed" %(compare_tag)
            pyl.figtext(0.15, 0.78, legendtext)
            # add name tags to filters
            if filter_tags == 'side':
                xtags = n.zeros(len(filterlist), dtype=float)
                xtags = xtags + 0.15
                ytags = n.arange(len(filterlist), 0, -1.0, dtype=float)
                ytags = ytags*0.04 + 0.35
            else:
                xtags = (0.17, 0.27, 0.42, 0.585, 0.677, 0.82, 0.82, 0.82)
                ytags = (0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.60, 0.57)
            index= 0
            colorindex = 0
            for f in filterlist: 
                pyl.figtext(xtags[index], ytags[index], f, color=colors[colorindex], va='top')
                index = index+1
                colorindex = colorindex+1 
                if colorindex==len(colors):
                    colorindex=0
            # set x/y limits
            pyl.xlim(xmin=xlim[0], xmax=xlim[1])
            pyl.ylim(ymin=ylimphi[0], ymax=ylimphi[1])
            pyl.xlabel("Wavelength (nm)")
            pyl.ylabel("Phi")
            pyl.grid()
            if savefig:
                figname = figroot + "_phi." + figformat
                pyl.savefig(figname, format=figformat)
        return


