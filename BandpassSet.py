"""
  Written by: Lynne Jones - UW - 4/28/10
   Questions or comments, email : ljones.uw@gmail.com

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

_stdX = 1.2

import os
import numpy as n
import pylab as pyl
import Bandpass
import Sed


figformat = 'png'

class BandpassSet:
    """ Set up a dictionary of a set of bandpasses (multi-filters). Be able to do things with them."""
    def __init__(self, filterlist=('u', 'g', 'r', 'i', 'z', 'y'), 
                 rootdir="thruputs", rootname="total_", rootsuffix=".dat", verbose=True):
        """Initialize filter set with filters in filterlist, in directory dir with root name root"""
        bandpass = {}
        for filter in filterlist:
            filename = os.path.join(rootdir, rootname+filter+rootsuffix)
            # read filter throughput and set up Sb/Phi and zeropoint
            if verbose:
                print "Initializing filter %s" %(filename)
            bandpass[filter] = Bandpass.Bandpass()
            bandpass[filter].readThroughput(filename)
            bandpass[filter].sbTophi()
        self.bandpass = bandpass
        self.filterlist = filterlist        
        return 

    def writePhis(self):
        """Write all phi values and wavelength to stdout"""
        # Print header.
        headerline = "#Wavelen(nm) "
        for filter in self.filterlist:
            headerline = headerline + "  phi_"  + filter
        print headerline
        # print data
        for i in range(0, len(self.bandpass[self.filterlist[0]].wavelen), 1):
            outline = "%.2f " %(self.bandpass[self.filterlist[0]].wavelen[i])
            for filter in self.filterlist:
                outline = outline + " %.6g " %(self.bandpass[filter].phi[i])
            print outline
        return

    def calcFilterEffWave(self, verbose=True):
        """Calculate the effective wavelengths for all filters."""
        effsb = {}
        effphi = {}
        for filter in self.filterlist:
            effphi[filter], effsb[filter] = self.bandpass[filter].calcEffWavelen()
        self.effsb = effsb
        self.effphi = effphi
        if verbose:
            print "Filter   Eff_Sb     Eff_phi"
            for filter in self.filterlist:
                print filter, self.effsb[filter], effphi[filter]
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
        for filter in filterlist:
            # Calculate minimum and maximum wavelengths for bandpass.
            minwavelen = bandpass[filter].wavelen.min()
            maxwavelen = bandpass[filter].wavelen.max()
            # Set defaults for dropoff points.
            drop_peak_blue[filter] = minwavelen
            drop_peak_red[filter] = maxwavelen
            drop_perc_blue[filter] = minwavelen
            drop_perc_red[filter] = maxwavelen
            # Find out what current wavelength grid is being used.
            wavelenstep=(bandpass[filter].wavelen[1] - bandpass[filter].wavelen[0])
            # Find peak throughput.
            maxthruput[filter] = bandpass[filter].sb.max()
            # Find the nearest spot on the wavelength grid used for filter, for edge lookup.
            condition = (abs(bandpass[filter].wavelen - effsb[filter]) < wavelenstep/2.0)
            waveleneffsb = bandpass[filter].wavelen[condition]
            # Now find where Sb drops below 'drop_peak_thruput' of max thruput for the first time.
            # Calculate wavelength where dropoff X percent of max level.
            wavelength = waveleneffsb
            while wavelength < maxwavelen:
                sb = bandpass[filter].sb[abs(bandpass[filter].wavelen-wavelength)<wavelenstep/2.0]
                if sb <= (drop_peak*maxthruput[filter]):
                    drop_peak_red[filter] = wavelength
                    break
                wavelength = wavelength + wavelenstep
            wavelength = waveleneffsb
            while wavelength > minwavelen:
                sb = bandpass[filter].sb[abs(bandpass[filter].wavelen-wavelength)<wavelenstep/2.0]
                if sb <= (drop_peak*maxthruput[filter]):
                    drop_peak_blue[filter] = wavelength
                    break
                wavelength = wavelength - wavelenstep      
            # Calculate wavelength where dropoff X percent,  absolute value
            wavelength = waveleneffsb
            while wavelength < maxwavelen:
                sb = bandpass[filter].sb[abs(bandpass[filter].wavelen-wavelength)<wavelenstep/2.0]
                if sb <= (drop_percent):
                    drop_perc_red[filter] = wavelength
                    break
                wavelength = wavelength + wavelenstep
            wavelength = waveleneffsb
            while wavelength > minwavelen:
                sb = bandpass[filter].sb[abs(bandpass[filter].wavelen-wavelength)<wavelenstep/2.0]
                if sb <= (drop_percent):
                    drop_perc_blue[filter] = wavelength
                    break
                wavelength = wavelength - wavelenstep  
        # Print output to screen.
        if verbose:
            print "Filter  MaxThruput EffWavelen  %s_max(blue)  %s_max(red)  %s_abs(blue)  %s_abs(red)" \
                %(drop_peak, drop_peak, drop_percent, drop_percent)
            for filter in self.filterlist:
                print "%4s   %10.4f %10.4f  %12.2f  %12.2f  %12.2f  %12.2f" \
                    % (filter, maxthruput[filter], 
                       effsb[filter],  
                       drop_peak_blue[filter], 
                       drop_peak_red[filter],
                       drop_perc_blue[filter], 
                       drop_perc_red[filter])
        # Set values (dictionaries keyed by filterlist).
        self.drop_peak_red = drop_peak_red
        self.drop_peak_blue = drop_peak_blue
        self.drop_perc_red = drop_perc_red
        self.drop_perc_blue = drop_perc_blue
        return 

    def calcFilterLeaks(self, makeplot=True, savefig=False, figroot = "filter"):
        """ Calculate throughput leaks beyond peak_droplo and peak_drophi wavelengths. 
        
        According to SRD these leaks must be below 1% of peak value in any 10nm interval,
        and less than 5% of total transmission over all wavelengths beyond where thruput<0.1peak.
        Assumes wavelength is in nanometers! (because of nm requirement).
        Generates plots for each filter, as well as calculation of fleaks. """
        # go through each filter, calculate filter leaks
        filterlist = self.filterlist
        bandpass = self.bandpass
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
        colors = ('m', 'b', 'g', 'y', 'r', 'k', 'c')
        colorindex = 0
        for filter in filterlist: 
            print "====="
            print "Analyzing %s filter" %(filter)
            # find wavelength range in use.
            minwavelen = bandpass[filter].wavelen.min()
            maxwavelen = bandpass[filter].wavelen.max()
            # find out what current wavelength grid is being used
            wavelenstep= bandpass[filter].wavelen[1] - bandpass[filter].wavelen[0]
            # find the wavelength in the wavelength grid which is closest to effsb
            condition = (abs(bandpass[filter].wavelen - effsb[filter]) < wavelenstep/2.0)
            waveleneffsb = bandpass[filter].wavelen[condition]
            # calculate peak transmission
            peaktrans = bandpass[filter].sb.max()
            # calculate total transmission in proper bandpass
            condition = ((bandpass[filter].wavelen>drop_peak_blue[filter]) & 
                         (bandpass[filter].wavelen<drop_peak_red[filter]))
            temporary = bandpass[filter].sb[condition]
            totaltrans = temporary.sum()
            # calculate total transmission outside drop_peak wavelengths of peak
            condition = ((bandpass[filter].wavelen>=drop_peak_red[filter]) | 
                         (bandpass[filter].wavelen<=drop_peak_blue[filter]))
            temporary = bandpass[filter].sb[condition]
            sumthruput_outside_bandpass = temporary.sum()
            print "Total transmission through filter: %s" %(totaltrans)
            print "Transmission outside of filter edges (drop_peak): %f" %(sumthruput_outside_bandpass)
            print "Ratio of total out-of-band to in-band transmission: %f%s" \
                %(100.*sumthruput_outside_bandpass / totaltrans, "%")
            infotext = "Out-of-band/in-band transmission %.2f%s" \
                %(sumthruput_outside_bandpass/totaltrans*100., '%')
            if sumthruput_outside_bandpass > 0.05 * totaltrans: 
                print "  Does not meet SRD - This is more than 5%s of throughput outside the bandpass %s" %('%', filter)
            else:
                print "  Meets SRD - This is less than 5%s of total throughput coming through outside bandpass" %('%')
            # calculate transmission in each 10nm interval.
            sb_10nm = n.zeros(len(bandpass[filter].sb), dtype='float')
            gapsize_10nm = 10.0 # wavelen gap in nm
            meet_SRD=True
            maxsb_10nm = 0.
            maxwavelen_10nm = 0.
            for i in range(0, len(sb_10nm), 1):
                # calculate 10 nm 'smoothed' transmission
                wavelen = bandpass[filter].wavelen[i]
                condition = ((bandpass[filter].wavelen >= wavelen-gapsize_10nm/2) & 
                             (bandpass[filter].wavelen<wavelen+gapsize_10nm/2) & 
                             ((bandpass[filter].wavelen<drop_peak_blue[filter]-gapsize_10nm/2) 
                              | (bandpass[filter].wavelen>drop_peak_red[filter]+gapsize_10nm/2)))
                sb_10nm[i]  = bandpass[filter].sb[condition].mean()
                # now see if in relevant 'non-bandpass' area
                if ((bandpass[filter].wavelen[i]<drop_peak_blue[filter]) 
                    | (bandpass[filter].wavelen[i]>drop_peak_red[filter])):
                    value = sb_10nm[i] / peaktrans *100.
                    if value > maxsb_10nm:
                        maxsb_10nm = value
                        maxwavelen_10nm = bandpass[filter].wavelen[i]
                    if (sb_10nm[i] > 0.01 * peaktrans):
                        meet_SRD = False
            if meet_SRD==False:
                print "Does not meet SRD - %s has at least one region not meeting the 10nm SRD filter leak requirement (max is %f%s of peak transmission at %.1f A)" %(filter, maxsb_10nm, "%", maxwavelen_10nm)
            if makeplot:
                # make plot for this filter
                pyl.figure()
                # set colors for filter in plot 
                color = colors[colorindex] + "-"
                colorindex = colorindex + 1
                if colorindex == len(colors):
                    colorindex = 0
                pyl.plot(bandpass[filter].wavelen, bandpass[filter].thruput, color)
                pyl.plot(bandpass[filter].wavelen, sb_10nm, 'r-', linewidth=2)
                pyl.axvline(drop_peak_blue[filter], color='b', linestyle=':')
                pyl.axvline(drop_peak_red[filter], color='b', linestyle=':')
                pyl.axhline(0.01*peaktrans, color='b', linestyle=':')
                legendstring = filter + " filter thruput, 10nm average thruput in red\n"
                legendstring = legendstring + "  Peak throughput is %.2f%s\n" %(peaktrans, '%')
                legendstring = legendstring + "  Total throughput (in band) is %.2f\n" %(totaltrans)
                legendstring = legendstring + "  " + infotext
                pyl.figtext(0.25, 0.76, legendstring)
                pyl.xlabel("Wavelength (nm)")
                pyl.ylabel("Throughput (%)")
                pyl.yscale('log')
                pyl.title(filter)
                pyl.ylim(1e-6, 1)
                pyl.xlim(xmin=3000, xmax=11000)
                if savefig:
                    figname = figroot + "_" + filter + "_fleak."+ figformat            
                    pyl.savefig(figname, format=figformat)
        # end of filters
        return

    def plotFilters(self, rootdir=".", throughput=True, phi=False, 
                    plotdropoffs=False, ploteffsb=True, compare=None, savefig=False, 
                    figroot='filters', xlim=(300, 1100), ylimthruput=(0, 1), ylimphi=(0, 0.002), 
                    filter_tags='normal', leg_tag='LSST', compare_tag='', atmos=True, 
                    linestyle='-', linewidth=2, newfig=True):
        """ Plot the filter throughputs and phi's, with limits xlim/ylimthruput/ylimphi. 
        
        Optionally add comparison (another teleThruputGroup or None) throughput and phi curves.
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
        colors = ('c', 'b', 'chartreuse', 'g', 'y', 'r', 'burlywood', 'm') 
        #colors = ('r', 'b', 'r', 'b', 'r', 'b', 'r', 'b')
        if (throughput):
            if newfig:
                pyl.figure()
            # plot throughputs
            colorindex = 0
            for filter in filterlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                    colorindex=0
                pyl.plot(bandpass[filter].wavelen, bandpass[filter].sb, 
                         color=color, linestyle=linestyle, linewidth=linewidth)
            # add effective wavelengths (optional)
            if ploteffsb:
                vertline = n.arange(0, 1.2, 0.1)
                temp = vertline*0.0 + 1.0
                colorindex = 0
                for filter in filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1            
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(effsb[filter]*temp, vertline, color=color, linestyle='-')        
            # add dropoff limits if desired (usually only good with reduced x/y limits) (optional)
            if (plotdropoffs): 
                colorindex = 0 
                for filter in filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex+1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(drop_peak_red[filter]*temp, vertline, color=color, linestyle='--') 
                    pyl.plot(drop_peak_blue[filter]*temp, vertline, color=color, linestyle='--')
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
                    pyl.plot(compare.bandpass[filter].wavelen, compare.bandpass[filter].sb, 
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
                ytags = n.arange(len(filterlist), 0, -1.0, dtype=float)
                ytags = ytags*0.04 + 0.35
            else: # 'normal' tagging
                xtags = (0.17, 0.27, 0.42, 0.585, 0.677, 0.8, 0.8, 0.8)
                ytags = (0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.69, 0.65)
            index= 0
            colorindex = 0
            for filter in filterlist: 
                pyl.figtext(xtags[index], ytags[index], filter, color=colors[colorindex], 
                            va='top', size='x-large')
                index = index+1
                colorindex = colorindex + 1
                if colorindex == len(colors):
                    colorindex = 0
            # set x/y limits
            pyl.xlim(xmin=xlim[0], xmax=xlim[1])
            pyl.ylim(ymin=ylimthruput[0], ymax=ylimthruput[1])
            pyl.xlabel("Wavelength (nm)")
            pyl.ylabel("Throughput (%)")
            pyl.grid()
            if savefig:
                figname = figroot + "_thruputs." + figformat
                pyl.savefig(figname, format=figformat)
        if (phi):
            if newfig:
                pyl.figure()
            # plot LSST 'phi' curves
            colorindex = 0
            for filter in filterlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                        colorindex = 0
                pyl.plot(bandpass[filter].wavelen, bandpass[filter].phi, color=color, 
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
                    pyl.plot(effphi[filter]*temp, vertline, color=color, linestyle='-')        
            # plot comparison throughputs (optional)
            if compare!=None:
                colorindex = 0
                for filter in compare.filterlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.plot(compare.bandpass[filter].wavelen, compare.bandpass[filter].phi,
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
            for filter in filterlist: 
                pyl.figtext(xtags[index], ytags[index], filter, color=colors[colorindex], va='top')
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


