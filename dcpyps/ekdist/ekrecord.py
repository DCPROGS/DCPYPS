import os
import math
import numpy as np

from dcpyps.ekdist import ekscn

class SingleChannelRecord(object):
    """
    A wrapper over a list of time intervals 
    from idealised single channel record.
    """
    def __init__(self, verbose=False):
        self.verbose = verbose
        if self.verbose:
            print('A new record initialised.')
        self.origin = None
        self.is_loaded = False
        self.record_type = None
        
        self.badopen=-1
            
    def load_SCN_file(self, infiles):
        """Load shut and open intervals from SCN file."""
        #TODO: check if infiles is valid entry: single file or a list of files 
        if isinstance(infiles, str):
            if os.path.isfile(infiles):
                infile = infiles
        elif isinstance(infiles, list):
            if os.path.isfile(infiles[0]):
                infile = infiles[0]
        #TODO: enable taking several scan files and join in a single record.
        # Just a single file could be loaded at present.
        self.header = ekscn.read_header(infile, self.verbose)
        self.itint, iampl, self.iprop = ekscn.read_data(
            infile, self.header)
        self.iampl = iampl.astype(float) * self.header['calfac2']
        self.origin = "Intervals loaded from SCN file: " + infile
        self.is_loaded = True
        self._tres = 0.0
        self.rtint, self.rampl, self.rprop = self.itint, self.iampl, self.iprop
        self._set_periods()
        
    def __repr__(self):
        """String representation of SingleChannelRecord instance."""
        if not self.is_loaded:
            str_repr = "Empty record" 
        else:
            str_repr = self.origin
            str_repr += "\nTotal number of intervals = {0:d}".format(
                len(self.itint))
            str_repr += ('\nResolution for HJC calculations = ' + 
                '{0:.1f} microseconds'.format(self._tres*1e6))
            str_repr += "\nNumber of resolved intervals = {0:d}".format(
                len(self.rtint))
            str_repr += "\nNumber of time periods = {0:d}".format(
                len(self.ptint))
            str_repr += '\n\nNumber of open periods = {0:d}'.format(len(self.opint))
            str_repr += ('\nMean and SD of open periods = {0:.9f} +/- {1:.9f} ms'.
                format(np.average(self.opint)*1000, np.std(self.opint)*1000))
            str_repr += ('\nRange of open periods from {0:.9f} ms to {1:.9f} ms'.
                format(np.min(self.opint)*1000, np.max(self.opint)*1000))
            str_repr += ('\n\nNumber of shut intervals = {0:d}'.format(len(self.shint)))
            str_repr += ('\nMean and SD of shut periods = {0:.9f} +/- {1:.9f} ms'.
                format(np.average(self.shint)*1000, np.std(self.shint)*1000))
            str_repr += ('\nRange of shut periods from {0:.9f} ms to {1:.9f} ms'.
                format(np.min(self.shint)*1000, np.max(self.shint)*1000))
        return str_repr
    
    def _set_resolution(self, tres=0.0):
        self._tres = tres
        self._impose_resolution()
        self._set_periods()
    def _get_resolution(self):
        return self._tres
    tres = property(_get_resolution, _set_resolution)
    
    def _impose_resolution(self):
        """
        Impose time resolution.
        First interval to start has to be resolvable, usable and preceded by
        an resolvable interval too. Otherwise its start will be defined by
        unresolvable interval and so will be unreliable.
        (1) A concantenated shut period starts with a good, resolvable
            shutting and ends when first good resolvable opening found.
            Length of concat shut period = sum of all durations before the
            resolved opening. Amplitude of concat shut period = 0.
        (2) A concantenated open period starts with a good, resolvable opening
            and ends when first good resolvable interval is found that
            has a different amplitude (either shut or open but different
            amplitude). Length of concat open period = sum of all concatenated
            durations. Amplitude of concatenated open period = weighted mean
            amplitude of all concat intervals.
        First interval in each concatenated group must be resolvable, but may
        be bad (in which case all group will be bad).
        """

        # Find negative intervals and set them unusable
        self.iprop[self.itint < 0] = 8
        # Find first resolvable and usable interval.
        n = np.intersect1d(np.where(self.itint > self._tres),
            np.where(self.iprop < 8))[0]
        
        # Initiat lists holding resolved intervals and their amplitudes and flags
        rtint, rampl, rprops = [], [], []
        # Set variables holding current interval values
        ttemp, otemp = self.itint[n], self.iprop[n]
        if (self.iampl[n] == 0):
            atemp = 0
        elif self.record_type == 'simulated':
            atemp = self.iampl[n]
        else:
            atemp = self.iampl[n] * self.itint[n]
        isopen = True if (self.iampl[n] != 0) else False
        n += 1

        # Iterate through all remaining intervals
        while n < (len(self.itint)):
            if self.itint[n] < self._tres: # interval is unresolvable

                if (len(self.itint) == n + 1) and self.iampl[n] == 0 and isopen:
                    rtint.append(ttemp)
#                    if self.record_type == 'simulated':
#                        rampl.append(atemp)
#                    else:
#                        rampl.append(atemp / ttemp)
                    rampl.append(atemp / ttemp)
                    rprops.append(otemp)
                    isopen = False
                    ttemp = self.itint[n]
                    atemp = 0
                    otemp = 8

                else:
                    ttemp += self.itint[n]
                    if self.iprop[n] >= 8: otemp = self.iprop[n]
                    if isopen: #self.iampl[n] != 0:
                        atemp += self.iampl[n] * self.itint[n]

            else:
                if (self.iampl[n] == 0): # next interval is resolvable shutting
                    if not isopen: # previous interval was shut
                        ttemp += self.itint[n]
                        if self.iprop[n] >= 8: otemp = self.iprop[n]
                    else: # previous interval was open
                        rtint.append(ttemp)
                        if self.record_type == 'simulated':
                            rampl.append(atemp)
                        else:
                            rampl.append(atemp / ttemp)
                        if (self.badopen > 0 and rtint[-1] > self.badopen):
                            rprops.append(8)
                        else:
                            rprops.append(otemp)
                        ttemp = self.itint[n]
                        otemp = self.iprop[n]
                        isopen = False
                else: # interval is resolvable opening
                    if not isopen:
                        rtint.append(ttemp)
                        rampl.append(0)
                        rprops.append(otemp)
                        ttemp, otemp = self.itint[n], self.iprop[n]
                        if self.record_type == 'simulated':
                            atemp = self.iampl[n]
                        else:
                            atemp = self.iampl[n] * self.itint[n]
                        isopen = True
                    else: # previous was open
                        if self.record_type == 'simulated':
                            ttemp += self.itint[n]
                            if self.iprop[n] >= 8: otemp = self.iprop[n]
                        elif (math.fabs((atemp / ttemp) - self.iampl[n]) <= 1.e-5):
                            ttemp += self.itint[n]
                            atemp += self.iampl[n] * self.itint[n]
                            if self.iprop[n] >= 8: otemp = self.iprop[n]
                        else:
                            rtint.append(ttemp)
                            rampl.append(atemp / ttemp)
                            if (self.badopen > 0 and rtint[-1] > self.badopen):
                                rprops.append(8)
                            else:
                                rprops.append(otemp)
                            ttemp, otemp = self.itint[n], self.iprop[n]
                            atemp = self.iampl[n] * self.itint[n]

            n += 1
        # end of while

        # add last interval
        if isopen:
            rtint.append(-1)
        else:
            rtint.append(ttemp)
        rprops.append(8)
        if isopen:
            if self.record_type == 'simulated':
                rampl.append(atemp)
            else:
                rampl.append(atemp / ttemp)
        else:
            rampl.append(0)

        self.rtint, self.rampl, self.rprop = rtint, rampl, rprops

    def _set_periods(self):
        """
        Separate open and shut intervals from the entire record.
        There may be many small amplitude transitions during one opening,
        each of which will count as an individual opening, so generally
        better to look at 'open periods'.
        Look for start of a group of openings i.e. any opening that has
        defined duration (i.e. usable).  A single unusable opening in a group
        makes its length undefined so it is excluded.
        NEW VERSION -ENSURES EACH OPEN PERIOD STARTS WITH SHUT-OPEN TRANSITION
        Find start of a group (open period) -valid start must have a good shut
        time followed by a good opening -if a bad opening is found as first (or
        any later) opening then the open period is abandoned altogether, and the
        next good shut time sought as start for next open period, but for the
        purposes of identifying the nth open period, rejected ones must be counted
        as an open period even though their length is undefined.
        """

        pint, pamp, popt = [], [], []
        # Remove first and last intervals if shut
        while self.rampl[-1] == 0:
            self.rtint = self.rtint[:-1]
            self.rampl = self.rampl[:-1]
            self.rprop = self.rprop[:-1]
        if self.rampl[0] == 0:
            self.rtint = self.rtint[1:]
            self.rampl = self.rampl[1:]
            self.rprop = self.rprop[1:]

        oint, oamp, oopt = self.rtint[0], self.rampl[0] * self.rtint[0], self.rprop[0]
        n = 1
        while n < len(self.rtint):
            if self.rampl[n] != 0:
                oint += self.rtint[n]
                oamp += self.rampl[n] * self.rtint[n]
                if self.rprop[n] >= 8: oopt = 8

                if n == (len(self.rtint) - 1):
                    pamp.append(oamp/oint)
                    pint.append(oint)
                    popt.append(oopt)
            else:
                # found two consequent gaps
                if oamp == 0 and self.rampl[n] == 0 and oopt < 8:
                    pint[-1] += self.rtint[n]
                # skip bad opening
                #elif (self.badopen > 0 and oint > self.badopen) or (oopt >= 8):
                elif (oopt >= 8):
                    popt[-1] = 8
                    oint, oamp, oopt = 0.0, 0.0, 0
#                    if n != (len(self.rint) - 2):
#                        n += 1
                else: # shutting terminates good opening
                    pamp.append(oamp/oint)
                    pint.append(oint)
                    popt.append(oopt)
                    oint, oamp, oopt = 0.0, 0.0, 0
                    pamp.append(0.0)
                    pint.append(self.rtint[n])
                    popt.append(self.rprop[n])
            n += 1

        self.ptint, self.pampl, self.pprop = pint, pamp, popt
        self.opint = self.ptint[0::2]
        self.opamp = self.pampl[0::2]
        self.oppro = self.pprop[0::2]
        self.shint = self.ptint[1::2]
        self.shamp = self.pampl[1::2]
        self.shpro = self.pprop[1::2]
