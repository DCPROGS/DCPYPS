#! /usr/bin/python
import sys
import math
import numpy as np
from dcpyps import dcio

class SCRecord(object):
    """
    A wrapper over a list of time intervals from idealised single channel
    record.
    
    iprops- numpy array of integer*1- holds properties of ith duration and amplitude
    (integer*1 has range -128 to +127 (bit 7 set gives -128; can use bits 0-6)
    0 = all OK;
    1 = amplitude dubious = bit 0;
    2 = amplitude fixed = bit 1;
    4 = amplitude of opening constrained (see fixamp) = bit 2;
    8 = duration unusable = bit 3; etc
    """

    def __init__(self, filenames=None, conc=None, tres=0.0, tcrit=None,
        badopen=-1, onechan=None, badend=None, 
        itint=None, iampl=None, iprops=None):

        self.record_type = None
        self.filenames = filenames
        self.badopen = badopen
        if filenames:
            self.load_from_file(filenames)
        else:
            if itint is not None: self.itint = itint
            if iampl is not None: self.iampl = iampl
            if iprops is not None: self.iprops = iprops
            
        if tres: 
            self._set_resolution(tres)
        else:
            self._set_resolution(np.amin(self.itint))
        if tcrit == None:
            self._set_tcrit(np.amax(self.itint))
        else:
            self._set_tcrit(tcrit)
        self.set_chs(self._tcrit)
        self.conc = conc
        self.onechan = onechan # opening from one channel only?
        self.badend = badend # bad shutting can terminate burst?

    def _set_resolution(self, tres):
        self._tres = tres
        self._impose_resolution()
        self._set_periods()
    def _get_resolution(self):
        return self._tres
    tres = property(_get_resolution, _set_resolution)

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
        while self.ramp[-1] != 0:
            self.rint.pop()
            self.ramp.pop()
            self.ropt.pop()
        if self.ramp[0] == 0:
            self.rint.pop(0)
            self.ramp.pop(0)
            self.ropt.pop(0)

        n = 1
        oint, oamp, oopt = self.rint[0], self.ramp[0] * self.rint[0], self.ropt[0]
#        print len(self.rint)
        while n < len(self.rint):
#            print n, 'int=', self.rint[n], 'amp=', self.ramp[n]
            if self.ramp[n] != 0:
                oint += self.rint[n]
                oamp += self.ramp[n] * self.rint[n]
                if self.ropt[n] >= 8: oopt = 8

                if n == (len(self.rint) - 1):
                    pamp.append(oamp/oint)
                    pint.append(oint)
                    popt.append(oopt)
#                    print '\t period= ', pint[-1], 'amplitude=', pamp[-1]
            else:
                # found two consequent gaps
                if oamp == 0 and self.ramp[n] == 0 and oopt < 8:
                    pint[-1] += self.rint[n]
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
#                    print '\t period= ', pint[-1], 'amplitude=', pamp[-1]
                    oint, oamp, oopt = 0.0, 0.0, 0
                    pamp.append(0.0)
                    pint.append(self.rint[n])
                    popt.append(self.ropt[n])
#                    print '\t period= ', pint[-1], 'amplitude=', pamp[-1]
            n += 1

        self.pint, self.pamp, self.popt = pint, pamp, popt
        self.opint = self.pint[0::2]
        self.opamp = self.pamp[0::2]
        self.oppro = self.popt[0::2]
        self.shint = self.pint[1::2]
        self.shamp = self.pamp[1::2]
        self.shpro = self.popt[1::2]
    def _get_periods(self):
        return self.pint
    periods = property(_get_periods, _set_periods)

    def _set_tcrit(self, tcrit):
        self._tcrit = tcrit
        self._set_bursts()
    def _get_tcrit(self):
        return self._tcrit
    tcrit = property(_get_tcrit, _set_tcrit)

    def load_from_file(self, filenames):
        self.filenames = filenames
        #TODO: enable taking several scan files and join in a single record.
        # Just a single file could be loaded at present.
        ioffset, nint, calfac, header = dcio.scn_read_header(filenames[0])
        self.itint, iampl, self.iprops = dcio.scn_read_data(
            filenames[0], header)
        self.iampl = iampl.astype(float) * calfac
        if header['iscanver'] == -103:
            self.record_type = 'simulated'
            
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
        self.iprops[self.itint < 0] = 8
        # Find first resolvable and usable interval.
        n = np.intersect1d(np.where(self.itint > self._tres),
            np.where(self.iprops < 8))[0]

        rtint, rampl, rprops = [], [], []
        ttemp, otemp = self.itint[n], self.iprops[n]
        if (self.iampl[n] == 0):
            atemp = 0
        elif self.record_type == 'simulated':
            atemp = self.iampl[n]
        else:
            atemp = self.iampl[n] * self.itint[n]
        isopen = True if (self.iampl[n] != 0) else False
        n += 1

        # Start looking for unresolvable intervals.
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
                    if self.iprops[n] >= 8: otemp = self.iprops[n]
                    if isopen: #self.iampl[n] != 0:
                        atemp += self.iampl[n] * self.itint[n]

            else:
                if (self.iampl[n] == 0): # next interval is resolvable shutting
                    if not isopen: # previous interval was shut
                        ttemp += self.itint[n]
                        if self.iprops[n] >= 8: otemp = self.iprops[n]
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
                        otemp = self.iprops[n]
                        isopen = False
                else: # interval is resolvable opening
                    if not isopen:
                        rtint.append(ttemp)
                        rampl.append(0)
                        rprops.append(otemp)
                        ttemp, otemp = self.itint[n], self.iprops[n]
                        if self.record_type == 'simulated':
                            atemp = self.iampl[n]
                        else:
                            atemp = self.iampl[n] * self.itint[n]
                        isopen = True
                    else: # previous was open
                        if self.record_type == 'simulated':
                            ttemp += self.itint[n]
                            if self.iprops[n] >= 8: otemp = self.iprops[n]
                        elif (math.fabs((atemp / ttemp) - self.iampl[n]) <= 1.e-5):
                            ttemp += self.itint[n]
                            atemp += self.iampl[n] * self.itint[n]
                            if self.iprops[n] >= 8: otemp = self.iprops[n]
                        else:
                            rtint.append(ttemp)
                            rampl.append(atemp / ttemp)
                            if (self.badopen > 0 and rtint[-1] > self.badopen):
                                rprops.append(8)
                            else:
                                rprops.append(otemp)
                            ttemp, otemp = self.itint[n], self.iprops[n]
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

        self.rint, self.ramp, self.ropt = rtint, rampl, rprops

    def print_all_record(self):
        for i in range(len(self.itint)):
            print (i, self.itint[i], self.iampl[i], self.iprops[i])

    def print_resolved_intervals(self):
        print('\n#########\nList of resolved intervals:\n')
        for i in range(len(self.rint)):
            print (i+1, self.rint[i]*1000, self.ramp[i], self.ropt[i])
            if (self.ramp[i] == 0) and (self.rint[i] > (self.tcrit)):
                print ('\n')
        print('\n###################\n\n')
        
    def print_resolved_periods(self):
        print ('tcrit=', self.tcrit)
        print('\n#########\nList of resolved periods:\n')
        for i in range(len(self.pint)):
            print (i+1, self.pint[i], self.pamp[i], self.popt[i])
            if self.pamp[i] == 0 and self.pint[i] > self.tcrit:
                print ('\n')
        print('\n###################\n\n')
        
    def _set_bursts(self):
        """
        Cut entire single channel record into clusters using critical shut time
        interval (tcrit).
        Default definition of cluster:
        (1) Doesn't require a gap > tcrit before the 1st cluster in each record;
        (2) Unusable shut time is a valid end of cluster;
        (3) Open probability of a cluster is calculated without considering
        last opening.
        """
        self._bursts = Bursts()
        burst = Burst()
        tcrit = math.fabs(self._tcrit)
        i = 0
        badend = False
        while i < (len(self.pint) - 1):
            

            if self.pamp[i] != 0:
                if not badend:
                    if self.popt[i] < 8:
                        burst.add_interval(self.pint[i], self.pamp[i])
                    else:
                        badend = True
            else: # found gap
                if ((self.pint[i] < tcrit) and (self.popt[i] < 8)):
                    if not badend:
                        burst.add_interval(self.pint[i], self.pamp[i])
#                elif self.pint[i] >= tcrit and self.popt[i] < 8:
                else:
                    if ((burst.get_openings_number() > 0) and 
                        (burst.get_openings_number() * 2 == len(burst.intervals) + 1) 
                        and (not badend)):
                        self._bursts.add_burst(burst)
                    burst = Burst()
                    badend = False
#                elif self.popt[i] >= 8:
#                    badend = True

            i += 1
        if self.pamp[i] != 0:
            burst.add_interval(self.pint[i], self.pamp[i])
            self._bursts.add_burst(burst)
        if burst.intervals and not badend:
            self._bursts.add_burst(burst)
    def _get_bursts(self):
        return self._bursts
    bursts = property(_get_bursts, _set_bursts)
     
    def set_conc(self, conc):
        self.conc = conc
        
    def set_chs(self, tcrit):
        if tcrit >= 0:
            self.chs = True # CHS vectors: yes
        else:
            self.chs = False
            
    def set_onechan(self, onechan):
        self.onechan = onechan # opening from one channel only?
        
    def set_badend(self, badend):
        self.badend = badend # bad shutting can terminate burst?
        
    def __repr__(self):
        
        str_repr = '\n\n Data loaded from file: '
        if self.filenames:
            str_repr += self.filenames[0]
        else:
            str_repr += "no file name; probably this is simulated record."
        if self.conc:
            str_repr += ('\nConcentration of agonist = {0:.3f} microMolar'.
                format(self.conc*1e6))
        else:
            str_repr += '\nConcentration unknown.'
        str_repr += ('\nResolution for HJC calculations = ' + 
            '{0:.1f} microseconds'.format(self._tres*1e6))
        if self._tcrit:
            str_repr += ('\nCritical gap length to define end of group (tcrit) ' +
                '= {0:.3f} milliseconds'.format(self._tcrit*1e3))
        str_repr += ('\n\t(defined so that all openings in a group prob ' + 
            'come from same channel)')
        if self.chs:
            str_repr += ('\nInitial and final vectors for bursts calculated as' +
                'in Colquhoun, Hawkes & Srodzinski, (1996, eqs 5.8, 5.11).\n')
        else:
            str_repr += ('\nInitial and final vectors for are calculated as ' +
                'for steady state openings and shuttings (this involves a ' +
                'slight approximation at start and end of bursts that are ' +
                'defined by shut times that have been set as bad).\n')

        
#        if self.record_type:
#            str_repr += '\n'
#            if self.record_type == 'simulated':
#                str_repr += '\nSimulated data loaded from file: '
#            elif self.record_type == 'recorded':
#                str_repr += '\nRecorded data loaded from file: '
#            str_repr += self.filenames[0]
#        else:
#            str_repr += '\nData not loaded...'

        if self._tres:
            str_repr += '\nNumber of resolved intervals = {0:d}'.format(len(self.rint))
            str_repr += '\nNumber of resolved periods = {0:d}'.format(len(self.opint) + len(self.shint))
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
            str_repr += ('\nLast shut period = {0:.9f} ms\n'.format(self.shint[-1]*1000))
        else:
            str_repr += '\nTemporal resolution not imposed...\n'

        if self._tcrit:
            str_repr += ('\nNumber of bursts = {0:d}'.format(self.bursts.count()))
            blength = self.bursts.get_length_list()
            openings = self.bursts.get_opening_num_list()
            
            if self.bursts.count() > 1:
                str_repr += ('\nAverage length = {0:.9f} ms'.
                    format(np.average(blength)*1000))
                str_repr += ('\nRange: {0:.3f}'.format(min(blength)*1000) +
                    ' to {0:.3f} millisec'.format(max(blength)*1000))
                #openings = self.get_openings_burst_list()
                str_repr += ('\nAverage number of openings= {0:.9f}\n'.
                    format(np.average(openings)))
            else:
                str_repr += ('\nBurst length = {0:.9f} ms'.
                    format(blength[0] * 1000))
                str_repr += ('\nNumber of openings= {0:.9f}\n'.
                    format(openings[0]))
        else:
            str_repr += '\nBursts not separated...\n'

        return str_repr

    def printout(self, output=sys.stdout):
        output.write('%s' % self)


class Burst(object):
    """

    """

    def __init__(self):
        self.setup()

    def setup(self):
        self.intervals = []
        self.amplitudes = []

    def add_interval(self, interval, amplitude):
        self.intervals.append(interval)
        self.amplitudes.append(amplitude)

    def concatenate_last(self, interval, amplitude):
        try:
            self.intervals[-1] += interval
        except:
            self.intervals.append(interval)
            self.amplitudes.append(amplitude)

    def get_open_intervals(self):
        return self.intervals[0::2]

    def get_shut_intervals(self):
        return self.intervals[1::2]

    def get_mean_amplitude(self):
        return np.average(self.amplitudes[0::2])

    def get_openings_number(self):
        return len(self.get_open_intervals())

    def get_openings_average_length(self):
        return np.average(self.get_open_intervals())

    def get_shuttings_average_length(self):
        return np.average(self.get_shut_intervals())

    def get_total_open_time(self):
        return np.sum(self.get_open_intervals())

    def get_total_shut_time(self):
        return np.sum(self.get_shut_intervals())

    def get_length(self):
        return np.sum(self.intervals)

    def get_popen(self):
        """
        Calculate Popen.
        """
        return self.get_total_open_time() / np.sum(self.intervals)

    def get_popen1(self):
        """
        Calculate Popen by excluding very last opening. Equal number of open
        and shut intervals are taken into account.
        """
        return ((self.get_total_open_time() - self.intervals[-1]) /
            (np.sum(self.intervals) - self.intervals[-1]))

    def get_running_mean_popen(self, N):

        if len(self.intervals)-1 > 2*N:
            openings = self.get_open_intervals()
            shuttings = self.get_shut_intervals()
            meanP = []
            for i in range(len(openings) - N):
                meanP.append(np.sum(openings[i: i+N]) /
                    (np.sum(openings[i: i+N]) + np.sum(shuttings[i: i+N])))
            return meanP
        else:
            return self.get_popen()

    def __repr__(self):
        ret_str = ('Group length = {0:.3f} ms; '.
            format(self.get_length() * 1000) +
            'number of openings = {0:d}; '.format(self.get_openings_number()) +
            'Popen = {0:.3f}'.format(self.get_popen()))
        if self.get_openings_number() > 1:
            ret_str += ('\n\t(Popen omitting last opening = {0:.3f})'.
            format(self.get_popen1()))
        ret_str += ('\n\tTotal open = {0:.3f} ms; total shut = {1:.3f} ms'.
            format(self.get_total_open_time() * 1000,
            self.get_total_shut_time() * 1000))
        return ret_str

class Bursts(object):
    """

    """

    def __init__(self):
        self.bursts = []
    def add_burst(self, burst):
        self.bursts.append(burst)

    def intervals(self):
        """
        Get all intervals in the record.
        """
        return [b.intervals for b in self.bursts]

    def get_op_lists(self):
        list = []
        for b in self.bursts:
            list.append(b.get_open_intervals())
        return list

    def get_sh_lists(self):
        list = []
        for b in self.bursts:
            list.append(b.get_shut_intervals())
        return list

    def all(self):
        return self.bursts

    def count(self):
        return len(self.bursts)

    def get_length_list(self):
        return [b.get_length() for b in self.bursts]

    def get_length_mean(self):
        return np.average(self.get_length_list())

    def get_opening_num_list(self):
        return [b.get_openings_number() for b in self.bursts]

    def get_opening_num_mean(self):
        return np.average(self.get_length_list())

    def get_opening_length_mean_list(self):
        return [np.average(b.get_open_intervals()) for b in self.bursts]
    
    def get_popen_list(self):
        return [b.get_popen1() for b in self.bursts]

    def get_mean_ampl_list(self):
        return [b.get_mean_amplitude() for b in self.bursts]
    
    def get_popen_mean(self):
        return np.average(self.get_popen_list())

    def get_long(self, minop):
        long = Bursts()
        for b in self.bursts:
            if b.get_openings_number() >= minop:
                long.add_burst(b)
        return long
    
    def remove_long_open_times(self, top):
        """Remove bursts which contain open intervals longer than specified 
        value.
        """
        cleaned = Bursts()
        for b in self.bursts:
            if max(b.get_open_intervals()) <= top:
                cleaned.add_burst(b)
        return [b.intervals for b in cleaned.bursts]


class PooledRecords(object):
    """
    A pool of single channel records of type SCRecord.
    """

    def __init__(self, records):
        """
        Parameters
        ----------
        records : a list of SCRecord type objects
        """
        self.recs = records
        
    
#TODO: move next 5 functions into Pooledrecords class

def all_popen(recs, min_op):
    all_popen = []
    for record in recs:
        clusters = record.bursts.get_long(min_op)
        all_popen.extend(clusters.get_popen_list())
    return all_popen

def all_mean_ampl(recs, min_op):
    malist = []
    for record in recs:
        clusters = record.bursts.get_long(min_op)
        malist.extend(clusters.get_mean_ampl_list())
    return malist

def opav(recs, min_op):
    opav = []
    for record in recs:
        clusters = record.bursts.get_long(min_op)
        opav.extend(clusters.get_opening_length_mean_list())
    return opav

def all_op_lists(recs, min_op):
    all = []
    for record in recs:
        clusters = record.bursts.get_long(min_op)
        all.extend(clusters.get_op_lists())
    return all

def all_sh_lists(recs, min_op):
    all = []
    for record in recs:
        clusters = record.bursts.get_long(min_op)
        all.extend(clusters.get_sh_lists())
    return all




