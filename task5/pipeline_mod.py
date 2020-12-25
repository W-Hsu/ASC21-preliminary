"""
A simple pipelien for demostrating presto
Weiwei Zhu
2015-08-14
Max-Plank Institute for Radio Astronomy
zhuwwpku@gmail.com
"""
import os, sys, glob, re
import presto.sifting as sifting
from commands import getoutput
import numpy as np

import threading
w_thread_count = 12

log_thread = None
log_buffer = []
def write_log(log_filename):
    with open(log_filename, 'wt') as logfile:
        for i in log_buffer:
            logfile.write(i)

#Tutorial_Mode = True
Tutorial_Mode = False

rootname = 'Sband'
maxDM = 80 #max DM to search
Nsub = 32 #32 subbands
Nint = 64 #64 sub integration
Tres = 0.5 #ms
zmax = 0

filename = sys.argv[1]
if len(sys.argv) > 2:
    maskfile = sys.argv[2]
else:
    maskfile = None

# whsu:
# tutorial mode code, useless.
def query(question, answer, input_type):
    print "Based on output of the last step, answer the following questions:"
    Ntry = 3
    while not input_type(raw_input("%s:" % question)) == answer and Ntry > 0:
        Ntry -= 1
        print "try again..."
    if Ntry == 0:print "The correct answer is:", answer

#"""

print '''

====================Read Header======================

'''

#try:
#myfil = filterbank(filename)

readheadercmd = 'readfile %s' % filename
print readheadercmd
output = getoutput(readheadercmd)
print output
header = {}
for line in output.split('\n'):
    items = line.split("=")
    if len(items) > 1:
        header[items[0].strip()] = items[1].strip()

#print header
#except:
    #print 'failed at reading file %s.' % filename
    #sys.exit(0)


print '''

============Generate Dedispersion Plan===============

'''

try:
    Nchan = int(header['Number of channels'])
    tsamp = float(header['Sample time (us)']) * 1.e-6
    BandWidth = float(header['Total Bandwidth (MHz)'])
    fcenter = float(header['Central freq (MHz)'])
    Nsamp = int(header['Spectra per file'])

    if Tutorial_Mode:
        query("Input file has how many frequency channel?", Nchan, int)
        query("what is the total bandwidth?", BandWidth, float)
        query("what is the size of each time sample in us?", tsamp*1.e6, float)
        query("what's the center frequency?", fcenter, float)
        print 'see how these numbers are used in the next step.'
        print ''

    ddplancmd = 'DDplan.py -d %(maxDM)s -n %(Nchan)d -b %(BandWidth)s -t %(tsamp)f -f %(fcenter)f -s %(Nsub)s -o DDplan.ps' % {
            'maxDM':maxDM, 'Nchan':Nchan, 'tsamp':tsamp, 'BandWidth':BandWidth, 'fcenter':fcenter, 'Nsub':Nsub}
    print ddplancmd
    ddplanout = getoutput(ddplancmd)
    print ddplanout
    planlist = ddplanout.split('\n')
    ddplan = []
    planlist.reverse()
    for plan in planlist:
        if plan == '':
            continue
        elif plan.strip().startswith('Low DM'):
            break
        else:
            ddplan.append(plan)
    ddplan.reverse()
except:
    print 'failed at generating DDplan.'
    sys.exit(0)


if Tutorial_Mode:
    calls = 0
    for line in ddplan:
        ddpl = line.split()
        calls += int(ddpl[7])
    query("According to the DDplan, how many times in total do we have to call prepsubband?", calls, int)
    print 'see how these numbers are used in the next step.'
    print ''

print '''

================Dedisperse Subbands==================

'''

# Parallelized #1
# save log contents
log_buffer = []
# save inputs
parallel1_input = []

def parallel1(i, ddpl, lowDM, hiDM, dDM, DownSamp, NDMs, calls, Nout):
    while i<len(parallel1_input):
        print "calculatin' " + str(i)
        dml = parallel1_input[i]
        lodm = dml[0]
        print dml
        subDM = np.mean(dml)
        if maskfile:
            prepsubband = "prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -mask ../%s -o %s %s" % (subDM, Nsub, subdownsamp, maskfile, rootname, '../'+filename)
        else:
            prepsubband = "prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" % (subDM, Nsub, subdownsamp, rootname, '../'+filename)
        print prepsubband
        log_buffer[i] = getoutput(prepsubband)

        subnames = rootname+"_DM%.2f.sub[0-9]*" % subDM
        prepsubcmd = "prepsubband -nsub %(Nsub)d -lodm %(lowdm)f -dmstep %(dDM)f -numdms %(NDMs)d -numout %(Nout)d -downsamp %(DownSamp)d -o %(root)s %(subfile)s" % {
                'Nsub':Nsub, 'lowdm':lodm, 'dDM':dDM, 'NDMs':NDMs, 'Nout':Nout, 'DownSamp':datdownsamp, 'root':rootname, 'subfile':subnames}
        print prepsubcmd
        log_buffer[i] = log_buffer[i] + getoutput(prepsubcmd)
        i = i+w_thread_count

cwd = os.getcwd()
try:
    if not os.access('subbands', os.F_OK):
        os.mkdir('subbands')
    os.chdir('subbands')
    for line in ddplan:
        ddpl = line.split()
        lowDM = float(ddpl[0])
        hiDM = float(ddpl[1])
        dDM = float(ddpl[2])
        DownSamp = int(ddpl[3])
        NDMs = int(ddpl[6])
        calls = int(ddpl[7])
        Nout = Nsamp/DownSamp 
        Nout -= (Nout % 500)
        dmlist = np.split(np.arange(lowDM, hiDM, dDM), calls)

        #copy from $PRESTO/python/Dedisp.py
        subdownsamp = DownSamp/2
        datdownsamp = 2
        if DownSamp < 2: subdownsamp = datdownsamp = 1

        parallel1_input = []
        log_buffer = [None] * len(dmlist)
        threadlist = []
        for i, dml in enumerate(dmlist):
            parallel1_input.append(dml)
        for i in range(w_thread_count):
            threadlist.append(threading.Thread(target=parallel1, args=(i, ddpl, lowDM, hiDM, dDM, DownSamp, NDMs, calls, Nout)))
        for i in threadlist:
            i.start()
        for i in threadlist:
            i.join()

        # write log in a new thread
        log_thread = threading.Thread(target=write_log, args=('dedisperse.log',))
        log_thread.start()

    os.system('rm *.sub*')
    os.chdir(cwd)

except:
    print 'failed at prepsubband.'
    os.chdir(cwd)
    sys.exit(0)

print '''

================fft-search subbands==================

'''

# parallelized #2 #3
# save log contents
log_buffer = []
# save inputs
parallel2_input = []
parallel3_input = []

def parallel2(i):
    while i<len(parallel2_input):
        df = parallel2_input[i]
        fftcmd = "realfft %s" % df
        print fftcmd
        log_buffer[i] = getoutput(fftcmd)
        i = i+w_thread_count

def parallel3(i):
    while i<len(parallel3_input):
        fftf = parallel3_input[i]
        searchcmd = "accelsearch -zmax %d %s"  % (zmax, fftf)
        print searchcmd
        log_buffer[i] = getoutput(searchcmd)
        i=i+w_thread_count

try:
    os.chdir('subbands')

    parallel2_input = glob.glob("*.dat")
    log_buffer = [None] * len(parallel2_input)
    threadlist = []
    for i in range(w_thread_count):
        threadlist.append(threading.Thread(target=parallel2, args=(i,)))
    for i in threadlist:
        i.start()
    for i in threadlist:
        i.join()
    # make sure the log-writing from last step is completed
    log_thread.join()
    log_thread = threading.Thread(target=write_log, args=('fft.log',))
    log_thread.start()
    
    parallel3_input = glob.glob("*.fft")
    log_buffer = [None] * len(parallel3_input)
    threadlist = []
    for i in range(w_thread_count):
        threadlist.append(threading.Thread(target=parallel3, args=(i,)))
    for i in threadlist:
        i.start()
    for i in threadlist:
        i.join()
    # make sure the log-writing from last step is completed
    log_thread.join()
    log_thread = threading.Thread(target=write_log, args=('accelsearch.log',))
    log_thread.start()
    
    os.chdir(cwd)
except:
    print 'failed at fft search.'
    os.chdir(cwd)
    sys.exit(0)

#"""


def ACCEL_sift(zmax):
    '''
    The following code come from PRESTO's ACCEL_sift.py
    '''

    globaccel = "*ACCEL_%d" % zmax
    globinf = "*DM*.inf"
    # In how many DMs must a candidate be detected to be considered "good"
    min_num_DMs = 2
    # Lowest DM to consider as a "real" pulsar
    low_DM_cutoff = 2.0
    # Ignore candidates with a sigma (from incoherent power summation) less than this
    sifting.sigma_threshold = 4.0
    # Ignore candidates with a coherent power less than this
    sifting.c_pow_threshold = 100.0

    # If the birds file works well, the following shouldn't
    # be needed at all...  If they are, add tuples with the bad
    # values and their errors.
    #                (ms, err)
    sifting.known_birds_p = []
    #                (Hz, err)
    sifting.known_birds_f = []

    # The following are all defined in the sifting module.
    # But if we want to override them, uncomment and do it here.
    # You shouldn't need to adjust them for most searches, though.

    # How close a candidate has to be to another candidate to                
    # consider it the same candidate (in Fourier bins)
    sifting.r_err = 1.1
    # Shortest period candidates to consider (s)
    sifting.short_period = 0.0005
    # Longest period candidates to consider (s)
    sifting.long_period = 15.0
    # Ignore any candidates where at least one harmonic does exceed this power
    sifting.harm_pow_cutoff = 8.0

    #--------------------------------------------------------------

    # Try to read the .inf files first, as _if_ they are present, all of
    # them should be there.  (if no candidates are found by accelsearch
    # we get no ACCEL files...
    inffiles = glob.glob(globinf)
    candfiles = glob.glob(globaccel)
    # Check to see if this is from a short search
    if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
        dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
    else:
        dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
    dms = map(float, dmstrs)
    dms.sort()
    dmstrs = ["%.2f"%x for x in dms]

    # Read in all the candidates
    cands = sifting.read_candidates(candfiles)

    # Remove candidates that are duplicated in other ACCEL files
    if len(cands):
        cands = sifting.remove_duplicate_candidates(cands)

    # Remove candidates with DM problems
    if len(cands):
        cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

    # Remove candidates that are harmonically related to each other
    # Note:  this includes only a small set of harmonics
    if len(cands):
        cands = sifting.remove_harmonics(cands)

    # Write candidates to STDOUT
    if len(cands):
        cands.sort()#sifting.cmp_sigma)
        #for cand in cands[:1]:
            #print cand.filename, cand.candnum, cand.p, cand.DMstr
        #sifting.write_candlist(cands)
    return cands


print '''

================sifting candidates==================

'''

#try:
cwd = os.getcwd()
os.chdir('subbands')
cands = ACCEL_sift(zmax)
os.chdir(cwd)
#except:
    #print 'failed at sifting candidates.'
    #os.chdir(cwd)
    #sys.exit(0)


print '''

================folding candidates==================

'''

# save log contents
log_buffer = []
# save inputs
parallel4_input = []

def parallel4(i):
    while i<len(parallel4_input):
        cand = parallel4_input[i]
        foldcmd = "prepfold -n %(Nint)d -nsub %(Nsub)d -dm %(dm)f -p %(period)f %(filfile)s -o %(outfile)s -noxwin -nodmsearch" % {
                'Nint':Nint, 'Nsub':Nsub, 'dm':cand.DM,  'period':cand.p, 'filfile':filename, 'outfile':rootname+'_DM'+cand.DMstr} #full plots
        print foldcmd
        log_buffer[i] = getoutput(foldcmd)
        i = i+w_thread_count

try:
    cwd = os.getcwd()
    os.chdir('subbands')
    os.system('ln -s ../%s %s' % (filename, filename))

    parallel4_input = cands
    log_buffer = [None] * len(parallel4_input)
    threadlist = []
    for i in range(w_thread_count):
        threadlist.append(threading.Thread(target=parallel4, args=(i,)))
    for i in threadlist:
        i.start()
    for i in threadlist:
        i.join()
    # make sure the log-writing from last step is completed
    log_thread.join()
    log_thread = threading.Thread(target=write_log, args=('folding.log',))
    log_thread.start()

    os.chdir(cwd)
    log_thread.join()
except:
    print 'failed at folding candidates.'
    os.chdir(cwd)
    sys.exit(0)

