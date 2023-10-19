'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------------------
Utils - Datetime
------------------------------------------------------------------------------

Module that provides tools for finding starting and ending times as Datetime
objects of a dataset given its filename, and for finding the time lapse in
seconds between two Datetime objects.

'''

# Loading from default packages
import datetime
#from datetime import datetime

def dataDatetimes(p):
    ''' (str) -> (datetime, datetime)

    Takes as input RCM or S1 data file path or name and returns
    start and end times as datetime objects using the file name.

    Keyword arguments: \\
    p -- data file path or name

    >>> start, end = dataDatetimes('2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv')
    >>> print(start.strftime("%b %d %Y %H:%M:%S"))
    Mar 01 2020 03:36:44
    >>> print(end.strftime("%b %d %Y %H:%M:%S"))
    Mar 03 2020 03:20:21
    '''


    # Create datetime objects for data starting and ending times
    # using the times specified in the file names
    if p[-3:]== 'trk':
        txt = p.split('_')
        start = datetime.strptime("%s%s" % (txt[9],txt[10]), '%Y%m%d%H%M%S')
        end   = datetime.strptime("%s%s" % (txt[19],txt[20]), '%Y%m%d%H%M%S')
    else:
        start = datetime.datetime(int(p[-35:-31]), # Year
                                  int(p[-31:-29]), # Month
                                  int(p[-29:-27]), # Day
                                  int(p[-27:-25]), # Hour
                                  int(p[-25:-23]), # Minute
                                  int(p[-23:-21])) # Second

        end   = datetime.datetime(int(p[-20:-16]), # Year
                                  int(p[-16:-14]), # Month
                                  int(p[-14:-12]), # Day
                                  int(p[-12:-10]), # Hour
                                  int(p[-10:-8]),  # Minute
                                  int(p[-8:-6]))   # Second

    return start, end


def dT(se):
    ''' tuple(datetime, datetime) -> float

    Takes as input start and end times as datetime objects
    and returns delta time (seconds)

    Keyword arguments: \\
    se -- (start, end) time as datetime objects

    >>> start, end = dataDatetimes('2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv')
    >>> print(dT( (start, end) ))
    171817.0
    '''

    # Unpack the input start and end times
    s, e = se

    return (e-s).total_seconds()


