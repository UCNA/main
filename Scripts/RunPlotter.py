# Fetches run numbers from calib DB and plots vs time
# Also contains a function to fix the runs mislabeled as "2010"
# usage:
# To plot runs: python RunPlotter.py
# To fix runs, open a python session do 'from RunPlotter import *' and 'fixRunDB()'

# Can use "GetTimeForRunlist" to convert a runlist to a timelist
# (for, eg., plotting a quantity vs. time instead of run)

# SS 12/08/14

from ucnacore.EncalDB import * 
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from datetime import timedelta
from datetime import date

def fixRunDB(undoBool = False):
    # fixes some sort of bug where many runs are listed
    # as having the year 2010
    # first run of 2011: 16581 
    # first run of 2012: 18712
    # first run of 2013: 22267 (or maybe 22274)
    # set 'undoBool = True' to undo 

    conn = open_connection()
    if undoBool == False:
        conn.execute("SELECT run_number FROM run WHERE YEAR(start_time) = 2010 ORDER BY run_number")
        factor = 1
    if undoBool == True:
        conn.execute("SELECT run_number FROM run")
        factor = -1 
    results = conn.fetchall()

    for r in results:
        fixRun(r[0], 0, conn, factor)
        
    return 0
    
def fixRun(run, years_to_add = 0, conn = 0, factor = 1):  #Auxilliary script for fixRunDB
    # use factor = -1 to undo changes, 
    # for fixing a single run, use conn = 0 to automatically open a new connection

    print "Fixing run " + str(run)
    if not conn:
        conn = open_connection()
    
    if years_to_add == 0:
        if run < 16581:
            years_to_add = 0*factor
        if run > 16580 and run < 18712:
            years_to_add = 1*factor
        if run > 18711 and run < 22267:
            years_to_add = 2*factor
        if run > 22266:
            years_to_add = 3*factor
    
#    print "Adding " + str(years_to_add) + " years"
        
    conn.execute("UPDATE run SET start_time=ADDDATE(start_time, INTERVAL %i YEAR), end_time=ADDDATE(end_time, INTERVAL %i YEAR) WHERE run_number = %i"%(years_to_add, years_to_add, run))

    return 0

def findMidTime(start_time, end_time):
    timedif = end_time - start_time          # a 'timedelta' object
    _midtime_f = timedif.total_seconds()/2. # a float
    midtime =  timedelta(seconds = _midtime_f)  # a 'timedelta' again
    
    return start_time + midtime

def getTimeForRun(runnumb, conn = 0):
    if not conn:
        conn = open_connection()
    conn.execute("SELECT start_time, end_time FROM run WHERE run_number = %s"%runnumb)
    result = conn.fetchall()
#    print runnumb
    if len (result) > 0:
 #       print result[0][0]
 #       print result[0][1]
        return findMidTime(result[0][0], result[0][1])
    else:
 #       prevrun = str(int(runnumb) - 1)
  #      nextrun = str(int(runnumb) + 1)
   #     return findMidTime(getTimeForRun(prevrun), getTimeForRun(nextrun))
        return date(2005,01,01)

def getTimeForRunlist(runlist):
    conn = open_connection()
    timelist = list()
    for run in runlist:
        timelist.append(getTimeForRun(run, conn))

    return timelist
    
def createtimelist(conn):
    conn.execute("SELECT run_number, start_time, end_time FROM run")

    results = conn.fetchall()

    runlist = list()
    startlist = list()
    endlist = list()
    midlist = list()
    for r in results:
        runlist.append(r[0])
        startlist.append(r[1])
        endlist.append(r[2])
        timedif = r[2] - r[1]           # a 'timedelta' object
        _midtime_f = timedif.total_seconds()/2. # a float
        midtime =  timedelta(seconds = _midtime_f)  # a 'timedelta' again
        midlist.append( r[1] + midtime )

    return midlist, runlist, startlist, endlist

def PlotRuns(conn):
    midlist, runlist, startlist, endlist = createtimelist(conn)

    fig, ax1 = plt.subplots()
#    ax1.xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y %H:%M:%S"))
    ax1.xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
    ax1.plot_date(midlist, runlist) 
    
    # plot run duration using error-bar style lines
    runplus = list()
    runminus = list()
    for run in runlist:
        runplus.append(run+0.5)
        runminus.append(run-0.5)
#    print runplus
    ax1.vlines(startlist,runminus, runplus)
    ax1.hlines(runlist,startlist,endlist)
    
    # Draw the bottom part of the error bars
    ax1.vlines(endlist,runminus, runplus)
#    ax1.hlines(startlist,midlist-.25,midlist+.25)
    
    ax1.set_xlabel("Time") 
    ax1.set_ylabel("Run Number")
                                  
    plt.show()

def writeRunsToFile(conn):
    midlist, runlist, startlist, endlist = createtimelist(conn)
    outfile = open('RunsByDate.txt', 'w')
    for i in range(len(midlist)):
        outstring  = str(runlist[i])   + "\t"
        outstring += str(startlist[i]) + "\t "
        outstring += str(midlist[i])   + "\t "
        outstring += str(endlist[i])   + "\n "
        outfile.write(outstring)

    outfile.close()
    return 0 

if __name__ == "__main__":
    
    conn = open_connection()
#    PlotRuns(conn)
    writeRunsToFile(conn)
    
    
    
    
