# Some helper functions for parsing the run log that might
# have general applicability

def identifySourceSegments(use_empirical_list_bool=0, runlogname = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Aux/UCNA Run Log 2012.txt"):
    if use_empirical_list_bool == 0:
        newcycle = 0
        runlist = list()
        f = open(runlogname, "r")
        for line in f:
            # Reset cycle if found first run of cycle
            if newcycle:
                if line.find("*2") > -1:
                    print "First Run: " + line[1:6]
                    print "----"
                    newcycle = 0
                    runlist.append(int(line[1:6]))
                
            # identify new cycle
            if line.find("@cycle") > -1:
                if line.find("#") < 0:
                    newcycle = 1
                    print "----"
                    print line
                
    if use_empirical_list_bool == 1:
        # empirically determined list of source runs that look similar
        # commented number indicates end of that source run segment
        runlist = [21086, # -21100
                   21274, # -21291
                   21291, # -21330
                   21679, # -21702
                   21704, # -21718
                   21914, # -21939
                   22215, # -22238
                   22294, # -22306
                   22437, # -22452
                   22767, # -22793
                   22924] # -22930

    return runlist

def makeVLinesforSegments(axisIn, startrun = 0):
    segmentlist = identifySourceSegments()
    vlinelist = [run - 1 for run in segmentlist if run > startrun]  # otherwise vlines overlap markers :-p                                       
    axisIn.vlines(vlinelist, 0.0, 2.0, linestyles = 'dashed')

    return axisIn


def plotRunsByDate(runlogname = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Aux/UCNA Run Log 2012.txt"):
    
    runlist = list()
    f = open(runlogname, "r")

    for line in f:
        # Reset cycle if found first run of cycle
        if line.find("*2") > -1:
            run = int(line[1:6])
            print run
            runlist.append(run)

    print runlist

    return runlist

#def identifySegmentofRun(runnum):
#    print "test"
#    return 0


if __name__ == "__main__":
    plotRunsByDate()
