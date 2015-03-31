import os
import sys
import time
from datetime import datetime, timedelta

#class for personalized errors 
class MyError( Exception ): pass


#True if file exists and is not empty; false otherwise
def fileExistsAndNotEmpty(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

#True id DIR exists and is not empty; false otherwise
def dirExistsAndNotEmpty(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False


#Returns True is file is older than days specified, false otherwise
def fileIsOlderThan(fpath=None, numDays=1):
    try:
        BlahDaysAgo =  datetime.now()-timedelta(days=numDays)
        filetime = datetime.fromtimestamp(os.path.getctime(fpath))
    except OSError as ex:
        raise MyError("No file!: %s" %(ex,))
        
    else:
        return True if filetime < BlahDaysAgo else False
        


        
        


