# some utilities used by led pulser scripts
# Simon Slutsky
# 11/08/13

from math import sqrt

class ValwError:
    def __init__(self, val, err):
        self.value = val
        self.error = err

# strip zeros and get average w error for arrayIn
def get_average_w_error(valArrayIn, errArrayIn):
    if len(valArrayIn) != len(errArrayIn):
        print "Input arrays have different lengths"
        badavgwerror = ValwError(0,0)
        return badavgwerror
    
    average = 0.0
    err = 0.0
    errsquared = 0.0
    cnt = 0;

    for row in range(len(valArrayIn)):
        if valArrayIn[row]:
            cnt +=1 
            average += valArrayIn[row] 
            err  = errArrayIn[row]
            errsquared += err*err

    average = average/cnt
    err = sqrt(errsquared)/cnt
    avgwerror = ValwError(average, err)
    
    return avgwerror

if __name__ == "__main__":
    print "Happy Maine!"
