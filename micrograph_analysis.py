#!/usr/bin/python
# 2016 Matt Iadanza - University of Leeds - Astbury Centre for Structural Molecular Biology

#   This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    # 
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    # 
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <http://www.gnu.org/licenses/>.

vers = '1.0'
import sys
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('ERROR: matplotlib is not installed ')
try:
    import numpy as np
except ImportError:
    sys.exit('ERROR: numpy is not installed ')
global dfa
dfa = 0

#------------------------------get arguments -----------------------------------------#
class Arg(object):
    _registry = []
    def __init__(self, flag, value, req):
        self._registry.append(self)
        self.flag = flag
        self.value = value
        self.req = req
errmsg = '''USAGE: micrograph_analysis.py --i <gCTF starfile>'''

def make_arg(flag, value, req):
    Argument = Arg(flag, value, req)
    if Argument.req == True:
        if Argument.flag not in sys.argv:
            print(errmsg)
            sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
    if Argument.value == True:
        try:
            test = sys.argv[sys.argv.index(Argument.flag)+1]
        except ValueError:
            if Argument.req == True:
                print(errmsg)
                sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
            elif Argument.req == False:
                return False
        except IndexError:
                print(errmsg)
                sys.exit("ERROR: argument '{0}' requires a value".format(Argument.flag))
        else:
            if Argument.value == True:
                Argument.value = sys.argv[sys.argv.index(Argument.flag)+1]
        
    if Argument.value == False:
        if Argument.flag in sys.argv:
            Argument.value = True
        else:
            Argument.value = False
    return Argument.value
#-----------------------------------------------------------------------------------#


#---- get the stats of the files - make some pretty graphs
def get_stats_make_graphs(alldata):
    global dfa
    global dpi
    data = []
    for i in alldata:
        if len(i.split()) > 3:
            data.append(i.split())
        if '_rlnDefocusU' in i:
            dfucol = int(i.split('#')[-1])-1
        if '_rlnDefocusV' in i:
            dfvcol = int(i.split('#')[-1])-1
        if '_rlnDefocusAngle' in i:
            dfacol = int(i.split('#')[-1])-1
        if '_rlnMicrographName' in i:
            namecol = int(i.split('#')[-1])-1
        if '_rlnCtfMaxResolution' in i:
            rescol = int(i.split('#')[-1])-1
            
    v,u,a,names,res = [],[],[],[],[]
    for i in data:
        u.append(float(i[dfucol]))
        v.append(float(i[dfvcol]))
        a.append(float(i[dfacol]))
        names.append(i[namecol])
        res.append(float(i[rescol]))
    amin = min(a)
    amax = max(a)
    maxfactor = 1/amax
    
    #-- make microsgraphs dictionary
    micsdic = {}
    count = 0
    for i in names:
        micsdic[i] = (u[count],v[count],a[count],res[count]) 
        count+=1
    #--- make defocus plot
    scaleda = []
    for i in a:
        scaleda.append(i*maxfactor)
    
    # -- make dfplot
    plt.subplot(311)
    plt.scatter(u, v, s=10, c=scaleda)        
    plt.xlabel('DefocusU',fontsize=10)
    plt.ylabel('DefocusV',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)
    
    #-- make astig plots
    astig = []
    count = 0
    for i in u:
        astig.append(abs(i - v[count]))
        count +=1    
    plt.subplot(312)
    n, bins, patches = plt.hist(astig, 100, facecolor='blue', alpha=0.75)
    plt.xlabel('Astigmatism',fontsize=10)
    plt.ylabel('Micrographs',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)
    
    dfumean = np.mean(u)
    dfvmean = np.mean(v)
    #-- make resolution plot
    plt.subplot(313)
    n, bins, patches = plt.hist(res, 100, facecolor='green', alpha=0.75)
    plt.xlabel('Resolution',fontsize=10)
    plt.ylabel('Micrographs',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)

    plt.tight_layout()
    plt.savefig('micrograph_analysis_{0}.png'.format(dfa))
    plt.close()
    dfa +=1 
    return (micsdic)
    

#------------------------------------------------------------------------
print("**** gCTF Micropgraph Analysis v{0} ****".format(vers))
print("2016 | Astbury Centre for Structural Molecular Biology | University of Leeds")
thefile = make_arg('--i',True,True)
alldata = open(thefile,'r').readlines()
micrographs = get_stats_make_graphs(alldata)

print('look at the pretty graphs in micrograph_analysis_0.png')
cull = raw_input('do you want to cull the micrographs (y/n)?')
if cull not in ('yes','Y','y','YES','Yes'):
    sys.exit('Finished, Goodbye')
print('leave if you do not wish to cull on criterion')
rescut = float(raw_input('Resolution cutoff: ') or 100000000)
acut = float(raw_input('Astigmatism cutoff: ') or 1000000000)
dfcut = float(raw_input('Defocus Cutoff: ') or 1000000000)

# find the bad micrographs
badmicrographs_df = []
badmicrographs_as = []
badmicrographs_res = []
for i in micrographs:
    if micrographs[i][0] > dfcut or micrographs[i][1] > dfcut:
        badmicrographs_df.append(i)
    if abs(micrographs[i][0] - micrographs[i][1]) > acut:
        badmicrographs_as.append(i)
    if micrographs[i][3] > rescut:
        badmicrographs_res.append(i)

output = open('{0}_culled.star'.format(thefile.split('.')[0]),'w')
for i in alldata:
    if len(i.split()) < 3:
        output.write('{0}\n'.format(i.strip('\n')))
        if '_rlnMicrographName' in i:
            namecol = int(i.split('#')[-1])-1
    else:
        if i.split()[namecol] not in badmicrographs_as and i.split()[namecol] not in badmicrographs_df and i.split()[namecol] not in badmicrographs_res:
            output.write('{0}\n'.format(i.strip('\n')))
output.close()
diagfile = open('bad_micrographs.txt','w')
diagfile.write('## Micrographs thrown out for defocus > {0}\n'.format(dfcut))
for i in badmicrographs_df:
    diagfile.write('{0}\n'.format(i))
diagfile.write('## Micrographs thrown out for astigmatism > {0}\n'.format(acut))
for i in badmicrographs_as:
    diagfile.write('{0}\n'.format(i))
diagfile.write('## Micrographs thrown out for resolution > {0}\n'.format(rescut))
for i in badmicrographs_res:
    diagfile.write('{0}\n'.format(i))

culledfile = '{0}_culled.star'.format(thefile.split('.')[0])
alldata = open(culledfile,'r').readlines()
micrographs = get_stats_make_graphs(alldata)
print('Finished! New graphs are in micrograph_analysis_1.png')
print('the new file is {0}'.format(culledfile))
print('list of culled micrographs is in bad_micrographs.txt')
