#!/usr/bin/env python

import os
import sys
import optparse
import commands

EOSPREFIX='root://eoscms//eos/cms'

"""
processes locally or submits to batch the RECO step for a ROOT file or a directory of ROOT files
"""
def main():

    #tag production with current hash from git
    githash=commands.getstatusoutput('git describe --always')[1]

    recoSWdir=os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))

    # configure
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',     dest='output' , 
                      help='Output directory (local or eos) [default=%default]', 
                      default='/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/RECO/%s/' % githash)
    parser.add_option('-i', '--input',      dest='input',   
                      help='Input directory or file (local or eos) [default=%default]',
                      default='/store/group/dpg_ecal/alca_ecalcalib/TimingTB_H2_Jul2015/raw/DataTree/3351')
    parser.add_option('-c', '--cfg',      dest='cfg',   
                      help='Configuration file for reconstructin [default=%default]',
                      default='%s/test/reco_cfg.json' % recoSWdir)
    parser.add_option('--base',      dest='base',
                      help='Directory with RECO software installation [default=%default]',
                      default=recoSWdir)
    parser.add_option('-q', '--queue',      dest='queue', 
                      help='Batch queue [default=%default]',
                      default='local')
    (opt, args) = parser.parse_args()

    scriptRealPath=os.path.realpath(sys.argv[0])
    if opt.queue!='local':
        os.system('bsub -q %s %s -o %s -i %s --base %s -c %s -q local'%(opt.queue,scriptRealPath,opt.output,opt.input,opt.base,opt.cfg))
        return 0

    #build list of input files
    inputList=opt.input.split(',')
    fileList=[]
    for i in inputList:
        if i.endswith('.root'):
            fileList.append(opt.input)
        else:
            listCmd='ls %s | grep -ir .root' % i
            if '/store/' in i :
                listCmd='cmsLs %s | awk \'{print $5}\'' % i
            fileList=commands.getstatusoutput(listCmd)[1].split()

    #build list of tasks
    inputUrl=''
    for f in fileList :
        if '/store/' in f and not EOSPREFIX in f:  f='%s/%s' % (EOSPREFIX,f)
        inputUrl += f +','
    inputUrl=inputUrl[:-1]

    #local output
    localOutputPrefix=''
    baseInputDir=os.path.dirname(f)
    localOutputPrefix=''       if len(baseInputDir)==0 else os.path.basename(os.path.dirname(f))
    localOutputPostfix='.root' if len(fileList)>1      else '_'+os.path.basename(fileList[0])
    localOutput='/tmp/RECO_%s%s' % ( localOutputPrefix,localOutputPostfix )

    #prepare output
    prepareCmd = 'cmsMkdir' if '/store/' in opt.output else 'mkdir -p'
    commands.getstatusoutput('%s %s' % (prepareCmd,opt.output) )
    stageCmd   = 'cmsStage -f' if '/store/' in opt.output else 'cp'
        
    #run 
    os.system('cd %s && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s/lib && ./bin/RunH4treeReco %s %s %s && %s %s %s && rm %s && cd -' 
                  % (opt.base,
                     opt.base,
                     inputUrl,opt.cfg,localOutput,
                     stageCmd,localOutput,opt.output,
                     localOutput)
                  )
        
    return 0


if __name__ == '__main__':
    main()
