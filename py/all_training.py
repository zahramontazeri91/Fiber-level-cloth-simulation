# -*- coding: utf-8 -*-
"""
Generate training data for teeth case for 24 configurations
@author: zahra
"""
import os
os.chdir('F:/YarnGeneration/x64/Release')

yarnType = 'yarn8'
info = yarnType + '/train/spacing0.5x/10 7000 14000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing0.5x/00011 7000 16000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing0.5x/10100 7000 15000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing0.5x/11110 7000 15000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))

info = yarnType + '/train/spacing1.0x/10 7000 14500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.0x/00011 7000 17000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.0x/10100 7000 15500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.0x/11110 7000 16000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))

info = yarnType + '/train/spacing1.5x/10 7000 15000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.5x/00011 7000 17500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.5x/10100 7000 16000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing1.5x/11110 7000 16500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))

info = yarnType + '/train/spacing2.0x/10 7000 16000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.0x/00011 7000 18000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.0x/10100 7000 17000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.0x/11110 7000 17500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))

info = yarnType + '/train/spacing2.5x/10 7000 16500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.5x/00011 7000 18500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.5x/10100 7000 17500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing2.5x/11110 7000 18500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))

info = yarnType + '/train/spacing3.0x/10 7000 16500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing3.0x/00011 7000 19000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing3.0x/10100 7000 19500 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
info = yarnType + '/train/spacing3.0x/11110 7000 19000 0 1 -k 500 -v 150 -t 1 -s 2'
os.system('YarnGeneration 1 %s' %(info))
