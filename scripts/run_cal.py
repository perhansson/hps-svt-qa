
# daqHardReset          
# daqSoftReset          
# daqRefreshState       
# daqSetDefaults        
# daqLoadSettings       file
# daqSaveSettings       file
# daqOpenData           file
# daqCloseData          
# daqSetRunParameters   rate count
# daqSetRunState        state
# daqGetRunState        
# daqResetCounters      
# daqSendCommand        command
# daqReadStatus         
# daqGetStatus          
# daqReadConfig         
# daqVerifyConfig       
# daqSetConfig          variable arg
# daqGetConfig          variable
# daqGetSystemStatus    
# daqGetUserStatus      
# daqGetError           
# daqSendXml            xml_string
# daqDisableTimeout
# daqOpen               system id

#--- imports ---#
import pythonDaq
import time
import os
import sys
import getopt

def doRun(path):
	while True:
		os.system('rm -f '+path)
		# Reset counters (including error)
		pythonDaq.daqResetCounters()
		pythonDaq.daqSetConfig('DataFile',path)
		configStatus = 0
		while configStatus != 1:
			configStatus = int(pythonDaq.daqVerifyConfig())
			print "Verifying config: " + str(configStatus)
		pythonDaq.daqSendCommand('OpenDataFile','')
		pythonDaq.daqSetRunState('Running')
		while pythonDaq.daqGetRunState() == "Running":
			time.sleep(1)
			print "Running ..."
		pythonDaq.daqCloseData()
		errorCount = int(pythonDaq.daqGetStatus("ErrorCount"),0)
		fileCount = pythonDaq.daqGetStatus("DataFileCount").split()[0]
		print "DataFileCount: " + fileCount
		print "ErrorCount: " + str(errorCount)
		if errorCount!=0:
			print "Bad run; trying again"
		else:
			break
	print "Run complete and saved to " + path
	
# Set some default values
run_number = -1
half_module_number = -1
svt_test = False
filter_enabled = False
hybrid_type = 0
cal_type = 0

# Read in the command line arguments
options, remainder = getopt.getopt(sys.argv[1:], 't:c:sfh', ['cal-type','hybrid-type','svt','enable-filter','help',])

# Parse the command line arguments
for opt, arg in options:
	if opt in ('-s', '--svt'):
		svt_test = True
	elif opt in ('-f', '--enable-filter'):
		filter_enabled = True	
	elif opt in ('-t', '--hybrid-type'):
		hybrid_type = int(arg)
	elif opt in ('-c', '--calibration-type'):
		cal_type = int(arg)
	elif opt in ('-h', '--help'):
		print "\nUsage: run_calibration.py <file prefix>"
		print "Arguments: "
		print "\t-t, --hybrid-type: Specify hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)"
		print "\t-c, --cal-type: Specify calibration type (0 for baseline only, 1 for response, 2 for shape on cal group 0, 3 for shape on all cal groups)"
		print "\t-s, --svt: Specifies the device to be tested is the SVT"
		print "\t-f, --enable-filter: Enable FIR filter"
		print "\n"
		sys.exit(0)

if (len(remainder) != 1):
	print "need output file prefix"
	sys.exit(-1)

pythonDaq.daqOpen("hps",1);

print "Starting calibration run ... "


# Configure the DAQ
pythonDaq.daqHardReset()
if svt_test:
	pythonDaq.daqLoadSettings('/u1/software/daq/config/local_defaults.xml') 
	pythonDaq.daqLoadSettings('/u1/software/daq/config/coda_groupC.xml')
	pythonDaq.daqSetConfig("cntrlFpga:TholdEnable", "False")
else:
	#pythonDaq.daqLoadSettings('/u1/software/daq/config/calibration_config.xml');
	pythonDaq.daqLoadSettings('/u1/software/daq/scripts/meeg/devboard.xml');
	#pythonDaq.daqSetConfig("cntrlFpga:AdcClkInvert", "False");

if filter_enabled:
	pythonDaq.daqSetConfig("cntrlFpga:FiltEnable",  "True") 
	print "Filter has been enabled" 
else:
	pythonDaq.daqSetConfig("cntrlFpga:FiltEnable",  "False") 
	print "Filter has been disabled" 

if hybrid_type==1:
	pythonDaq.daqSetConfig("cntrlFpga:HybridType",  "Old") 
	print "Configured for old (test run) hybrid" 
elif hybrid_type==2:
	pythonDaq.daqSetConfig("cntrlFpga:HybridType",  "New") 
	print "Configured for new (2014 run) hybrid" 
else:
	pythonDaq.daqSetConfig("cntrlFpga:HybridType",  "Old") 
	print "WARNING: no hybrid type set; use -t to specify old or new hybrid"
	print "Configured for old (test run) hybrid" 

pythonDaq.daqSoftReset();
pythonDaq.daqSetRunParameters("1000Hz",2000);
#pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:Ical", str(22))
#pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:Isha", str(100))
#pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:Vfs", str(80))
 
if (remainder[0][0] == '/'):
	data_path = remainder[0]
else:
	data_path = os.environ['PWD'] + '/' + remainder[0]

print 'Saving output to '  + data_path

pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:CalibInhibit","True")
pythonDaq.daqSetConfig("cntrlFpga:ApvTrigType", "DoubleTrig" )
pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:CalGroup","0")
file_path = data_path + '_baseline_dtrig.bin'
print("Running baseline (DoubleTrig) - dummy run")
doRun(file_path)
print("Running baseline (DoubleTrig)")
doRun(file_path)
print "Baseline run complete"

if cal_type>=1:
	if cal_type==1:
		print "Single delay, all cal groups"
		num_delay = 1
		num_calgroup = 8
	elif cal_type==2:
		print "All delays, single cal group"
		num_delay = 8
		num_calgroup = 1
	else:
		print "All delays, all cal groups"
		num_delay = 8
		num_calgroup = 8
	# Run baseline and calibrations
	pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:CalibInhibit","False")
	pythonDaq.daqSetConfig("cntrlFpga:ApvTrigType", "DoubleCalib" )
	for calgroup in range (0,num_calgroup):
		pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:CalGroup",str(calgroup))
		print("Running calibration group: " + str(calgroup))
		for delay in range (1,1+num_delay):
			print("Running delay: " + str(delay))
			pythonDaq.daqSetConfig("cntrlFpga:hybrid:apv25:Csel","Dly_" + str(delay) + "x3_125ns")
			file_path = data_path + '_cal_g' + str(calgroup) + "_d" + str(delay) + ".bin"
			doRun(file_path)

	   
	print "Calibration run complete"

