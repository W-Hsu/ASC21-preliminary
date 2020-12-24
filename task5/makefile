LOGFILE = presto_log

.PHONY: 1
1:
	$(MAKE) clean
	cd TestData1 && time python pipeline_mod.py GBT_Lband_PSR.fil > $(LOGFILE) 2>&1

.PHONY: 2
2:
	$(MAKE) clean
	cd TestData1 && time python pipeline_mod.py GBT_Lband_PSR.fil > $(LOGFILE) 2>&1

.PHONY: clean
clean:
	rm -rf TestData1/DDplan.ps TestData1/subbands/ TestData1/$(LOGFILE) > /dev/null 2>&1
	rm -rf TestData2/DDplan.ps TestData2/subbands/ TestData2/$(LOGFILE) > /dev/null 2>&1

