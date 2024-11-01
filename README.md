# Temperature tests the SiPMs for ASHIPH SND

To get the code, you need to run the command: <br />
```
git clone https://github.com/ivovtin/temp_test_ashiph_sipm.git
```

Just write data from three channels of the digitizer: <br />
```
. start_run.sh
```

Write data from three channels of the digitizer in cycle: <br />
```
. start_meas.sh
```

Run processing the collected data: <br />
```
. start.sh
```

Draw stability temperature and LED amplitude from time: <br />
```
python3 draw_res_led.py
```

Draw stability temperature and LED amplitude from time with cut bad regions: <br />
```
python3 cut_bads_regions_draw_res_led.py
```