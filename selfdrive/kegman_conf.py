import json
import os

class kegman_conf():
  def __init__(self, CP=None):
    self.conf = self.read_config()
    if CP is not None:
      self.init_config(CP)

  def init_config(self, CP):
    write_conf = False
    if self.conf['tuneGernby'] != "1":
      self.conf['tuneGernby'] = str(1)
      write_conf = True
	
    # only fetch Kp, Ki, Kf sR and sRC from interface.py if it's a PID controlled car
    if CP.lateralTuning.which() == 'pid':
      if self.conf['Kp'] == "0.2":
        self.conf['Kp'] = str(round(CP.lateralTuning.pid.kpV[0],3))
        write_conf = True
      if self.conf['Ki'] == "0.01":
        self.conf['Ki'] = str(round(CP.lateralTuning.pid.kiV[0],3))
        write_conf = True
      if self.conf['Kf'] == "0.000071":
        self.conf['Kf'] = str('{:f}'.format(CP.lateralTuning.pid.kf))
        write_conf = True
    
    if self.conf['steerRatio'] == "15.7":
      self.conf['steerRatio'] = str(round(CP.steerRatio,3))
      write_conf = True
    
    if self.conf['steerRateCost'] == "1.0":
      self.conf['steerRateCost'] = str(round(CP.steerRateCost,3))
      write_conf = True

    if self.conf['steerActuatorDelay'] == "0.2":
      self.conf['steerActuatorDelay'] = str(round(CP.steerActuatorDelay, 3))
      write_conf = True

    if self.conf['steerLimitTimer'] == "2.5":
      self.conf['steerLimitTimer'] = str(round(CP.steerLimitTimer, 3))
      write_conf = True

    if write_conf:
      self.write_config(self.config)

  def read_config(self):
    self.element_updated = False

    if os.path.isfile('/data/kegman.json'):
      with open('/data/kegman.json', 'r') as f:
        self.config = json.load(f)

      if "battPercOff" not in self.config:
        self.config.update({"battPercOff":"25"})
        self.config.update({"carVoltageMinEonShutdown":"11800"})
        self.config.update({"brakeStoppingTarget":"0.25"})
        self.element_updated = True

      if "tuneGernby" not in self.config:
        self.config.update({"tuneGernby":"1"})
        self.config.update({"Kp":"0.2"})
        self.config.update({"Ki":"0.01"})
        self.element_updated = True

      if "liveParams" not in self.config:
        self.config.update({"liveParams":"1"})
        self.element_updated = True
	
      if "steerRatio" not in self.config:
        self.config.update({"steerRatio":"15.7"})
        self.config.update({"steerRateCost":"-1.0"})
        self.config.update({"steerActuatorDelay":"0.2"})
        self.element_updated = True
	
      if "leadDistance" not in self.config:
        self.config.update({"leadDistance":"5"})
        self.element_updated = True
	
      if "deadzone" not in self.config:
        self.config.update({"deadzone":"0.0"})
        self.element_updated = True
	
      if "1barBP0" not in self.config:
        self.config.update({"1barBP0":"-0.3"})
        self.config.update({"1barBP1":"2.0"})
        self.config.update({"2barBP0":"-0.3"})
        self.config.update({"2barBP1":"2.25"})
        self.config.update({"3barBP0":"0.0"})
        self.config.update({"3barBP1":"2.5"})
        self.element_updated = True


      if "1barMax" not in self.config:
        self.config.update({"1barMax":"2.1"})
        self.config.update({"2barMax":"2.1"})
        self.config.update({"3barMax":"2.1"})
        self.element_updated = True
	
      if "1barHwy" not in self.config:
        self.config.update({"1barHwy":"0.4"})
        self.config.update({"2barHwy":"0.3"})
        self.config.update({"3barHwy":"0.1"})
        self.element_updated = True
	
      if "slowOnCurves" not in self.config:
        self.config.update({"slowOnCurves":"1"})
        self.element_updated = True
	
      if "Kf" not in self.config:
        self.config.update({"Kf":"-1"})
        self.element_updated = True
	
      if "sR_boost" not in self.config:
        self.config.update({"sR_boost":"4.0"})
        self.config.update({"sR_BP0":"2.0"})
        self.config.update({"sR_BP1":"15.0"})
        self.config.update({"sR_time":"3.0"})
        self.element_updated = True

      if "ALCnudgeLess" not in self.config:
        self.config.update({"ALCnudgeLess":"1"})
        self.config.update({"ALCminSpeed":"8.6"})
        self.element_updated = True

      if "ALCtimer" not in self.config:
        self.config.update({"ALCtimer":"0.5"})
        self.element_updated = True

      if "CruiseDelta" not in self.config:
        self.config.update({"CruiseDelta":"5"})
        self.element_updated = True

      if "CruiseEnableMin" not in self.config:
        self.config.update({"CruiseEnableMin":"15"})
        self.element_updated = True

      if "AutoHold" not in self.config:
        self.config.update({"AutoHold":"1"})
        self.element_updated = True

      if "steerLimitTimer" not in self.config:
        self.config.update({"steerLimitTimer":"2.5"})
        self.element_updated = True

      if "epsModded" not in self.config:
        self.config.update({"epsModded":"0"})
        self.element_updated = True

      if "accelerationMode" not in self.config:
        self.config.update({"accelerationMode":"1"})
        self.element_updated = True

      if self.element_updated:
        print("updated")
        self.write_config(self.config)

    else:
      self.config = {"cameraOffset":"0.01", "lastTrMode":"1", "battChargeMin":"60", "battChargeMax":"70", \
                     "wheelTouchSeconds":"180000", "accelerationMode":"1","battPercOff":"25", "carVoltageMinEonShutdown":"11800", \
                     "brakeStoppingTarget":"0.25", "tuneGernby":"1", \
                     "Kp":"0.2", "Ki":"0.01", "liveParams":"1", "leadDistance":"5", "deadzone":"0.0", \
                     "1barBP0":"-0.3", "1barBP1":"2.0", "2barBP0":"-0.3", "2barBP1":"2.25", "3barBP0":"0.0", \
                     "3barBP1":"2.5", "1barMax":"2.1", "2barMax":"2.1", "3barMax":"2.1", \
                     "1barHwy":"0.4", "2barHwy":"0.3", "3barHwy":"0.1", \
                     "steerRatio":"15.7", "steerRateCost":"1.0", "steerActuatorDelay":"0.2", "slowOnCurves":"1", "Kf":"0.000071", \
                     "sR_boost":"4.0", "sR_BP0":"2.0", "sR_BP1":"15.0", "sR_time":"3.0", \
                     "ALCnudgeLess":"1", "ALCminSpeed":"8.6", "ALCtimer":"0.5", "AutoHold":"1", "CruiseDelta":"5", \
                     "CruiseEnableMin":"15", "steerLimitTimer":"2.5", "epsModded": "0"}


      self.write_config(self.config)
    return self.config

  def write_config(self, config):
    try:
      with open('/data/kegman.json', 'w') as f:
        json.dump(self.config, f, indent=2, sort_keys=True)
        os.chmod("/data/kegman.json", 0o764)
    except IOError:
      os.mkdir('/data')
      with open('/data/kegman.json', 'w') as f:
        json.dump(self.config, f, indent=2, sort_keys=True)
        os.chmod("/data/kegman.json", 0o764)
