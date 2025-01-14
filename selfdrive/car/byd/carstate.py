from cereal import car
from opendbc.can.parser import CANParser
from opendbc.can.can_define import CANDefine
from common.numpy_fast import mean
from selfdrive.config import Conversions as CV
from selfdrive.car.interfaces import CarStateBase
from selfdrive.car.byd.values import DBC, CAR, HUD_MULTIPLIER

class CarState(CarStateBase):
  def __init__(self, CP):
    super().__init__(CP)
    can_define = CANDefine(DBC[CP.carFingerprint]['pt'])
    self.shifter_values = can_define.dv["DRIVE_STATE"]['GEAR']
    self.set_distance_values = can_define.dv['ACC_HUD_ADAS']['SET_DISTANCE']
    self.is_cruise_latch = False
    self.prev_angle = 0
    self.lss_state = 0
    self.lss_alert = 0
    self.tsr = 0
    self.ahb = 0
    self.passthrough = 0
    self.lka_on = 0
    self.HMA = 0
    self.pt2 = 0
    self.pt3 = 0
    self.pt4 = 0
    self.pt5 = 0
    self.lkas_rdy_btn = False

  def update(self, cp):
    ret = car.CarState.new_message()

    self.tsr = cp.vl["LKAS_HUD_ADAS"]['TSR']
    self.lka_on = cp.vl["LKAS_HUD_ADAS"]['STEER_ACTIVE_ACTIVE_LOW']
    self.lkas_rdy_btn = cp.vl["PCM_BUTTONS"]['LKAS_ON_BTN']
    self.abh = cp.vl["LKAS_HUD_ADAS"]['SET_ME_XFF']
    self.passthrough = cp.vl["LKAS_HUD_ADAS"]['SET_ME_X5F']
    self.HMA = cp.vl["LKAS_HUD_ADAS"]['HMA']
    self.pt2 = cp.vl["LKAS_HUD_ADAS"]['PT2']
    self.pt3 = cp.vl["LKAS_HUD_ADAS"]['PT3']
    self.pt4 = cp.vl["LKAS_HUD_ADAS"]['PT4']
    self.pt5 = cp.vl["LKAS_HUD_ADAS"]['PT5']
    self.counter_pcm_buttons = cp.vl["PCM_BUTTONS"]['COUNTER']

    # EV irrelevant messages
    ret.brakeHoldActive = False

    ret.wheelSpeeds = self.get_wheel_speeds(
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_FL'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_FR'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_BL'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_BL'], # TODO: why would BR make the value wrong? Wheelspeed sensor prob?
    )
    ret.vEgoRaw = mean([ret.wheelSpeeds.rr, ret.wheelSpeeds.rl, ret.wheelSpeeds.fr, ret.wheelSpeeds.fl])

    # unfiltered speed from CAN sensors
    ret.vEgo, ret.aEgo = self.update_speed_kf(ret.vEgoRaw)
    ret.vEgoCluster = ret.vEgo
    ret.standstill = ret.vEgoRaw < 0.01

    # safety checks to engage
    can_gear = int(cp.vl["DRIVE_STATE"]['GEAR'])

    ret.doorOpen = any([cp.vl["METER_CLUSTER"]['BACK_LEFT_DOOR'],
                     cp.vl["METER_CLUSTER"]['FRONT_LEFT_DOOR'],
                     cp.vl["METER_CLUSTER"]['BACK_RIGHT_DOOR'],
                     cp.vl["METER_CLUSTER"]['FRONT_RIGHT_DOOR']])

    ret.seatbeltUnlatched = cp.vl["METER_CLUSTER"]['SEATBELT_DRIVER'] == 0
    ret.gearShifter = self.parse_gear_shifter(self.shifter_values.get(can_gear, None))

    disengage = ret.doorOpen or ret.seatbeltUnlatched or ret.brakeHoldActive
    if disengage:
      self.is_cruise_latch = False

    # gas pedal
    ret.gas = cp.vl["PEDAL"]['GAS_PEDAL']
    ret.gasPressed = ret.gas > 0.01

    # brake pedal
    ret.brake = cp.vl["PEDAL"]['BRAKE_PEDAL']
    ret.brakePressed = bool(cp.vl["DRIVE_STATE"]["BRAKE_PRESSED"]) or ret.brake > 0.01

    # steer
    ret.steeringAngleDeg = cp.vl["STEER_MODULE_2"]['STEER_ANGLE_2']
    steer_dir = 1 if (ret.steeringAngleDeg - self.prev_angle >= 0) else -1
    self.prev_angle = ret.steeringAngleDeg
    ret.steeringTorque = cp.vl["STEERING_TORQUE"]['MAIN_TORQUE']
    ret.steeringTorqueEps = cp.vl["STEER_MODULE_2"]['DRIVER_EPS_TORQUE'] * steer_dir
    ret.steeringPressed = bool(abs(ret.steeringTorqueEps) > 6)
    ret.steerWarning = False
    ret.steerError = False       # TODO

    # TODO: get the real value
    ret.stockAeb = False
    ret.stockFcw = False
    ret.cruiseState.available = any([cp.vl["ACC_HUD_ADAS"]["ACC_ON1"], cp.vl["ACC_HUD_ADAS"]["ACC_ON2"]])

    distance_val = int(cp.vl["ACC_HUD_ADAS"]['SET_DISTANCE'])
    ret.cruiseState.setDistance = self.parse_set_distance(self.set_distance_values.get(distance_val, None))

    # engage and disengage logic, do we still need this?
    if (cp.vl["PCM_BUTTONS"]["SET_BTN"] != 0 or cp.vl["PCM_BUTTONS"]["RES_BTN"] != 0) and not ret.brakePressed:
      self.is_cruise_latch = True

    # this can override the above engage disengage logic
    if bool(cp.vl["ACC_CMD"]["ACC_REQ_NOT_STANDSTILL"]):
      self.is_cruise_latch = True

    # byd speedCluster will follow wheelspeed if cruiseState is not available
    if ret.cruiseState.available:
      ret.cruiseState.speedCluster = max(int(cp.vl["ACC_HUD_ADAS"]['SET_SPEED']), 30) * CV.KPH_TO_MS
    else:
      ret.cruiseState.speedCluster = 0

    ret.cruiseState.speed = ret.cruiseState.speedCluster
    ret.cruiseState.standstill = bool(cp.vl["ACC_CMD"]["STANDSTILL_STATE"])
    ret.cruiseState.nonAdaptive = False

    stock_acc_on =  bool(cp.vl["ACC_CMD"]["ACC_CONTROLLABLE_AND_ON"])
    if not ret.cruiseState.available or ret.brakePressed or not stock_acc_on:
      self.is_cruise_latch = False

    ret.cruiseState.enabled = self.is_cruise_latch

    # button presses
    ret.leftBlinker = bool(cp.vl["STALKS"]["LEFT_BLINKER"])
    ret.rightBlinker = bool(cp.vl["STALKS"]["RIGHT_BLINKER"])
    ret.genericToggle = bool(cp.vl["STALKS"]["GENERIC_TOGGLE"])
    ret.espDisabled = False

    # blindspot sensors
    if self.CP.enableBsm:
      # used for lane change so its okay for the chime to work on both side.
      ret.leftBlindspot = bool(cp.vl["BSM"]["LEFT_APPROACH"])
      ret.rightBlindspot = bool(cp.vl["BSM"]["RIGHT_APPROACH"])

    self.lss_state = cp.vl["LKAS_HUD_ADAS"]["LSS_STATE"]
    self.lss_alert = cp.vl["LKAS_HUD_ADAS"]["SETTINGS"]
    return ret


  @staticmethod
  def get_can_parser(CP):
    signals = [
      # sig_name, sig_address, default
      ("GEAR", "DRIVE_STATE", 1),
      ("BRAKE_PRESSED", "DRIVE_STATE", 0),
      ("SET_DISTANCE", "ACC_HUD_ADAS", 1),
      ("SET_SPEED", "ACC_HUD_ADAS", 0),
      ("ACC_ON1", "ACC_HUD_ADAS", 0),
      ("ACC_ON2", "ACC_HUD_ADAS", 0),
      ("WHEELSPEED_FR", "WHEEL_SPEED", 0.),
      ("WHEELSPEED_FL", "WHEEL_SPEED", 0.),
      ("WHEELSPEED_BR", "WHEEL_SPEED", 0.),
      ("WHEELSPEED_BL", "WHEEL_SPEED", 0.),
      ("GAS_PEDAL", "PEDAL", 0.),
      ("BRAKE_PEDAL", "PEDAL", 0.),
      ("SEATBELT_DRIVER", "METER_CLUSTER", 0),
      ("BACK_RIGHT_DOOR", "METER_CLUSTER", 0),
      ("BACK_LEFT_DOOR", "METER_CLUSTER", 0),
      ("FRONT_RIGHT_DOOR", "METER_CLUSTER", 0),
      ("FRONT_LEFT_DOOR", "METER_CLUSTER", 0),
      ("STEER_ANGLE_2", "STEER_MODULE_2", 0.),
      ("DRIVER_EPS_TORQUE", "STEER_MODULE_2", 0.),
      ("MAIN_TORQUE", "STEERING_TORQUE", 0.),
      ("GENERIC_TOGGLE", "STALKS", 0),
      ("RIGHT_BLINKER", "STALKS", 0),
      ("LEFT_BLINKER", "STALKS", 0),
      ("LEFT_APPROACH", "BSM", 0),
      ("RIGHT_APPROACH", "BSM", 0),
      ("STANDSTILL_STATE", "ACC_CMD", 0),
      ("ACC_REQ_NOT_STANDSTILL", "ACC_CMD", 0),
      ("ACC_CONTROLLABLE_AND_ON", "ACC_CMD", 0),
      ("SET_BTN", "PCM_BUTTONS", 0),
      ("RES_BTN", "PCM_BUTTONS", 0),
      ("LSS_STATE", "LKAS_HUD_ADAS", 0),
      ("TSR", "LKAS_HUD_ADAS", 0),
      ("HMA", "LKAS_HUD_ADAS", 0),
      ("PT2", "LKAS_HUD_ADAS", 0),
      ("PT3", "LKAS_HUD_ADAS", 0),
      ("PT4", "LKAS_HUD_ADAS", 0),
      ("PT5", "LKAS_HUD_ADAS", 0),
      ("SET_ME_X5F", "LKAS_HUD_ADAS", 0),
      ("SET_ME_XFF", "LKAS_HUD_ADAS", 0),
      ("LKAS_ON_BTN", "PCM_BUTTONS", 0),
      ("COUNTER", "PCM_BUTTONS", 0),
      ("STEER_ACTIVE_ACTIVE_LOW", "LKAS_HUD_ADAS", 0),
      ("SETTINGS", "LKAS_HUD_ADAS", 0)
    ]
    checks = []

    # todo: make it such that enforce_checks=True
    return CANParser(DBC[CP.carFingerprint]['pt'], signals, checks, 0, enforce_checks=False)
